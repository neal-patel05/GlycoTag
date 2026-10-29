"""Local, paginated GlyCosmos catalogue and conservative ExPASy composition masses.

WURCS reference: https://github.com/glycoinfo/WURCS/wiki/Representations-for-Basic-Structures
Only fixed compositions representable by GlycanMass are assigned a mass.
Unknown chemistry, repeats, cyclic glycans and reduced ends are never guessed.
"""
from collections import Counter
from contextlib import closing
from decimal import Decimal
from pathlib import Path
import json
import gzip
import hashlib
import math
import os
import re
import shutil
import sqlite3
import tempfile
import threading
from functools import lru_cache

from glycans import EXPASY_REDUCING_END
from masses import WATER_DECIMAL

DATABASE = Path(__file__).parent / 'data/glycosmos.sqlite3'
SNAPSHOT = DATABASE.with_suffix('.sqlite3.gz')
IMAGE_DATABASE = DATABASE.with_name('glycosmos-images.sqlite3')
IMAGE_MANIFEST = DATABASE.with_name('glycosmos-image-ids.json.gz')
IMAGE_ORIGIN = 'https://image.glycosmos.org/snfg/png/'
_cache_path = None
_cache_lock = threading.Lock()
RESIDUES = {k: Decimal(v) for k, v in dict(
    Hex='162.0528', HexNAc='203.0794', Deoxyhex='146.0579', Pent='132.0423',
    HexA='176.03209', KDN='250.0689', NeuAc='291.0954', NeuGc='307.0903',
    Phos='79.9663', Sulph='79.9568').items()}

# GlyCosmos CSV's Motif name(s) annotations, matched by motif accession.
# Attachment categories: https://glycosmos.org/glyconavi/motifs
# Do not classify by mass, monosaccharide composition, or N/O protein sites.
N_MOTIFS = frozenset(('G00026MO', 'G00028MO', 'G00029MO', 'G00030MO'))
O_MOTIFS = frozenset(f'G{i:05}MO' for i in range(31, 45)) | frozenset((
    'G82576YO',  # O-Fucose
    'G70323CJ', 'G59126YU', 'G81516LH',  # O-Mannose variants
    'G06953DR', 'G06898ZN',  # branched O-mannose / HNK-1 on O-mannose
    'G01614ZM', 'G65562ZE',  # di-sialyl / sialyl T antigen
    'G57321FI', 'G36123IU',  # Tn / sialyl Tn antigen
))


def linkage_types(motifs):
    accessions = {item.rsplit(':', 1)[-1].strip() for item in motifs.split(',')}
    return [kind for kind, known in (('N', N_MOTIFS), ('O', O_MOTIFS))
            if accessions & known]


@lru_cache(maxsize=1)
def linkage_index():
    groups = {'all': set(), 'N': set(), 'O': set()}
    with closing(connect()) as db:
        for identifier, motifs in db.execute('SELECT id,motifs FROM glycans'):
            groups['all'].add(identifier)
            for kind in linkage_types(motifs):
                groups[kind].add(identifier)
    groups['both'] = groups['N'] | groups['O']
    return {key: frozenset(value) for key, value in groups.items()}


def linkage_ids(kind):
    if not isinstance(kind, str) or kind not in ('all', 'N', 'O', 'both'):
        raise ValueError('Choose all, N, O, or both GlyCosmos linkage classes.')
    return linkage_index()[kind]


# Common search words resolve to the same source taxonomy IDs as scientific
# names. They never infer extra glycan/organism relationships.
SPECIES_ALIASES = {'9606': ('human', 'humans'),
                   '9913': ('bovine', 'cow', 'cows', 'cattle')}


def normalize_species(query):
    if not isinstance(query, str) or len(query) > 200:
        raise ValueError('Species must be a name or taxonomy ID of at most 200 characters.')
    return ' '.join(query.split()).casefold()


@lru_cache(maxsize=1)
def taxonomy_index():
    with closing(connect()) as db:
        taxa = {r['id']: dict(id=r['id'], name=r['name'], ids=set()) for r in db.execute('SELECT * FROM taxa')}
        for glycan_id, taxon_id in db.execute('SELECT glycan_id,taxon_id FROM glycan_taxa'):
            taxa[taxon_id]['ids'].add(glycan_id)
    for taxon in taxa.values():
        taxon['ids'] = frozenset(taxon['ids'])
    return taxa


def matching_species(query):
    query = normalize_species(query)
    taxa = list(taxonomy_index().values())
    if not query:
        return []
    exact = [t for t in taxa if query == t['id'] or query == t['name'].casefold()
             or query in SPECIES_ALIASES.get(t['id'], ())]
    if exact:
        return exact
    if query.isdigit():
        return []  # Taxonomy IDs are exact, never prefix matches.
    return sorted([t for t in taxa if query in t['name'].casefold()
                   or any(query in alias for alias in SPECIES_ALIASES.get(t['id'], ()))],
                  key=lambda t: (-len(t['ids']), t['name']))


@lru_cache(maxsize=32)
def species_ids(query):
    return frozenset(i for taxon in matching_species(query) for i in taxon['ids'])


def species_label(taxon):
    return dict(id=taxon['id'], name=taxon['name'], count=len(taxon['ids']),
                aliases=list(SPECIES_ALIASES.get(taxon['id'], ())))


def residue_composition(unit):
    parts = unit.split('_')
    skeleton = parts[0].split('-')[0]
    kind = None
    for pattern, name in [(r'[auo][12x]{4}h', 'Hex'), (r'ha[12x]{3}h', 'Hex'),
                          (r'[auo][12x]{4}m', 'Deoxyhex'), (r'[auo][12x]{3}h', 'Pent'),
                          (r'[auo][12x]{4}A', 'HexA'), (r'A[aU]d[12x]{5}h', 'KDN')]:
        if re.fullmatch(pattern, skeleton):
            kind = name
            break
    if kind is None:
        raise ValueError('Sugar chemistry or reducing end is outside the supported ExPASy inputs.')
    extra = Counter()
    for mod in parts[1:]:
        if re.fullmatch(r'[12]-[1-9?]', mod):  # ring closure only
            continue
        if re.fullmatch(r'2\*NCC/3=O', mod) and kind == 'Hex':
            kind = 'HexNAc'
        elif mod == '5*NCC/3=O' and kind == 'KDN':
            kind = 'NeuAc'
        elif mod == '5*NCCO/3=O' and kind == 'KDN':
            kind = 'NeuGc'
        elif re.fullmatch(r'[1-9?]\*OPO/3O/3=O', mod):
            extra['Phos'] += 1
        elif re.fullmatch(r'[1-9?]\*OSO/3=O/3=O', mod):
            extra['Sulph'] += 1
        else:
            raise ValueError('Substituent is outside the supported ExPASy inputs.')
    extra[kind] += 1
    return extra


def composition_mass(wurcs, reported_mass=None):
    match = re.fullmatch(r'WURCS=2\.0/(\d+),(\d+),(\d+)(\+?)/((?:\[[^]]+\])+)/([\d-]+)/([^~]*)', wurcs)
    if not match:
        raise ValueError('Missing or variable composition / repeating structure; no fixed ExPASy mass.')
    unique, count, links = map(int, match.group(1, 2, 3))
    if count < 1 or links >= count or '*' in match.group(7):
        raise ValueError('Cyclic or modified inter-residue bonds are not supported by this mass model.')
    units = re.findall(r'\[([^]]+)\]', match.group(5))
    indices = [int(i) for i in match.group(6).split('-')]
    if len(units) != unique or len(indices) != count or any(i < 1 or i > unique for i in indices):
        raise ValueError('Incomplete composition; no fixed ExPASy mass.')
    composition = Counter()
    unit_compositions = [residue_composition(u) for u in units]
    for index in indices:
        composition.update(unit_compositions[index-1])
    free = sum(RESIDUES[s]*n for s, n in composition.items()) + EXPASY_REDUCING_END
    # An independent, coarsely rounded source mass guards against chemistry
    # outside this intentionally restricted WURCS-to-composition mapping.
    if reported_mass and abs(float(free)-reported_mass) > max(.02, count*.0001):
        raise ValueError('Composition mass disagrees with the source; review this structure manually.')
    return dict(composition=composition, free_mass=float(free), attached_mass=float(free-WATER_DECIMAL))


def connect():
    global _cache_path
    with _cache_lock:
        if _cache_path is None:
            if not SNAPSHOT.exists():
                raise ValueError('GlyCosmos library is missing. Run scripts/build_glycosmos.py with the source CSV.')
            digest = hashlib.sha256(SNAPSHOT.read_bytes()).hexdigest()[:20]
            cache = Path(tempfile.gettempdir()) / ('glypeptid-' + digest + '.sqlite3')
            if not cache.exists():
                descriptor, name = tempfile.mkstemp(prefix='glypeptid-', suffix='.building')
                try:
                    with os.fdopen(descriptor, 'wb') as target, gzip.open(SNAPSHOT, 'rb') as source:
                        shutil.copyfileobj(source, target)
                    os.replace(name, cache)
                finally:
                    Path(name).unlink(missing_ok=True)
            _cache_path = cache
    db = sqlite3.connect(_cache_path.as_uri() + '?mode=ro', uri=True)
    db.row_factory = sqlite3.Row
    return db


def serialize(row):
    result = dict(row)
    result['composition'] = json.loads(result['composition'])
    result['source'] = 'https://glycosmos.org/glycans/' + result['id']
    from glycan_details import details
    result.update(details(result))
    result['has_image'] = result['id'] in image_ids()
    result['linkage_types'] = linkage_types(result['motifs'])
    result['organisms'] = [dict(name=entry.rpartition(':')[0], id=entry.rpartition(':')[2])
                           for entry in result['taxonomy'].split(',') if entry]
    return result


def record(identifier):
    if not re.fullmatch(r'G\d{5}[A-Z]{2}', identifier):
        raise ValueError('Invalid GlyTouCan accession.')
    with closing(connect()) as db:
        row = db.execute('SELECT * FROM glycans WHERE id=?', (identifier,)).fetchone()
    return serialize(row) if row is not None else None


@lru_cache(maxsize=1)
def image_ids():
    path = IMAGE_DATABASE
    if not path.exists():
        if IMAGE_MANIFEST.exists():
            with gzip.open(IMAGE_MANIFEST, 'rt', encoding='utf-8') as stream:
                return frozenset(json.load(stream))
        return frozenset()
    with closing(sqlite3.connect(path.as_uri() + '?mode=ro', uri=True)) as db:
        return frozenset(row[0] for row in db.execute('SELECT id FROM images'))


def structure_image(identifier):
    if not re.fullmatch(r'G\d{5}[A-Z]{2}', identifier):
        raise ValueError('Invalid GlyTouCan accession.')
    path = IMAGE_DATABASE
    if not path.exists():
        return None
    with closing(sqlite3.connect(path.as_uri() + '?mode=ro', uri=True)) as db:
        row = db.execute('SELECT png FROM images WHERE id=?', (identifier,)).fetchone()
    return row[0] if row else None


def hosted_image_url(identifier):
    """Use official images when the optional local archive is not installed."""
    if not re.fullmatch(r'G\d{5}[A-Z]{2}', identifier):
        raise ValueError('Invalid GlyTouCan accession.')
    return IMAGE_ORIGIN + identifier if not IMAGE_DATABASE.exists() and identifier in image_ids() else None


@lru_cache(maxsize=1)
def mass_index():
    with closing(connect()) as db:
        return {r['id']: (r['free_mass'], r['attached_mass'])
                for r in db.execute('SELECT id,free_mass,attached_mass FROM glycans ORDER BY id')}


def selection_groups(spec, apply_to, site_types, mode, *, require_linkage=True):
    """Group peptide candidates by annotated linkage, or free glycans without sites."""
    if not isinstance(spec, dict) or not isinstance(spec.get('all', False), bool):
        raise ValueError('Invalid GlyCosmos selection.')
    if apply_to not in ('both', 'N', 'O'):
        raise ValueError('Choose N, O, or both compatible protein site types.')
    index = mass_index()
    ids, excluded = spec.get('ids', []), spec.get('excluded', [])
    for values in (ids, excluded):
        if not isinstance(values, list) or any(not isinstance(i, str) or i not in index for i in values):
            raise ValueError('Unknown GlyCosmos accession in selection.')
    linkage = spec.get('linkage', 'all')
    allowed = linkage_ids(linkage)
    species = normalize_species(spec.get('species', ''))
    if species:
        allowed = allowed & species_ids(species)
    chosen = (set(allowed) if spec.get('all') else set(ids) & allowed) - set(excluded)
    kinds = set(site_types) if apply_to == 'both' else {apply_to} & set(site_types)
    groups = {kind: {} for kind in kinds}
    supported = matched = unclassified = incompatible = 0
    matched_masses = set()
    annotations = linkage_index()
    for identifier in sorted(chosen):
        free, attached = index[identifier]
        if free is None:
            continue
        supported += 1
        annotated_kinds = {kind for kind in ('N', 'O') if identifier in annotations[kind]}
        compatible = kinds & annotated_kinds if require_linkage else kinds
        if not compatible:
            if not annotated_kinds:
                unclassified += 1
            else:
                incompatible += 1
            continue
        matched += 1
        mass = attached if mode == 'mono_attached' else free
        matched_masses.add(mass)
        for kind in compatible:
            group = groups[kind].setdefault(mass, dict(ids=[], names=[], custom=False, library_count=0))
            group['library_count'] += 1
            if len(group['ids']) < 20:
                group['ids'].append(identifier)
                group['names'].append(identifier + ' (GlyCosmos)')
    return groups, dict(selected=len(chosen), supported=supported,
            skipped=len(chosen)-supported, unique_masses=len(matched_masses), apply_to=apply_to,
            matched=matched, linkage_skipped=supported-matched, unclassified=unclassified,
            incompatible=incompatible, linkage_enforced=require_linkage,
            linkage=linkage, species=species,
            species_matches=[species_label(t) for t in matching_species(species)])


def browse(params):
    query = params.get('q', [''])[0].strip()[:200]
    linkage = params.get('linkage', ['all'])[0]
    allowed = linkage_ids(linkage)
    species = normalize_species(params.get('species', [''])[0])
    organism_ids = species_ids(species) if species else None
    if organism_ids is not None:
        allowed = allowed & organism_ids
    try:
        page = int(params.get('page', ['1'])[0])
        if page < 1:
            raise ValueError()
    except (ValueError, TypeError):
        raise ValueError('Page must be a positive integer.') from None
    where, args = [], []
    if query:
        where.append('(id LIKE ? ESCAPE \'\\\' OR iupac LIKE ? ESCAPE \'\\\' OR motifs LIKE ? ESCAPE \'\\\' OR composition LIKE ? ESCAPE \'\\\')')
        term = '%' + query.replace('\\', '\\\\').replace('%', '\\%').replace('_', '\\_') + '%'
        args.extend([term]*4)
    if params.get('supported', ['0'])[0] == '1':
        where.append('free_mass IS NOT NULL')
    if linkage != 'all' or species:
        where.append('has_linkage(id)')
    clause = ' WHERE ' + ' AND '.join(where) if where else ''
    with closing(connect()) as db:
        db.create_function('has_linkage', 1, allowed.__contains__)
        total = db.execute('SELECT COUNT(*) FROM glycans' + clause, args).fetchone()[0]
        rows = db.execute('SELECT * FROM glycans' + clause + " ORDER BY (iupac != '') DESC, id LIMIT 30 OFFSET ?", args + [(page-1)*30]).fetchall()
        metadata = json.loads(db.execute('SELECT value FROM metadata').fetchone()[0])
    matches = matching_species(species)
    options = matches if species else sorted(taxonomy_index().values(), key=lambda t: (-len(t['ids']), t['name']))
    return dict(glycans=[serialize(r) for r in rows], total=total, page=page,
                pages=math.ceil(total/30), metadata=metadata, image_count=len(image_ids()),
                linkage=linkage, species=species,
                species_matches=[species_label(t) for t in matches],
                species_options=[species_label(t) for t in options[:25]],
                linkage_counts={key: len(ids & organism_ids if organism_ids is not None else ids)
                                for key, ids in linkage_index().items()})


def resolve_selections(selections, *, require_linkage=True):
    if not isinstance(selections, list) or len(selections) > 500:
        raise ValueError('Select at most 500 GlyCosmos entries per search.')
    entries, seen = [], set()
    if not selections:
        return entries
    with closing(connect()) as db:
        for selection in selections:
            if not isinstance(selection, dict) or selection.get('type') not in ('N', 'O') or not isinstance(selection.get('id'), str):
                raise ValueError('Each GlyCosmos selection needs an accession and an N/O site type.')
            identifier, kind = selection['id'], selection['type']
            if (identifier, kind) in seen:
                continue
            seen.add((identifier, kind))
            row = db.execute('SELECT * FROM glycans WHERE id=?', (identifier,)).fetchone()
            if row is None or row['free_mass'] is None:
                raise ValueError(f'{identifier}: no supported ExPASy mass is available.')
            g = serialize(row)
            if require_linkage and kind not in g['linkage_types']:
                annotated = '/'.join(g['linkage_types']) + '-linked' if g['linkage_types'] else 'no N/O linkage annotation'
                raise ValueError(f'{identifier}: {annotated}; cannot assign it to an {kind}-linked protein site.')
            g.update(type=kind, name=identifier + ' (GlyCosmos)')
            entries.append(g)
    return entries
