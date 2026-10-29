"""Display source-derived monosaccharide names, independently of mass support."""
from collections import Counter
from pathlib import Path
import json
import re

DATA = Path(__file__).parent / 'data'
NAMES = json.loads((DATA / 'residue-names.json').read_text())['names']
UNITS = json.loads((DATA / 'residue-units.json').read_text())


def normalize(unit):
    # Anomer and ring size do not change the parent sugar's name. Preserve
    # stereochemistry, alditol ends and every chemical substituent.
    parts = unit.split('_')
    backbone = re.sub(r'-[12][abx]', '', parts[0])
    if backbone.startswith('u'):
        backbone = 'a' + backbone[1:]
    if backbone.startswith('AU'):
        backbone = 'Aa' + backbone[2:]
    return '_'.join([backbone] + [p for p in parts[1:] if not re.fullmatch(r'[12]-[1-9?]', p)])


PARENT_UNITS = {}
for unit, item in UNITS.items():
    # Ring-specific source abbreviations are left to exact matches.
    label = item['label']
    if label in NAMES:
        PARENT_UNITS.setdefault(normalize(unit), item)


def full_name(label):
    if label in NAMES:
        return NAMES[label]
    for symbol in sorted(NAMES, key=len, reverse=True):
        if len(symbol) < 3:
            continue
        if label.startswith(symbol) and label[len(symbol):len(symbol)+1] in ('f', '-', '2','3','4','5','6','7','8','9'):
            return f'{NAMES[symbol]} derivative ({label})'
    return label  # Retain the source notation when expansion is unknown.


def subunits(wurcs):
    match = re.fullmatch(r'WURCS=2\.0/[^/]+/((?:\[[^]]+\])+)/([^/]+)/.*', wurcs)
    if not match:
        return [], 'Subunit identities are not defined by this source record.'
    units = re.findall(r'\[([^]]+)\]', match.group(1))
    indices = re.findall(r'\d+', match.group(2))
    counts = Counter(int(i)-1 for i in indices)
    items = []
    for index, count in counts.items():
        if not 0 <= index < len(units):
            continue
        unit = units[index]
        named = UNITS.get(unit) or PARENT_UNITS.get(normalize(unit))
        if named:
            items.append(dict(symbol=named['label'], name=full_name(named['label']), count=count,
                              source_id=named['source_id'], wurcs=unit))
        else:
            # A generic formula class can be reported without inventing a
            # specific stereoisomer or dropping an unrecognized modification.
            from glycosmos import residue_composition
            try:
                composition = residue_composition(unit)
                aliases = dict(Deoxyhex='dHex', Pent='Pen', KDN='Kdn', NeuAc='Neu5Ac', NeuGc='Neu5Gc', Phos='P', Sulph='S')
                label = ' + '.join(aliases.get(s,s) for s in composition)
                name = ' + '.join(NAMES.get(aliases.get(s,s),s) for s in composition)
                items.append(dict(symbol=label, name=name + ' (specific identity unspecified)', count=count, wurcs=unit))
            except ValueError:
                items.append(dict(symbol='Unresolved residue', name='No supported residue name; WURCS notation retained', count=count, wurcs=unit))
    note = 'Counts describe the encoded repeat unit, not a fixed whole-glycan composition.' if '~' in wurcs or '<' in match.group(2) else ''
    return items, note


def details(record):
    items, note = subunits(record['wurcs'])
    return dict(name=record['iupac'] or 'Glycan ' + record['id'],
                name_kind='IUPAC condensed' if record['iupac'] else 'Accession label; no IUPAC name supplied',
                subunits=items, subunit_note=note,
                image_url='/api/glycosmos/image/' + record['id'])
