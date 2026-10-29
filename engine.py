"""Dependency-free, single-glycan contiguous peptide mass matching."""
import bisect
import math
import re
import heapq
from collections import Counter
from glycans import BY_ID, model_mass
from adducts import selected_adducts, hypotheses
from glycosmos import resolve_selections, selection_groups, image_ids

PROTON = 1.007276466621
from masses import (WATER, WATER_DECIMAL, MONO, FREE_AA, RESIDUE_DECIMAL,
                    FREE_DECIMAL, RESIDUE_UNITS, SCALE)
EXAMPLE = 'MRLAVGALLVCAVLGLCLAVPDKTVRWCAVSEHEATKCQSFRDHMKSVIPSDGPSVACVKKASYLDCIRAIAANEADAVTLDAGLVYDAYLAPNNLKPVVAEFYGSKEDPQTFYYAVAVVKKDSGFQMNQLRGKKSCHTGLGRSAGWNIPIGLLYCDLPEPRKPLEKAVANFFSGSCAPCADGTDFPQLCQLCPGCGCSTLNQYFGYSGAFKCLKDGAGDVAFVKHSTIFENLANKADRDQYELLCLDNTRKPVDEYKDCHLAQVPSHTVVARSMGGKEDLIWELLNQAQEHFGKDKSKEFQLFSSPHGKDLLFKDSAHGFLKVPPRMDAKMYLGYEYVTAIRNLREGTCPEAPTDECKPVKWCALSHHERLKCDEWSVNSVGKIECVSAETTEDCIAKIMNGEADAMSLDGGFVYIAGKCGLVPVLAENYNKSDNCEDTPEAGYFAIAVVKKSASDLTWDNLKGKKSCHTAVGRTAGWNIPMGLLYNKINHCRFDEFFSEGCAPGSKKDSSLCKLCMGSGLNLCEPNNKEGYYGYTGAFRCLVEKGDVAFVKHQTVPQNTGGKNPDPWAKNLNEKDYELLCLDGTRKPVEEYANCHLARAPNHAVVTRKDKEACVHKILRQQQHLFGSNVTDCSGNFCLFRSETKDLLFRDDTVCLAKLHDRNTYEKYLGEEYVKAVGNLRKCSTSSLLEACTFRRP'

def sequence(text):
    lines = str(text).strip().splitlines()
    if sum(line.startswith('>') for line in lines) > 1:
        raise ValueError('Paste one protein / FASTA record at a time.')
    seq = re.sub(r'\s+', '', ''.join(x for x in lines if not x.startswith('>'))).upper()
    if not seq or len(seq) > 10000:
        raise ValueError('Sequence must contain 1–10,000 amino acids.')
    invalid = sorted(set(seq) - set(MONO))
    if invalid:
        raise ValueError('Unsupported sequence characters: ' + ', '.join(invalid))
    return seq

def positive(value, name):
    try:
        value = float(value)
    except (TypeError, ValueError):
        raise ValueError(name + ' must be a number.') from None
    if not math.isfinite(value) or value <= 0:
        raise ValueError(name + ' must be finite and greater than zero.')
    return value

def integer(value, name, maximum):
    n = positive(value, name)
    if n != int(n) or n > maximum:
        raise ValueError(f'{name} must be an integer from 1 to {maximum}.')
    return int(n)

AA_NAMES = dict(zip('ARNDCEQGHILKMFPSTWYV', [
    'Alanine','Arginine','Asparagine','Aspartic acid','Cysteine','Glutamic acid',
    'Glutamine','Glycine','Histidine','Isoleucine','Leucine','Lysine','Methionine',
    'Phenylalanine','Proline','Serine','Threonine','Tryptophan','Tyrosine','Valine']))


def mass_breakdown(peptide, mode, glycan):
    """Use the same residue constants as the search; keep all corrections explicit."""
    table = RESIDUE_DECIMAL
    water_per_aa = WATER_DECIMAL*0
    amino_acids = []
    for aa,count in sorted(Counter(peptide).items()):
        amino_acids.append(dict(aa=aa,name=AA_NAMES[aa],count=count,
            input_mass=float(table[aa]),water_per_aa=float(water_per_aa),
            residue_mass=MONO[aa],subtotal=float(count*RESIDUE_DECIMAL[aa]),
            input_mass_text=format(table[aa], ".11f"),
            water_per_aa_text=format(water_per_aa, ".11f"),
            residue_mass_text=format(RESIDUE_DECIMAL[aa], ".11f"),
            subtotal_text=format(count*RESIDUE_DECIMAL[aa], ".11f")))
    input_total = float(sum(table[aa]*count for aa,count in Counter(peptide).items()))
    aa_total = sum(RESIDUE_UNITS[aa] for aa in peptide)/SCALE
    terminal_water = WATER
    attachment_water = WATER if mode in ('user','mono_free','free') else 0.0
    return dict(amino_acids=amino_acids,aa_input_total=input_total,
        aa_water_loss=float(len(peptide)*water_per_aa),aa_total=aa_total,
        terminal_water=terminal_water,attachment_water=attachment_water,
        peptide_mass=aa_total+terminal_water,glycan_mass=glycan,
        attached_glycan_mass=glycan-attachment_water,
        calculated_mass=aa_total+terminal_water+glycan-attachment_water)


def analyze(data):
    seq = sequence(data.get('sequence', ''))
    legacy_mode = data.get('mode', 'mono_free')
    if legacy_mode not in ('user','mono_free','mono_attached'):
        raise ValueError('Unknown mass model.')
    glycan_mass_type = data.get('glycan_mass_type', 'attached' if legacy_mode == 'mono_attached' else 'free')
    if glycan_mass_type not in ('free','attached'):
        raise ValueError('Glycan mass type must be free or attached.')
    mode = 'mono_attached' if glycan_mass_type == 'attached' else 'mono_free'
    ion_mode = data.get('ion_mode', 'positive')
    if ion_mode not in ('positive','negative'):
        raise ValueError('Ion mode must be positive or negative.')
    adduct_entries = selected_adducts(data, ion_mode)
    unit = data.get('unit','ppm')
    if unit not in ('ppm','Da','mz'):
        raise ValueError('Tolerance must be ppm, m/z, or neutral-mass Da.')
    tolerance = positive(data.get('tolerance',20), 'Tolerance')
    lo = integer(data.get('min_length',1), 'Minimum length',10000)
    hi = integer(data.get('max_length',60), 'Maximum length',10000)
    if hi < lo or lo > len(seq):
        raise ValueError('Check peptide length range against the sequence length.')
    hi = min(hi,len(seq))
    sites = []
    for line in str(data.get('sites','')).strip().splitlines():
        parts = [x.strip() for x in line.split(',')]
        if len(parts) != 2 or parts[1].upper() not in ('N','O'):
            raise ValueError('Sites require one position,type per line, e.g. 51,O.')
        pos = integer(parts[0], 'Site position',len(seq))
        kind = parts[1].upper()
        if seq[pos-1] not in ('N' if kind == 'N' else 'ST'):
            raise ValueError(f'Site {pos} is {seq[pos-1]}, incompatible with {kind}-linked glycosylation (N or S/T). Check precursor numbering.')
        sites.append((pos,kind))
    sites = sorted(set(sites))
    if not sites:
        raise ValueError('Enter at least one glycosylation site.')
    selected = data.get('selected_glycans',[])
    if not isinstance(selected,list) or any(not isinstance(x,str) or x not in BY_ID for x in selected):
        raise ValueError('Unknown glycan library selection. Reload the page and try again.')
    selected = sorted(set(selected))
    extended = resolve_selections(data.get('glycosmos_glycans', []))
    bulk, coverage = {}, None
    if 'glycosmos_selection' in data:
        bulk, coverage = selection_groups(data['glycosmos_selection'], data.get('glycosmos_apply_to', 'both'),
                                          [kind for _,kind in sites], mode)
    glycans = {}
    glycan_info = {}
    for kind in ('O','N'):
        values = re.split(r'[,\s]+', str(data.get('glycans_'+kind,'' )).strip())
        custom = set(positive(x,kind+' glycan mass') for x in values if x)
        if len(custom) > 20:
            raise ValueError('Use at most 20 custom glycan masses per type.')
        if glycan_mass_type == 'free' and any(value <= WATER for value in custom):
            raise ValueError('Free glycan mass must exceed the attachment water mass.')
        masses = set(custom)
        for value in custom:
            glycan_info[(kind,value)] = dict(ids=[],names=[],custom=True)
        for glycan in [BY_ID[identifier] for identifier in selected] + extended:
            identifier = glycan['id']
            if glycan['type'] != kind:
                continue
            value = model_mass(glycan,mode)
            masses.add(value)
            info = glycan_info.setdefault((kind,value),dict(ids=[],names=[],custom=False))
            info['ids'].append(identifier)
            info['names'].append(glycan['name'])
        glycans[kind] = sorted(masses)
        for value, group in bulk.get(kind, {}).items():
            info = glycan_info.setdefault((kind,value), dict(ids=[],names=[],custom=False))
            info['library_count'] = len(info['ids']) + group['library_count']
            info['ids'] = list(dict.fromkeys(info['ids'] + group['ids']))[:20]
            info['names'] = list(dict.fromkeys(info['names'] + group['names']))[:20]
            masses.add(value)
        glycans[kind] = sorted(masses)
    skipped_sites = [dict(site=pos, type=kind) for pos,kind in sites if not glycans[kind]]
    sites = [(pos,kind) for pos,kind in sites if glycans[kind]]
    if not sites:
        raise ValueError('No linkage-compatible glycans with supported masses are available for the selected protein sites. N-linked glycans require N sites; O-linked glycans require S/T sites. Select compatible glycans or enter a custom mass under its correct linkage type. GlyCosmos records without N/O annotations cannot be assigned to protein sites.')
    observations = []
    for line in str(data.get('observations','')).strip().splitlines():
        p = re.split(r'[,\s]+',line.strip())
        if len(p) != 2:
            raise ValueError('Measurements require m/z,charge per line (no header).')
        mz,z = positive(p[0],'m/z'), integer(p[1],'Charge',100)
        ions = hypotheses(adduct_entries, ion_mode, z)
        targets = [(ion, (z*mz-ion['mass_shift'])/ion['multimer']) for ion in ions]
        targets = [(ion, target) for ion, target in targets if target > 0]
        if adduct_entries is None and not targets:
            raise ValueError('Observed m/z and charge must give a positive neutral mass.')
        observations.append((mz,z,targets))
    if not 1 <= len(observations) <= 100:
        raise ValueError('Enter between 1 and 100 measurements.')
    estimate = sum(hi*(hi+1)//2 for _,k in sites)
    if estimate > 1500000:
        raise ValueError('Search is too large. Reduce maximum peptide length or site count (limit: 1.5 million peptide/site combinations).')
    # Exact integer prefix sums avoid floating-point accumulation and cancellation.
    prefix = [0]
    for aa in seq:
        prefix.append(prefix[-1]+RESIDUE_UNITS[aa])
    terminal = WATER if mode == 'mono_attached' else 0.0
    peptides = []
    for pos,kind in sites:
        for start in range(max(0,pos-hi),pos):
            for end in range(max(pos,start+lo),min(len(seq),start+hi)+1):
                base = (prefix[end]-prefix[start])/SCALE+terminal
                peptides.append((base,start,end,pos,kind))
    candidate_count = sum(len(glycans[kind]) for _,_,_,_,kind in peptides)
    if len(peptides)*sum(len(targets) for _,_,targets in observations) > 1500000:
        raise ValueError('Search is too large. Reduce peptide length, sites, measurements, or selected adducts.')
    results = []
    for mz,z,targets in observations:
        best, match_count = [], 0
        for ion, target in targets:
            multimer, shift = ion['multimer'], ion['mass_shift']
            window = mz*z*tolerance/1e6/multimer if unit == 'ppm' else tolerance*z/multimer if unit == 'mz' else tolerance
            for base,start,end,pos,kind in peptides:
                values = glycans[kind]
                glycan_target = target-base
                left = bisect.bisect_left(values,glycan_target-window)
                right = bisect.bisect_right(values,glycan_target+window)
                match_count += right-left
                middle = bisect.bisect_left(values,glycan_target)
                lower, upper = middle-1, middle
                # Merge both sides in increasing error, retaining only a
                # bounded global top 25 instead of materializing the product.
                for _ in range(min(25,len(values))):
                    if lower >= 0 and (upper >= len(values) or abs(values[lower]-glycan_target) <= abs(values[upper]-glycan_target)):
                        glycan = values[lower]; lower -= 1
                    else:
                        glycan = values[upper]; upper += 1
                    mass = base+glycan
                    error = abs((multimer*mass+shift)/z-mz)
                    if len(best) == 25 and error > -best[0][0]:
                        break
                    # Negative keys make heap[0] the worst retained result.
                    key = (-error,-mass,-start,-end,-pos,-ord(kind),-glycan,tuple(-ord(c) for c in ion['id']))
                    item = (key[0],key,(mass,start,end,pos,kind,glycan),ion,target,window)
                    if len(best) < 25:
                        heapq.heappush(best,item)
                    elif key > best[0][1]:
                        heapq.heapreplace(best,item)
        nearest = sorted(best,key=lambda item: (-item[0],item[2],item[3]['id']))
        rows = []
        for _,_,candidate,ion,target,window in nearest:
            mass,start,end,pos,kind,glycan = candidate
            predicted_mz = (ion['multimer']*mass+ion['mass_shift'])/z
            delta_mz = predicted_mz-mz
            delta_mass = mass-target
            rows.append(dict(sequence=seq[start:end],start=start+1,end=end,site=pos,type=kind,
                glycan=glycan,mass=mass,neutral_mass=mass,predicted_mz=predicted_mz,
                adduct=ion['name'],adduct_id=ion['id'],adduct_mass_shift=ion['mass_shift'],
                multimer=ion['multimer'],target_neutral_mass=target,
                signed_delta_mz=delta_mz,absolute_delta_mz=abs(delta_mz),
                error_da=delta_mass,error_ppm=delta_mz/mz*1e6,absolute_ppm_error=abs(delta_mz/mz*1e6),
                within=abs(delta_mass)<=window,
                glycan_info=glycan_info[(kind,glycan)],breakdown=mass_breakdown(seq[start:end],mode,glycan)))
        target = targets[0][1] if len(targets) == 1 else None
        results.append(dict(mz=mz,charge=z,target_mass=target,target_neutral_mass=target,
                            match_count=match_count,rows=rows,adduct_count=len(targets),
                            note='' if targets else 'No selected adduct has this charge and a positive neutral target. Select a compatible adduct or use proton-only mode.'))
    reference_ids = {identifier for result in results for row in result['rows']
                     for identifier in row['glycan_info']['ids']}
    return dict(sequence=seq,candidate_count=candidate_count,results=results,mode=mode,skipped_sites=skipped_sites,
                glycosmos_image_ids=sorted(reference_ids & image_ids()),
                ion_mode=ion_mode,glycan_mass_type=glycan_mass_type,
                selected_glycans=selected,proton_mass=PROTON,water_mass=WATER,
                selected_adducts=data.get('selected_adducts'),
                glycosmos_glycans=[dict(id=g['id'],type=g['type']) for g in extended],
                glycosmos_selection=data.get('glycosmos_selection'),glycosmos_coverage=coverage,
                aa_mass_basis="NIST isotope-derived monoisotopic masses",aa_decimal_places=11)
