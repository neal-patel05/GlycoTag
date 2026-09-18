"""Dependency-free, single-glycan contiguous peptide mass matching."""
import bisect
import math
import re
from collections import Counter
from glycans import BY_ID, model_mass

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
    table = FREE_DECIMAL if mode == 'user' else RESIDUE_DECIMAL
    water_per_aa = WATER_DECIMAL if mode == 'user' else WATER_DECIMAL*0
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
    terminal_water = WATER if mode != 'user' else 0.0
    attachment_water = WATER if mode == 'mono_free' else 0.0
    return dict(amino_acids=amino_acids,aa_input_total=input_total,
        aa_water_loss=float(len(peptide)*water_per_aa),aa_total=aa_total,
        terminal_water=terminal_water,attachment_water=attachment_water,
        peptide_mass=aa_total+terminal_water,glycan_mass=glycan,
        calculated_mass=aa_total+terminal_water+glycan-attachment_water)


def analyze(data):
    seq = sequence(data.get('sequence', ''))
    mode = data.get('mode', 'user')
    if mode not in ('user','mono_free','mono_attached'):
        raise ValueError('Unknown mass model.')
    unit = data.get('unit','ppm')
    if unit not in ('ppm','Da'):
        raise ValueError('Tolerance must be ppm or neutral-mass Da.')
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
    glycans = {}
    glycan_info = {}
    for kind in ('O','N'):
        values = re.split(r'[,\s]+', str(data.get('glycans_'+kind,'' )).strip())
        custom = set(positive(x,kind+' glycan mass') for x in values if x)
        if len(custom) > 20:
            raise ValueError('Use at most 20 custom glycan masses per type.')
        masses = set(custom)
        for value in custom:
            glycan_info[(kind,value)] = dict(ids=[],names=[],custom=True)
        for identifier in selected:
            glycan = BY_ID[identifier]
            if glycan['type'] != kind:
                continue
            value = model_mass(glycan,mode)
            masses.add(value)
            info = glycan_info.setdefault((kind,value),dict(ids=[],names=[],custom=False))
            info['ids'].append(identifier)
            info['names'].append(glycan['name'])
        glycans[kind] = sorted(masses)
    if any(not glycans[k] for _,k in sites):
        raise ValueError('Enter a glycan mass for each selected site type.')
    observations = []
    for line in str(data.get('observations','')).strip().splitlines():
        p = re.split(r'[,\s]+',line.strip())
        if len(p) != 2:
            raise ValueError('Measurements require m/z,charge per line (no header).')
        mz,z = positive(p[0],'m/z'), integer(p[1],'Charge',100)
        target = z*mz
        observations.append((mz,z,target))
    if not 1 <= len(observations) <= 100:
        raise ValueError('Enter between 1 and 100 measurements.')
    estimate = sum(hi*(hi+1)//2 * len(glycans[k]) for _,k in sites)
    if estimate > 1500000:
        raise ValueError('Search is too large. Reduce maximum peptide length, site count, or glycan count (limit: 1.5 million candidates).')
    # Exact integer prefix sums avoid floating-point accumulation and cancellation.
    prefix = [0]
    for aa in seq:
        prefix.append(prefix[-1]+RESIDUE_UNITS[aa])
    terminal = WATER if mode == 'mono_attached' else 0.0
    candidates = []
    for pos,kind in sites:
        for start in range(max(0,pos-hi),pos):
            for end in range(max(pos,start+lo),min(len(seq),start+hi)+1):
                for glycan in glycans[kind]:
                    mass = (prefix[end]-prefix[start])/SCALE+terminal+glycan
                    candidates.append((mass,start,end,pos,kind,glycan))
    candidates.sort()
    masses = [c[0] for c in candidates]
    results = []
    for mz,z,target in observations:
        window = target*tolerance/1e6 if unit == 'ppm' else tolerance
        # Candidates are stored as raw peptide + glycan masses. Translate the
        # lookup center to that space; the input target itself stays unchanged.
        center = target+z*PROTON
        left = bisect.bisect_left(masses,center-window)
        right = bisect.bisect_right(masses,center+window)
        middle = bisect.bisect_left(masses,center)
        nearest = sorted(candidates[max(0,middle-25):middle+25],key=lambda c:(abs(c[0]-center),c))[:25]
        rows = []
        for mass,start,end,pos,kind,glycan in nearest:
            corrected_mass = mass-z*PROTON
            delta = corrected_mass-target
            rows.append(dict(sequence=seq[start:end],start=start+1,end=end,site=pos,type=kind,glycan=glycan,mass=mass,corrected_mass=corrected_mass,predicted_mz=mass/z-PROTON,error_da=delta,error_ppm=delta/target*1e6,within=abs(delta)<=window,
                             glycan_info=glycan_info[(kind,glycan)],breakdown=mass_breakdown(seq[start:end],mode,glycan)))
        results.append(dict(mz=mz,charge=z,target_mass=target,match_count=right-left,rows=rows))
    return dict(sequence=seq,candidate_count=len(candidates),results=results,mode=mode,
                selected_glycans=selected,proton_mass=PROTON,water_mass=WATER,
                aa_mass_basis="NIST isotope-derived monoisotopic masses",aa_decimal_places=11)
