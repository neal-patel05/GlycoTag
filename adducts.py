"""Fiehn Lab's complete calculator table (31 positive and 15 negative ions).

Use the corrected spreadsheet's numeric column, including M+H+2Na and
3M-H corrections. Other offsets retain the source's published precision;
pure protonation/deprotonation uses the engine proton constant. Solvent
cluster labels, including the source's unusual IsoProp+Na+H, are retained
as tabulated rather than silently assigning a different species.
"""
import json
from pathlib import Path

PROTON = 1.007276466621
SOURCE = json.loads((Path(__file__).parent / 'data/adducts.json').read_text())
CATALOG = SOURCE['adducts']
BY_ID = {entry['id']: entry for entry in CATALOG}


def catalog_payload():
    return dict(SOURCE, note='Offsets follow the corrected Fiehn Lab spreadsheet. '
                'Pure proton ions use 1.007276466621 Da. ACN = acetonitrile; '
                'IsoProp = isopropanol; FA = formic acid; Hac = acetic acid; '
                'TFA = trifluoroacetic acid; DMSO = dimethyl sulfoxide. '
                'The IsoProp+Na+H label and offset are retained as tabulated.')


def selected_adducts(data, mode):
    selected = data.get('selected_adducts')
    if selected is None:  # Backwards-compatible, arbitrary-charge proton ions.
        return None
    if not isinstance(selected, list) or any(not isinstance(x, str) or x not in BY_ID for x in selected):
        raise ValueError('Unknown adduct selection. Reload the page and try again.')
    entries = [BY_ID[x] for x in dict.fromkeys(selected) if BY_ID[x]['mode'] == mode]
    if not entries:
        raise ValueError(f'Select at least one {mode} ion adduct, or use proton-only mode.')
    return entries


def hypotheses(entries, mode, charge):
    sign = 1 if mode == 'positive' else -1
    if entries is None:
        return [dict(id='proton', name=f'[M{ "+" if sign == 1 else "-"}{charge}H]{charge}{ "+" if sign == 1 else "-"}',
                     mode=mode, charge=charge, multimer=1, mass_shift=sign*charge*PROTON)]
    result = []
    for entry in entries:
        if entry['charge'] != charge:
            continue
        shift = entry['mz_offset'] * charge
        if entry['name'] in ('[M+H]+', '[M+2H]2+', '[M+3H]3+',
                             '[M-H]-', '[M-2H]2-', '[M-3H]3-',
                             '[2M+H]+', '[2M-H]-', '[3M-H]-'):
            shift = sign*charge*PROTON
        result.append(dict(entry, mass_shift=shift))
    return result
