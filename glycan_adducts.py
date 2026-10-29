"""One best glycan-only or adduct-associated assignment per observed ion."""
import bisect
import math
import re

from adducts import selected_adducts, hypotheses
from engine import positive, integer
from glycans import BY_ID, model_mass
from glycosmos import resolve_selections, selection_groups, image_ids
from masses import WATER


def analyze_glycan_adducts(data):
    mode = data.get('ion_mode', 'positive')
    if mode not in ('positive', 'negative'):
        raise ValueError('Ion mode must be positive or negative.')
    mass_type = data.get('glycan_mass_type', 'free')
    if mass_type not in ('free', 'attached'):
        raise ValueError('Glycan mass type must be free or attached.')
    unit = data.get('unit', 'ppm')
    if unit not in ('ppm', 'Da', 'mz'):
        raise ValueError('Tolerance must be ppm, m/z, or neutral-mass Da.')
    tolerance = positive(data.get('tolerance', 20), 'Tolerance')
    entries = selected_adducts(data, mode)
    groups = {}

    def add(mass, ids=(), names=(), custom=False):
        group = groups.setdefault(mass, dict(ids=[], names=[], custom=False))
        group['ids'] = list(dict.fromkeys(group['ids'] + list(ids)))[:20]
        group['names'] = list(dict.fromkeys(group['names'] + list(names)))[:20]
        group['custom'] = group['custom'] or custom

    selected = data.get('selected_glycans', [])
    if not isinstance(selected, list) or any(not isinstance(i, str) or i not in BY_ID for i in selected):
        raise ValueError('Unknown glycan library selection. Reload the page and try again.')
    for g in [BY_ID[i] for i in dict.fromkeys(selected)] + resolve_selections(data.get('glycosmos_glycans', []), require_linkage=False):
        add(model_mass(g, 'mono_free'), [g['id']], [g['name']])
    coverage = None
    if 'glycosmos_selection' in data:
        # Free glycan signals have no protein attachment site. Keep every
        # selected linkage here, including records without N/O annotations.
        bulk, coverage = selection_groups(data['glycosmos_selection'], 'both', ['N', 'O'], 'mono_free', require_linkage=False)
        for mass, group in bulk.get('N', {}).items():
            add(mass, group['ids'], group['names'])
    for kind in ('N', 'O'):
        values = {positive(v, kind + ' glycan mass') for v in re.split(r'[,\s]+', str(data.get('glycans_' + kind, '')).strip()) if v}
        if len(values) > 20:
            raise ValueError('Use at most 20 custom glycan masses per type.')
        for mass in values:
            if mass_type == 'free' and mass <= WATER:
                raise ValueError('Free glycan mass must exceed the attachment water mass.')
            add(mass + (WATER if mass_type == 'attached' else 0), custom=True)
    if not groups:
        raise ValueError('Select at least one glycan with a supported mass, or enter a custom glycan mass.')
    masses = sorted(groups)
    lines = str(data.get('observations', '')).strip().splitlines()
    if not 1 <= len(lines) <= 100:
        raise ValueError('Enter between 1 and 100 measurements.')
    results, seen = [], set()
    for line in lines:
        parts = re.split(r'[,\s]+', line.strip())
        if len(parts) != 2:
            raise ValueError('Measurements require m/z,charge per line (no header).')
        mz, z = positive(parts[0], 'm/z'), integer(parts[1], 'Charge', 100)
        if (mz, z) in seen:
            continue
        seen.add((mz, z))
        calculated = mz * z - z
        if not math.isfinite(calculated) or calculated <= 0:
            raise ValueError('The requested formula (m/z × z) − z must give a positive finite MW.')
        # A glycan-only signal still requires charge carriers. Monomeric
        # protonation/deprotonation is the baseline, not an extra adduct.
        baseline = hypotheses(None, mode, z)[0]
        ions = [(baseline, 'glycan_only')]
        for ion in hypotheses(entries, mode, z):
            if ion['multimer'] == 1 and math.isclose(ion['mass_shift'], baseline['mass_shift'], rel_tol=0, abs_tol=1e-9):
                continue  # Do not report the same proton-only explanation twice.
            ions.append((ion, 'adduct_associated'))
        assignments = []
        for ion, assignment_type in ions:
            n, shift = ion['multimer'], ion['mass_shift']
            target = (mz * z - shift) / n
            pivot = bisect.bisect_left(masses, target)
            # Both neighbors suffice to find the global minimum for each
            # monotonic ion model, including ties across different glycans.
            for mass in masses[max(0, pivot - 1):pivot + 1]:
                difference = mass - calculated
                expected = z - shift - (n - 1) * mass
                residual = difference - expected
                predicted = (n * mass + shift) / z
                if not math.isfinite(predicted) or predicted <= 0:
                    continue
                error_mz = residual / z
                error_ppm = error_mz / mz * 1e6
                error_da = residual / n
                error = {'ppm':abs(error_ppm), 'mz':abs(error_mz), 'Da':abs(error_da)}[unit]
                assignments.append(dict(assignment_type=assignment_type,
                    glycan_mw=mass, mass_difference=difference, glycan_info=groups[mass],
                    adduct_id=ion['id'] if assignment_type == 'adduct_associated' else None,
                    adduct=ion['name'] if assignment_type == 'adduct_associated' else None,
                    ion=ion['name'], multimer=n, mass_shift=shift,
                    expected_difference=expected, difference_error=residual,
                    predicted_mz=predicted, error_mz=error_mz, error_ppm=error_ppm,
                    error_da=error_da, within=error <= tolerance))
        eligible = [a for a in assignments if a['within']] or assignments
        smallest = min(abs(a['error_mz']) for a in eligible)
        tied = [a for a in eligible if math.isclose(abs(a['error_mz']), smallest, rel_tol=0, abs_tol=1e-9)]
        # Deterministic representative for numerical ties; explicitly flag the
        # ambiguity rather than claiming a uniquely identified ion.
        best = min(tied, key=lambda a: (a['assignment_type'] != 'glycan_only',
                                       a['glycan_mw'], a['adduct_id'] or ''))
        results.append(dict(mz=mz, charge=z, calculated_mw=calculated,
                            best_assignment=best, tied_assignment_count=len(tied),
                            assignments_evaluated=len(masses)*len(ions)))
    ids = {i for result in results for i in result['best_assignment']['glycan_info']['ids']}
    return dict(results=results, ion_mode=mode, tolerance=tolerance, unit=unit,
                formula='(m/z × z) − z', difference_formula='glycan MW − calculated MW',
                glycan_mass_basis='free neutral glycan', unique_masses=len(masses),
                glycosmos_coverage=coverage, glycosmos_image_ids=sorted(ids & image_ids()),
                selected_adducts=data.get('selected_adducts'),
                glycosmos_selection=data.get('glycosmos_selection'))
