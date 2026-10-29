import unittest

from engine import analyze, MONO, PROTON
from glycan_adducts import analyze_glycan_adducts
from glycosmos import record, linkage_index, selection_groups
from masses import WATER
from test_expanded_library import request
import test_vercel


class AttachmentLinkageTests(unittest.TestCase):
    def test_wrong_linkage_cannot_be_assigned_even_at_exact_mass(self):
        for identifier, residue, site in [('G00026MO','S','O'), ('G00039MO','N','N')]:
            mz = MONO[residue] + record(identifier)['free_mass'] + PROTON
            for spec in ({'glycosmos_selection':{'ids':[identifier]}},
                         {'glycosmos_glycans':[{'id':identifier,'type':site}]}):
                with self.subTest(identifier=identifier, spec=spec), self.assertRaisesRegex(ValueError, 'linkage|linked'):
                    analyze(request(sequence=residue, sites=f'1,{site}', glycans_O='', observations=f'{mz},1', **spec))

    def test_both_site_types_keep_identical_mass_references_separate(self):
        # Two real records have exactly the same composition mass but different
        # annotated linkage. A shared mass bucket must not mix their IDs.
        n, o = 'G00017IP', 'G04143WH'
        mass = record(n)['free_mass']
        self.assertEqual(mass, record(o)['free_mass'])
        for convention in ('free', 'attached'):
            for polarity, sign in [('positive',1), ('negative',-1)]:
                observations = '\n'.join(f'{MONO[aa]+mass+sign*PROTON},1' for aa in 'NS')
                result = analyze(request(sequence='NS', sites='1,N\n2,O', glycans_O='',
                                         glycosmos_selection={'ids':[n,o]}, glycan_mass_type=convention,
                                         ion_mode=polarity, observations=observations))
                self.assertEqual(result['candidate_count'], 2)
                for measurement, kind, identifier in zip(result['results'], 'NO', (n,o)):
                    self.assertEqual(measurement['match_count'], 1)
                    self.assertEqual(measurement['rows'][0]['type'], kind)
                    for row in measurement['rows']:
                        expected = n if row['type']=='N' else o
                        self.assertEqual(row['glycan_info']['ids'], [expected])
                        self.assertEqual(row['glycan_info']['library_count'], 1)
                        self.assertAlmostEqual(row['glycan'], mass if convention=='free' else mass-WATER)

    def test_mixed_sites_skip_unavailable_type_without_reassigning_glycan(self):
        result = analyze(request(sequence='NS', sites='1,N\n2,O', glycans_O='',
                                 glycosmos_selection={'ids':['G00026MO']}))
        self.assertEqual(result['candidate_count'], 1)
        self.assertEqual(result['skipped_sites'], [{'site':2,'type':'O'}])
        self.assertTrue(all(row['type']=='N' for row in result['results'][0]['rows']))
        self.assertEqual(result['glycosmos_coverage']['matched'], 1)

    def test_site_selector_only_restricts_and_cannot_override_linkage(self):
        for kind, identifier in [('N','G00039MO'), ('O','G00026MO')]:
            with self.subTest(kind=kind), self.assertRaisesRegex(ValueError, 'No linkage-compatible'):
                analyze(request(sequence='NS', sites='1,N\n2,O', glycans_O='',
                                glycosmos_selection={'ids':[identifier]}, glycosmos_apply_to=kind))
        for kind in 'NO':
            result = analyze(request(sequence='NS', sites='1,N\n2,O', glycans_O='',
                                     glycosmos_selection={'ids':['G00026MO','G00039MO']}, glycosmos_apply_to=kind))
            self.assertEqual(result['glycosmos_coverage']['incompatible'], 1)
            self.assertTrue(all(row['type']==kind for row in result['results'][0]['rows']))

    def test_unclassified_mass_is_available_only_without_a_protein_assignment(self):
        identifier = 'G00001UD'
        self.assertEqual(record(identifier)['linkage_types'], [])
        data = request(glycans_O='', glycosmos_selection={'ids':[identifier]},
                       observations=f"{record(identifier)['free_mass']+PROTON},1")
        with self.assertRaisesRegex(ValueError, 'No linkage-compatible'):
            analyze(data)
        standalone = analyze_glycan_adducts(data)
        self.assertEqual(standalone['results'][0]['best_assignment']['glycan_info']['ids'], [identifier])
        self.assertTrue(standalone['results'][0]['best_assignment']['within'])
        self.assertEqual(standalone['glycosmos_coverage']['matched'], 1)
        data['glycosmos_selection']['ids'].append('G00039MO')
        peptide = analyze(data)
        self.assertEqual(peptide['glycosmos_coverage']['unclassified'], 1)
        self.assertEqual(peptide['glycosmos_coverage']['matched'], 1)
        self.assertEqual(peptide['results'][0]['rows'][0]['glycan_info']['ids'], ['G00039MO'])

    def test_full_library_groups_and_results_have_only_compatible_references(self):
        groups, coverage = selection_groups({'all':True}, 'both', ('N','O'), 'mono_free')
        for kind, masses in groups.items():
            for group in masses.values():
                self.assertTrue(set(group['ids']) <= linkage_index()[kind])
        self.assertEqual(coverage['supported'], coverage['matched']+coverage['linkage_skipped'])
        self.assertEqual(coverage['matched'], 28942)
        self.assertEqual(coverage['unclassified'], 87147)
        result = analyze(request(sequence='NS', sites='1,N\n2,O', glycans_O='',
                                 glycosmos_selection={'all':True}, observations='500,1\n1500,1\n3000,1'))
        for measurement in result['results']:
            for row in measurement['rows']:
                self.assertTrue(set(row['glycan_info']['ids']) <= linkage_index()[row['type']])

    def test_custom_mass_does_not_attach_wrong_library_reference_to_same_mass(self):
        mass = record('G00026MO')['free_mass']
        result = analyze(request(glycans_O=str(mass), observations=f'{MONO["S"]+mass+PROTON},1',
                                 glycosmos_selection={'ids':['G00026MO']}))
        row = result['results'][0]['rows'][0]
        self.assertEqual(row['type'], 'O')
        self.assertTrue(row['glycan_info']['custom'])
        self.assertEqual(row['glycan_info']['ids'], [])
        self.assertEqual(result['glycosmos_coverage']['incompatible'], 1)

    def test_api_enforces_bulk_and_legacy_selection_and_exports_coverage(self):
        client = test_vercel.DeploymentTests()
        for spec in ({'glycosmos_selection':{'ids':['G00026MO']}},
                     {'glycosmos_glycans':[{'id':'G00026MO','type':'O'}]}):
            code, _ = client.request('/api/analyze', request(glycans_O='', **spec))
            self.assertEqual(code, 400)
        code, data = client.request('/api/analyze', request(glycans_O='', glycosmos_selection={'ids':['G00026MO','G00039MO']}))
        self.assertEqual(code, 200)
        self.assertEqual(data['glycosmos_coverage']['matched'], 1)
        self.assertEqual(data['glycosmos_coverage']['incompatible'], 1)
        self.assertEqual(data['results'][0]['rows'][0]['glycan_info']['ids'], ['G00039MO'])
