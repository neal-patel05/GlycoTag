import unittest

from adducts import CATALOG
from glycan_adducts import analyze_glycan_adducts
from masses import WATER
import test_vercel


def request(**changes):
    data = dict(glycans_O='1000', observations='1022.989218,1',
                selected_adducts=[a['id'] for a in CATALOG], tolerance=1, unit='ppm')
    data.update(changes)
    return data


class GlycanAdductTests(unittest.TestCase):
    def test_requested_formula_and_signed_sodium_difference(self):
        r = analyze_glycan_adducts(request())['results'][0]
        self.assertAlmostEqual(r['calculated_mw'], 1021.989218)
        c = r['best_assignment']
        self.assertEqual(c['glycan_mw'], 1000)
        self.assertAlmostEqual(c['mass_difference'], -21.989218)
        a = c
        self.assertEqual(a['adduct'], '[M+Na]+')
        self.assertAlmostEqual(a['expected_difference'], -21.989218)
        self.assertAlmostEqual(a['difference_error'], 0)
        self.assertTrue(a['within'])

    def test_multiple_charges_negative_mode_and_multimer(self):
        # Independent examples from the source table, including the sign of
        # the requested subtraction in negative mode and non-unit charges.
        for name, mz, z, mode in [('[M+2Na]2+',522.989218,2,'positive'),
                                 ('[M-H]-',998.992723533379,1,'negative'),
                                 ('[3M-H]-',2998.992723533379,1,'negative')]:
            with self.subTest(adduct=name):
                r = analyze_glycan_adducts(request(observations=f'{mz},{z}', ion_mode=mode))['results'][0]
                self.assertAlmostEqual(r['calculated_mw'], mz*z-z)
                c = r['best_assignment']
                self.assertAlmostEqual(c['mass_difference'],1000-(mz*z-z))
                self.assertEqual(c['assignment_type'],'glycan_only' if name=='[M-H]-' else 'adduct_associated')
                self.assertEqual(c['adduct'],None if name=='[M-H]-' else name)
                self.assertTrue(c['within'])
                self.assertAlmostEqual(c['difference_error'],0,places=9)

    def test_best_assignment_competes_across_all_selected_glycans(self):
        r = analyze_glycan_adducts(request(glycans_O='1000,1022'))['results'][0]
        winner = r['best_assignment']
        self.assertEqual(winner['glycan_mw'],1000)
        self.assertEqual(winner['assignment_type'],'adduct_associated')
        self.assertEqual(winner['adduct'],'[M+Na]+')
        self.assertTrue(winner['within'])
        self.assertNotIn('candidates',r)

    def test_glycan_only_competes_even_when_only_sodium_selected(self):
        sodium = next(a['id'] for a in CATALOG if a['name']=='[M+Na]+')
        r = analyze_glycan_adducts(request(observations='1001.007276466621,1',selected_adducts=[sodium]))['results'][0]
        winner = r['best_assignment']
        self.assertEqual(winner['assignment_type'],'glycan_only')
        self.assertIsNone(winner['adduct'])
        self.assertTrue(winner['within'])
        self.assertEqual(r['tied_assignment_count'],1)

    def test_tie_is_flagged_with_one_representative_and_no_duplicate_proton_ions(self):
        # Sodium on one glycan and protonation of another give identical m/z.
        other = 1022.989218 - 1.007276466621
        r = analyze_glycan_adducts(request(glycans_O=f'1000,{other}'))['results'][0]
        self.assertEqual(r['tied_assignment_count'],2)
        self.assertEqual(r['best_assignment']['assignment_type'],'glycan_only')
        self.assertAlmostEqual(r['best_assignment']['glycan_mw'],other)
        self.assertIsNone(r['best_assignment']['adduct'])
        r = analyze_glycan_adducts(request(observations='1001.007276466621,1'))['results'][0]
        self.assertEqual(r['tied_assignment_count'],1)

    def test_repeated_observations_get_one_assignment_per_mz_and_charge(self):
        result = analyze_glycan_adducts(request(observations='1022.989218,1\n1022.9892180,1\n1022.989218,2'))
        self.assertEqual(len(result['results']),2)
        self.assertEqual([r['charge'] for r in result['results']],[1,2])

    def test_attached_masses_convert_to_free(self):
        r = analyze_glycan_adducts(request(glycans_O=str(1000-WATER),glycan_mass_type='attached'))['results'][0]
        self.assertAlmostEqual(r['best_assignment']['glycan_mw'],1000)
        self.assertTrue(r['best_assignment']['within'])

    def test_in_tolerance_assignment_preferred_for_multimer_da_tolerance(self):
        # A monomer has smaller m/z error but fails the neutral-monomer Da
        # tolerance. The dimer fits that tolerance and must win.
        proton = 1.007276466621
        monomer = 999.98 - proton
        dimer = (1000.03 - proton)/2
        ion = next(a['id'] for a in CATALOG if a['name']=='[2M+H]+')
        r = analyze_glycan_adducts(request(glycans_O=f'{monomer},{dimer}',
            observations='1000,1',selected_adducts=[ion],unit='Da',tolerance=.016))['results'][0]
        self.assertEqual(r['best_assignment']['adduct'],'[2M+H]+')
        self.assertTrue(r['best_assignment']['within'])
        self.assertAlmostEqual(r['best_assignment']['error_da'],.015)

    def test_glycosmos_filter_references_and_images_without_protein(self):
        result = analyze_glycan_adducts(request(glycans_O='',glycosmos_selection={'ids':['G00039MO'],'linkage':'O'},
                                               observations='447.1585726,1'))
        c = result['results'][0]['best_assignment']
        self.assertEqual(c['glycan_info']['ids'],['G00039MO'])
        self.assertEqual(result['glycosmos_image_ids'],['G00039MO'])
        self.assertEqual(c['adduct'],'[M+Na]+')
        self.assertTrue(c['within'])
        self.assertEqual(result['glycosmos_coverage']['selected'],1)
        result = analyze_glycan_adducts(request(glycans_O='',glycosmos_selection={'all':True,'linkage':'N'}))
        self.assertEqual(result['glycosmos_coverage']['selected'],22667)
        self.assertGreater(result['unique_masses'],0)

    def test_tolerance_and_incompatible_charge(self):
        for unit, tolerance, within in [('mz',.02,True),('mz',.001,False),
                                        ('Da',.02,True),('Da',.001,False),
                                        ('ppm',20,True),('ppm',1,False)]:
            r = analyze_glycan_adducts(request(observations='1022.999218,1',unit=unit,tolerance=tolerance))['results'][0]
            self.assertEqual(r['best_assignment']['within'],within)
        r = analyze_glycan_adducts(request(observations='251,4'))['results'][0]
        self.assertEqual(r['calculated_mw'],1000)
        self.assertEqual(r['best_assignment']['assignment_type'],'glycan_only')
        self.assertIsNone(r['best_assignment']['adduct'])
        self.assertFalse(r['best_assignment']['within'])

    def test_validation_and_hosted_endpoint(self):
        for changes in [dict(observations='1,1'),dict(observations='2,0'),dict(observations='nan,1'),
                        dict(observations='1e308,100'),dict(glycans_O=''),dict(ion_mode='invalid'),
                        dict(selected_adducts=[]),dict(glycan_mass_type='invalid')]:
            with self.subTest(changes=changes),self.assertRaises(ValueError):
                analyze_glycan_adducts(request(**changes))
        client = test_vercel.DeploymentTests()
        code, data = client.request('/api?route=glycan-adducts', request())
        self.assertEqual(code,200)
        self.assertAlmostEqual(data['results'][0]['calculated_mw'],1021.989218)


if __name__ == '__main__':
    unittest.main()
