import unittest
from engine import analyze, PROTON, WATER

class IonModeTests(unittest.TestCase):
    def data(self, **changes):
        data = dict(sequence='AS', sites='2,O', min_length=2, max_length=2,
                    glycan_mass_type='free', glycans_O='221.0899546',
                    observations='200,4', ion_mode='positive', unit='ppm', tolerance=20)
        data.update(changes)
        return data

    def test_constants(self):
        self.assertEqual(PROTON, 1.007276466621)
        self.assertEqual(WATER, 18.010564684)

    def test_both_polarities_mass_types_charges_and_errors(self):
        residues = 71.03711378515 + 87.03202840472
        peptide = residues + 18.010564684
        free = 221.0899546
        attached = free - 18.010564684
        neutral = peptide + attached
        for mode, sign in [('positive',1), ('negative',-1)]:
            for kind, glycan in [('free',free), ('attached',attached)]:
                for z in (1,2,4,7):
                    for offset in (0,0.01,-0.01):
                        with self.subTest(mode=mode,kind=kind,z=z,offset=offset):
                            predicted = neutral/z+sign*1.007276466621
                            observed = predicted+offset
                            result = analyze(self.data(ion_mode=mode,glycan_mass_type=kind,
                                glycans_O=str(glycan), observations=f'{observed},{z}'))['results'][0]
                            row = result['rows'][0]
                            self.assertEqual(result['charge'],z)
                            self.assertAlmostEqual(row['mass'],neutral,places=10)
                            self.assertAlmostEqual(result['target_neutral_mass'],observed*z-sign*z*PROTON,places=10)
                            self.assertAlmostEqual(row['predicted_mz'],predicted,places=10)
                            self.assertAlmostEqual(row['signed_delta_mz'],-offset,places=10)
                            self.assertAlmostEqual(row['absolute_delta_mz'],abs(offset),places=10)
                            self.assertAlmostEqual(row['error_ppm'],-offset/observed*1e6,places=7)
                            self.assertAlmostEqual(row['error_da'],-offset*z,places=10)
                            b = row['breakdown']
                            self.assertAlmostEqual(b['peptide_mass'],peptide,places=10)
                            self.assertEqual(b['aa_water_loss'],0)
                            self.assertEqual(b['attachment_water'],WATER if kind=='free' else 0)
                            self.assertAlmostEqual(b['attached_glycan_mass'],attached,places=10)

    def test_attached_composition_is_not_dehydrated_again(self):
        row = analyze(self.data(glycan_mass_type='attached',glycans_O='203.0794'))['results'][0]['rows'][0]
        self.assertAlmostEqual(row['mass'],71.03711378515+87.03202840472+18.010564684+203.0794,places=10)
        self.assertEqual(row['breakdown']['attachment_water'],0)

    def test_tolerance_units_and_nearest_order(self):
        for mode,sign in [('positive',1),('negative',-1)]:
            mz = (71.03711378515+87.03202840472+221.0899546)/4+sign*PROTON+0.01
            for unit,tolerance,within in [('mz',0.011,True),('mz',0.009,False),('Da',0.039,False),('Da',0.041,True),('ppm',0.011/mz*1e6,True),('ppm',0.009/mz*1e6,False)]:
                result=analyze(self.data(ion_mode=mode,observations=f'{mz},4',unit=unit,tolerance=tolerance))['results'][0]
                self.assertEqual(result['match_count'],int(within))
                self.assertEqual(result['rows'][0]['within'],within)

    def test_validation(self):
        for changes in [dict(ion_mode='bad'),dict(glycan_mass_type='bad'),dict(observations='200,-4'),dict(observations='200,0'),dict(observations='200,1.5'),dict(observations='0.5,4'),dict(glycans_O='18')]:
            with self.subTest(changes=changes),self.assertRaises(ValueError): analyze(self.data(**changes))
