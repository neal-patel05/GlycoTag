from decimal import Decimal
import unittest
from engine import analyze, mass_breakdown, PROTON
from masses import RESIDUE_DECIMAL, FREE_DECIMAL, RESIDUE_UNITS, WATER_DECIMAL, SCALE

class PrecisionTests(unittest.TestCase):
    def test_known_monoisotopic_values_and_all_twenty_residues(self):
        self.assertEqual(set(RESIDUE_DECIMAL),set('ARNDCEQGHILKMFPSTWYV'))
        for aa,mass in {'A':'71.03711378515','G':'57.02146372069',
                        'N':'114.04292744138','D':'115.02694302429',
                        'C':'103.00918495955','M':'131.04048508847',
                        'S':'87.03202840472','W':'186.07931295073'}.items():
            self.assertEqual(RESIDUE_DECIMAL[aa],Decimal(mass))
        self.assertEqual(WATER_DECIMAL,Decimal('18.01056468403'))
        for aa,mass in RESIDUE_DECIMAL.items():
            self.assertEqual(Decimal(RESIDUE_UNITS[aa])/SCALE,mass)
            self.assertEqual(FREE_DECIMAL[aa]-WATER_DECIMAL,mass)

    def test_breakdown_retains_all_AA_digits_in_every_mode(self):
        for mode in ('user','mono_free','mono_attached'):
            b=mass_breakdown('ARNDCEQGHILKMFPSTWYV',mode,674.2382)
            for aa in b['amino_acids']:
                self.assertEqual(Decimal(aa['residue_mass_text']),RESIDUE_DECIMAL[aa['aa']])
                self.assertEqual(len(aa['residue_mass_text'].split('.')[1]),11)
                self.assertEqual(Decimal(aa['input_mass_text'])-Decimal(aa['water_per_aa_text']),RESIDUE_DECIMAL[aa['aa']])

    def test_long_prefix_does_not_change_short_peptide_mass(self):
        masses=[]
        for seq,pos in [('S',1),('W'*9999+'S',10000)]:
            data=dict(sequence=seq,sites=f'{pos},O',observations='400,2',
                      glycans_O='674.2382',glycans_N='',min_length=1,max_length=1)
            for mode in ('user','mono_free'):
                row=analyze(dict(data,mode=mode))['results'][0]['rows'][0]
                self.assertEqual(row['mass'],float(RESIDUE_DECIMAL['S'])+674.2382)
                self.assertEqual(row['predicted_mz'],row['mass']/2-PROTON)
                masses.append(row['mass'])
        self.assertEqual(len(set(masses)),1)

if __name__=='__main__':unittest.main()
