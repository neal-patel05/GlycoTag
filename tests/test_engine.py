import unittest
from engine import analyze, sequence, PROTON, WATER, MONO, FREE_AA, EXAMPLE

class MatchingTests(unittest.TestCase):
    def data(self, **kwargs):
        d=dict(sequence='ASNTA',sites='2,O\n3,N',observations='500,2',glycans_O='674.2382',glycans_N='2222.7830',min_length=1,max_length=5,mode='user',ion_mode='negative',tolerance=20,unit='ppm')
        d.update(kwargs)
        return d

    def test_user_formula_both_sides_and_charge(self):
        mass=674.2382+sum(FREE_AA[x]-WATER for x in 'ASN')
        result=analyze(self.data(observations=f'{mass/3-PROTON},3'))['results'][0]
        best=result['rows'][0]
        self.assertEqual((best['sequence'],best['start'],best['end'],best['site']),('ASN',1,3,2))
        self.assertAlmostEqual(best['error_da'],0,places=9)
        self.assertTrue(best['within'])

    def test_mono_water_conventions(self):
        mass=674.2382+MONO['S']
        for mode,water in [('mono_free',0),('mono_attached',WATER)]:
            r=analyze(self.data(sequence='S',sites='1,O',mode=mode,observations=f'{mass+water-PROTON},1'))['results'][0]
            self.assertAlmostEqual(r['rows'][0]['error_da'],0,places=8)
            self.assertEqual(r['match_count'],1)

    def test_matches_bruteforce(self):
        data=self.data(tolerance=10000,unit='Da')
        result=analyze(data)
        expected=[]
        for p,k,g in [(2,'O',674.2382),(3,'N',2222.7830)]:
            for a in range(5):
                for b in range(a+1,6):
                    if a<p<=b:
                        expected.append(g+sum(FREE_AA[x]-WATER for x in 'ASNTA'[a:b]))
        self.assertEqual(result['candidate_count'],len(expected))
        self.assertEqual(result['results'][0]['match_count'],len(expected))
        target=2*500
        errors=sorted(abs(x-2*PROTON-target) for x in expected)
        for row,error in zip(result['results'][0]['rows'],errors):
            self.assertAlmostEqual(abs(row['error_da']),error)

    def test_boundaries_batch_and_no_match(self):
        r=analyze(self.data(sequence='SNS',sites='1,O\n3,O',max_length=1,observations='1.1,1\n1.2,2'))
        self.assertEqual(r['candidate_count'],2)
        self.assertEqual(len(r['results']),2)
        self.assertTrue(all(x['match_count']==0 for x in r['results']))

    def test_fasta(self):
        self.assertEqual(sequence('>protein\nas nt\nA'),'ASNTA')
        for s in ('>one\nAA\n>two\nSS','ASXB',''):
            with self.assertRaises(ValueError):sequence(s)

    def test_negative_ion_target_neutral_mass(self):
        raw_mass = 674.2382 + FREE_AA['S'] - WATER
        for charge in (1,2,5):
            mz = raw_mass/charge - PROTON + 0.01
            result = analyze(self.data(sequence='S',sites='1,O',
                observations=f'{mz},{charge}',unit='Da',tolerance=0.02))['results'][0]
            row = result['rows'][0]
            self.assertEqual(result['mz'],mz)
            self.assertEqual(result['target_mass'],mz*charge+charge*PROTON)
            self.assertAlmostEqual(row['mass'],raw_mass)
            self.assertAlmostEqual(row['predicted_mz'],raw_mass/charge-PROTON)
            self.assertAlmostEqual(row['error_da'],-0.01*charge)
            self.assertAlmostEqual(row['error_ppm'],-0.01/mz*1e6)
            if charge == 1:
                self.assertEqual(result['match_count'],1)
            elif charge == 5:
                self.assertEqual(result['match_count'],0)

    def test_validation(self):
        for change in [dict(sites='1,N'),dict(sites='6,O'),dict(observations='NaN,2'),dict(observations='500,2.5'),dict(tolerance=-1),dict(glycans_O='inf'),dict(min_length=4,max_length=2),dict(unit='bad'),dict(sites=''),dict(mode='bad')]:
            with self.subTest(change=change),self.assertRaises(ValueError):analyze(self.data(**change))

    def test_serotransferrin(self):
        self.assertEqual([EXAMPLE[p-1] for p in [51,432,491,630]],list('SNNN'))
        r=analyze(self.data(sequence=EXAMPLE,sites='51,O\n432,N\n491,N\n630,N',max_length=60))
        self.assertGreater(r['candidate_count'],1000)

if __name__=='__main__':unittest.main()
