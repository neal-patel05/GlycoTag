import json
from pathlib import Path
import unittest

from adducts import CATALOG, hypotheses
from engine import analyze, MONO, PROTON, EXAMPLE
from glycosmos import (browse, composition_mass, resolve_selections, RESIDUES,
                      structure_image, record, linkage_types, linkage_index, selection_groups)
from glycans import EXPASY_REDUCING_END
from masses import WATER
import test_vercel


def request(**changes):
    data = dict(sequence='S', sites='1,O', glycans_O='221.0899546', observations='400,1',
                min_length=1, max_length=1, tolerance=1, unit='ppm')
    data.update(changes)
    return data


class AdductTests(unittest.TestCase):
    def test_full_source_catalogue(self):
        self.assertEqual(sum(a['mode']=='positive' for a in CATALOG), 31)
        self.assertEqual(sum(a['mode']=='negative' for a in CATALOG), 15)
        self.assertEqual(len({a['id'] for a in CATALOG}), 46)
        mixed = next(a for a in CATALOG if a['name']=='[M+H+2Na]3+')
        self.assertAlmostEqual(mixed['mz_offset'], 15.662453333333334)

    def test_all_adducts_reconcile_both_mass_models(self):
        neutral = MONO['S'] + 221.0899546
        for entry in CATALOG:
            ion = hypotheses([entry], entry['mode'], entry['charge'])[0]
            mz = (ion['multimer']*neutral+ion['mass_shift'])/ion['charge']
            for kind, mass in [('free',221.0899546), ('attached',221.0899546-WATER)]:
                with self.subTest(adduct=entry['name'],kind=kind):
                    result = analyze(request(ion_mode=entry['mode'],selected_adducts=[entry['id']],
                                             observations=f'{mz},{ion["charge"]}',glycan_mass_type=kind,glycans_O=str(mass)))['results'][0]
                    row = result['rows'][0]
                    self.assertEqual(result['match_count'],1)
                    self.assertAlmostEqual(row['target_neutral_mass'],neutral,places=10)
                    self.assertAlmostEqual(row['predicted_mz'],mz,places=10)
                    self.assertAlmostEqual(row['breakdown']['calculated_mass'],neutral,places=10)
                    self.assertEqual(row['adduct'],entry['name'])

    def test_independent_sodium_and_trimer_calculations(self):
        neutral = MONO['S'] + 221.0899546
        sodium = next(a for a in CATALOG if a['name']=='[M+Na]+')
        trimer = next(a for a in CATALOG if a['name']=='[3M-H]-')
        for a,mz in [(sodium,neutral+22.989218),(trimer,3*neutral-PROTON)]:
            r = analyze(request(ion_mode=a['mode'],selected_adducts=[a['id']],observations=f'{mz},1'))['results'][0]['rows'][0]
            self.assertAlmostEqual(r['error_da'],0,places=10)

    def test_multimer_tolerances_are_scaled(self):
        a = next(a for a in CATALOG if a['name']=='[3M-H]-')
        mz = 3*(MONO['S']+221.0899546)-PROTON + .03
        for unit,tolerance,within in [('Da',.011,True),('Da',.009,False),('mz',.031,True),('mz',.029,False),('ppm',.031/mz*1e6,True),('ppm',.029/mz*1e6,False)]:
            r = analyze(request(ion_mode='negative',selected_adducts=[a['id']],observations=f'{mz},1',unit=unit,tolerance=tolerance))['results'][0]
            self.assertEqual(r['match_count'],int(within))
            self.assertEqual(r['rows'][0]['within'],within)

    def test_select_all_ranks_by_mz_and_respects_polarity_and_charge(self):
        r = analyze(request(selected_adducts=[a['id'] for a in CATALOG]))['results'][0]
        errors = [row['absolute_delta_mz'] for row in r['rows']]
        self.assertEqual(errors,sorted(errors))
        self.assertIsNone(r['target_neutral_mass'])
        positive = {a['id'] for a in CATALOG if a['mode']=='positive' and a['charge']==1}
        self.assertTrue(all(row['adduct_id'] in positive for row in r['rows']))
        r = analyze(request(selected_adducts=[a['id'] for a in CATALOG],observations='400,4'))['results'][0]
        self.assertEqual(r['rows'],[])
        self.assertIn('No selected adduct',r['note'])

    def test_invalid_selections(self):
        for ids in [[],['unknown'],['pos_1',1],'pos_1',['neg_1']]:
            with self.subTest(ids=ids),self.assertRaises(ValueError):
                analyze(request(selected_adducts=ids))


class GlyCosmosTests(unittest.TestCase):
    def test_source_linkage_filters(self):
        # Known source annotations, including uncommon O-linked motifs and
        # non-linkage "N-" names that must never be mistaken for N-glycans.
        self.assertEqual(linkage_types('N-acetyllactosamine:G00055MO,N-Glycolyl-GM2:G56478LD'), [])
        self.assertEqual(linkage_types('N-Glycan core basic:G00026MO'), ['N'])
        for identifier in ['G00039MO', 'G36123IU', 'G06898ZN', 'G06953DR', 'G82576YO']:
            self.assertEqual(record(identifier)['linkage_types'], ['O'])
        self.assertEqual(record('G98562EA')['linkage_types'], [])
        self.assertEqual(linkage_types('N-Glycan core basic:G00026MO,O-Glycan core 5:G00039MO'), ['N', 'O'])
        expected = {'all':257024, 'N':22667, 'O':8850, 'both':31517}
        for kind, count in expected.items():
            result = browse({'linkage':[kind]})
            self.assertEqual(result['total'],count)
            self.assertEqual(result['linkage_counts'],expected)
            if kind != 'all':
                self.assertTrue(all(set(g['linkage_types']) & (set('NO') if kind == 'both' else {kind})
                                    for g in result['glycans']))
        self.assertEqual(browse({'linkage':['N'],'q':['G00039MO']})['total'],0)
        filtered = browse({'linkage':['O'],'q':['G00039MO'],'supported':['1']})
        self.assertIn('G00039MO', [g['id'] for g in filtered['glycans']])
        self.assertTrue(all('O' in g['linkage_types'] and g['free_mass'] is not None
                            for g in filtered['glycans']))
        self.assertEqual(linkage_index()['both'],linkage_index()['N'] | linkage_index()['O'])
        with self.assertRaises(ValueError): browse({'linkage':['invalid']})

    def test_linkage_scoped_selections_and_exclusions(self):
        for kind, count in [('N',22667), ('O',8850), ('both',31517)]:
            groups, coverage = selection_groups({'all':True,'linkage':kind},'both',{'N','O'},'mono_free')
            self.assertEqual(coverage['selected'],count)
            self.assertEqual(coverage['linkage'],kind)
            self.assertTrue(all(set(group['ids']) <= linkage_index()[kind]
                                for group in groups['O'].values()))
        _, coverage = selection_groups({'all':True,'linkage':'O','excluded':['G00039MO','G00026MO']},'both',{'O'},'mono_free')
        self.assertEqual(coverage['selected'],8849)
        groups, coverage = selection_groups({'ids':['G00039MO','G00026MO','G98562EA'],'linkage':'O'},'O',{'O'},'mono_free')
        self.assertEqual(coverage['selected'],1)
        self.assertEqual(next(iter(groups['O'].values()))['ids'],['G00039MO'])
        with self.assertRaises(ValueError): selection_groups({'all':True,'linkage':[]},'both',{'O'},'mono_free')

    def test_result_image_references_and_record_endpoint(self):
        result = analyze(request(glycans_O='',glycosmos_selection={'ids':['G00039MO']}))
        self.assertEqual(result['glycosmos_image_ids'],['G00039MO'])
        client=test_vercel.DeploymentTests()
        code, data=client.request('/api?route=glycosmos/record/G00039MO')
        self.assertEqual(code,200)
        self.assertTrue(data['has_image'])
        self.assertEqual(data['image_url'],'/api/glycosmos/image/G00039MO')
        self.assertEqual(data['linkage_types'],['O'])
        code, data=client.request('/api?route=glycosmos/record/G00000CV')
        self.assertEqual(code,200)
        self.assertFalse(data['has_image'])
        code, _=client.request('/api?route=glycosmos/record/invalid')
        self.assertEqual(code,400)
        code, data=client.request('/api?route=glycosmos&linkage=N&q=G00039MO')
        self.assertEqual(code,200)
        self.assertEqual(data['total'],0)

    def test_names_subunits_and_local_image(self):
        g=next(g for g in browse({'q':['G00039MO']})['glycans'] if g['id']=='G00039MO')
        self.assertEqual(g['name'],'GalNAc(a1-3)GalNAc(a1-')
        self.assertEqual(g['subunits'][0]['symbol'],'GalNAc')
        self.assertEqual(g['subunits'][0]['count'],2)
        self.assertIn('galactosamine',g['subunits'][0]['name'])
        self.assertTrue(g['has_image'])
        self.assertTrue(structure_image(g['id']).startswith(b'\x89PNG'))
        self.assertIsNone(structure_image('G00000CV'))
        with self.assertRaises(ValueError): structure_image('../data')

    def test_select_all_search_and_individual_exclusions(self):
        d=request(sequence=EXAMPLE,sites='51,O\n432,N\n491,N\n630,N',glycans_O='',max_length=60,
                  glycosmos_selection={'all':True},glycosmos_apply_to='both')
        r=analyze(d)
        self.assertEqual(r['glycosmos_coverage']['selected'],257024)
        self.assertEqual(r['glycosmos_coverage']['supported'],116089)
        self.assertEqual(r['glycosmos_coverage']['skipped'],140935)
        self.assertGreater(r['candidate_count'],1500000)
        self.assertEqual(len(r['results'][0]['rows']),25)
        self.assertLess(len(json.dumps(r)),100000)
        d['glycosmos_selection']={'all':True,'excluded':['G00039MO']}
        r=analyze(d)
        self.assertEqual(r['glycosmos_coverage']['selected'],257023)
        self.assertEqual(r['glycosmos_coverage']['supported'],116088)

    def test_selection_includes_records_without_n_o_classification_or_mass(self):
        r=analyze(request(glycans_O='',glycosmos_selection={'ids':['G00039MO','G98562EA']}))
        self.assertEqual(r['glycosmos_coverage']['selected'],2)
        self.assertEqual(r['glycosmos_coverage']['skipped'],1)
        self.assertEqual(r['results'][0]['rows'][0]['glycan_info']['ids'],['G00039MO'])
        for spec in [{'all':'yes'},{'ids':['bad']},{'excluded':[1]}]:
            with self.assertRaises(ValueError): analyze(request(glycosmos_selection=spec))

    def test_indexed_search_matches_cartesian_bruteforce(self):
        seq='ASNTS'
        for polarity in ['positive','negative']:
            entries=[a for a in CATALOG if a['mode']==polarity]
            observed=500
            d=request(sequence=seq,sites='2,O\n3,N\n5,O',glycans_O='200,400,600',glycans_N='300,500',
                      max_length=5,observations='500,1',ion_mode=polarity,
                      selected_adducts=[a['id'] for a in entries],unit='mz',tolerance=25)
            r=analyze(d)['results'][0]
            expected=[]
            for ion in hypotheses(entries,polarity,1):
                for site,kind,values in [(2,'O',[200,400,600]),(3,'N',[300,500]),(5,'O',[200,400,600])]:
                    for start in range(site):
                        for end in range(site,6):
                            for glycan in values:
                                mass=sum(MONO[a] for a in seq[start:end])+glycan
                                mz=ion['multimer']*mass+ion['mass_shift']
                                expected.append((abs(mz-observed),mz))
            expected.sort()
            self.assertEqual(r['match_count'],sum(e<=25 for e,_ in expected))
            for row,(error,mz) in zip(r['rows'],expected[:25]):
                self.assertAlmostEqual(row['absolute_delta_mz'],error,places=9)
                self.assertAlmostEqual(row['predicted_mz'],mz,places=9)

    def test_expasy_all_ten_inputs(self):
        fixture=json.loads((Path(__file__).parent/'fixtures/expasy-expanded.json').read_text())
        actual=sum(RESIDUES[s]*n for s,n in fixture['composition'].items())+EXPASY_REDUCING_END
        self.assertEqual(float(actual),fixture['free_mass'])

    def test_wurcs_reference_composition_and_unknown_stereochemistry(self):
        wurcs='WURCS=2.0/3,5,4/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]/1-1-2-3-3/a4-b1_b4-c1_c3-d1_c6-e1'
        result=composition_mass(wurcs,910.33)
        self.assertEqual(result['composition'],{'HexNAc':2,'Hex':3})
        self.assertAlmostEqual(result['free_mass'],910.3277546)
        r=composition_mass('WURCS=2.0/2,3,0+/[uxxxxh][uxxxxh_2*NCC/3=O]/1-1-2/')
        self.assertEqual(r['composition'],{'Hex':2,'HexNAc':1})

    def test_rejects_unknown_and_nonreducing_chemistry(self):
        for w in ['', 'WURCS=2.0/1,1,0/[h2122h]/1/', 'WURCS=2.0/1,1,0/[a2122h-1b_1-5_2*N]/1/',
                  'WURCS=2.0/1,2,2/[a2122h-1b_1-5]/1-1/a1-b4_a4-b1',
                  'WURCS=2.0/1,2,2/[a2122h-1b_1-5]/1-1/a1-b4~n']:
            with self.subTest(w=w),self.assertRaises(ValueError): composition_mass(w)
        with self.assertRaises(ValueError): composition_mass('WURCS=2.0/1,1,0/[uxxxxh]/1/',999)

    def test_full_snapshot_pagination_filter_and_literal_search(self):
        r=browse({})
        self.assertEqual(r['total'],257024)
        self.assertEqual(len(r['glycans']),30)
        self.assertEqual(r['metadata']['supported'],116089)
        next_page=browse({'page':['2']})
        self.assertTrue({g['id'] for g in r['glycans']}.isdisjoint(g['id'] for g in next_page['glycans']))
        supported=browse({'supported':['1']})
        self.assertEqual(supported['total'],116089)
        self.assertTrue(all(g['free_mass'] for g in supported['glycans']))
        self.assertEqual(browse({'q':["%' OR 1=1 --"]})['total'],0)
        with self.assertRaises(ValueError): browse({'page':['0']})

    def test_real_glycan_matches_and_mass_conversion(self):
        selection=[dict(id='G00039MO',type='O')]
        g=resolve_selections(selection)[0]
        self.assertEqual(g['free_mass'],424.1693546)
        mz=MONO['S']+424.1693546+PROTON
        for kind in ['free','attached']:
            r=analyze(request(glycans_O='',glycosmos_glycans=selection,observations=f'{mz},1',glycan_mass_type=kind))['results'][0]
            self.assertEqual(r['match_count'],1)
            self.assertAlmostEqual(r['rows'][0]['error_da'],0)
            self.assertEqual(r['rows'][0]['glycan_info']['ids'],['G00039MO'])

    def test_rejects_missing_mass_and_bad_linkage(self):
        for entries in [[dict(id='G00039MO',type='X')],[dict(id='G98562EA',type='O')],[dict(id='bad',type='N')],['G00039MO']]:
            with self.subTest(entries=entries),self.assertRaises(ValueError): resolve_selections(entries)

    def test_hosted_query_parameters_preserved(self):
        client=test_vercel.DeploymentTests()
        code,data=client.request('/api?route=glycosmos&q=G00039MO&supported=1')
        self.assertEqual(code,200)
        self.assertTrue(any(g['id']=='G00039MO' for g in data['glycans']))
        self.assertTrue(all(g['free_mass'] for g in data['glycans']))
        code,data=client.request('/api?route=adducts')
        self.assertEqual(code,200)
        self.assertEqual(len(data['adducts']),46)
