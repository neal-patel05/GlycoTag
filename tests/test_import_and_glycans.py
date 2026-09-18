import copy
import json
from pathlib import Path
import unittest
from unittest.mock import patch
from engine import analyze, mass_breakdown, MONO, FREE_AA, PROTON, WATER
from glycans import CATALOG, BY_ID, SUGARS, model_mass
from uniprot import parse_record, import_entry

FIXTURE = json.loads((Path(__file__).parent/'fixtures'/'P02787-glycosylation.json').read_text())

class UniProtTests(unittest.TestCase):
    def test_real_annotations_include_atypical_and_evidence(self):
        record = parse_record(FIXTURE)
        self.assertEqual(record['sites_text'],'51,O\n432,N\n491,N\n630,N')
        self.assertIn('atypical',record['sites'][2]['descriptions'][0])
        self.assertEqual(record['sites'][2]['evidences'][0]['id'],'15536627')
        self.assertEqual(record['skipped'],[])

    def test_missing_annotations_are_not_predicted(self):
        record = copy.deepcopy(FIXTURE)
        record['features'] = []
        record['sequence']['value'] = 'NNSTSSSSNTT'
        self.assertEqual(parse_record(record)['sites_text'],'')

    def test_uncertain_ranges_wrong_isoforms_and_unsupported_residues(self):
        record = copy.deepcopy(FIXTURE)
        features = []
        for change in ['range','uncertain','wrong_residue','isoform','other_link']:
            f = copy.deepcopy(record['features'][0])
            if change == 'range': f['location']['end']['value'] = 55
            if change == 'uncertain': f['location']['start']['modifier'] = 'OUTSIDE'
            if change == 'wrong_residue': f['location']['start']['value'] = f['location']['end']['value'] = 1
            if change == 'isoform': f['location']['sequence'] = 'P02787-2'
            if change == 'other_link': f['description'] = 'C-linked mannose'
            features.append(f)
        record['features'] = features
        parsed = parse_record(record)
        self.assertEqual(parsed['sites'],[])
        self.assertEqual(len(parsed['skipped']),5)

    def test_isoform_never_inherits_wrong_coordinates(self):
        with patch('uniprot.download',side_effect=[json.dumps(FIXTURE),'>isoform\nNSST']):
            result = import_entry('P02787-2')
        self.assertEqual(result['sequence'],'NSST')
        self.assertEqual(result['sites'],[])
        self.assertTrue(result['warnings'])

    def test_same_sequence_isoform_retains_annotations(self):
        with patch('uniprot.download',side_effect=[json.dumps(FIXTURE),'>isoform\n'+FIXTURE['sequence']['value']]):
            self.assertEqual(len(import_entry('P02787-1')['sites']),4)

    def test_entry_id_and_validation(self):
        with patch('uniprot.download',return_value=json.dumps(FIXTURE)) as download:
            self.assertEqual(import_entry('trfe_human')['accession'],'P02787')
            download.assert_called_once_with('TRFE_HUMAN.json')
        with self.assertRaises(ValueError): import_entry('../invalid')

class GlycanTests(unittest.TestCase):
    def data(self, **kwargs):
        d=dict(sequence='S',sites='1,O',observations='500,2',glycans_O='',glycans_N='',selected_glycans=['o_core1'],mode='mono_free',min_length=1,max_length=1)
        d.update(kwargs)
        return d

    def test_reference_masses(self):
        self.assertAlmostEqual(BY_ID['o_sialyl_core1']['free_mass'],674.2382,places=4)
        self.assertAlmostEqual(BY_ID['o_disialyl_core1']['free_mass'],965.3336,places=4)
        self.assertAlmostEqual(BY_ID['n_a2g2s2']['free_mass'],2222.7830,places=4)
        self.assertAlmostEqual(BY_ID['n_a3g3s3']['free_mass'],2879.0106,places=4)
        self.assertAlmostEqual(BY_ID['o_core1']['attached_mass'],365.132196,places=5)

    def test_structure_composition_and_no_double_occupancy(self):
        for g in CATALOG:
            with self.subTest(g=g['id']):
                self.assertEqual(sum(g['composition'].values()),len(g['nodes']))
                self.assertAlmostEqual(g['free_mass']-g['attached_mass'],WATER)
                acceptors = set()
                for i,node in enumerate(g['nodes']):
                    self.assertIn(node['sugar'],SUGARS)
                    if i == 0: self.assertIsNone(node['parent'])
                    else:
                        self.assertLess(node['parent'],i)
                        attachment = (node['parent'],node['linkage'].split('–')[-1])
                        self.assertNotIn(attachment,acceptors)
                        acceptors.add(attachment)

    def test_library_conversion_preserves_neutral_mass(self):
        rows=[]
        for mode in ['mono_free','mono_attached']:
            row=analyze(self.data(mode=mode))['results'][0]['rows'][0]
            self.assertEqual(row['glycan_info']['ids'],['o_core1'])
            self.assertAlmostEqual(row['glycan'],model_mass(BY_ID['o_core1'],mode))
            rows.append(row)
        self.assertAlmostEqual(rows[0]['mass'],rows[1]['mass'])
        self.assertAlmostEqual(rows[0]['glycan']-rows[1]['glycan'],WATER)

    def test_custom_mass_unchanged_and_not_assigned_a_structure(self):
        for mode in ['user','mono_free','mono_attached']:
            r=analyze(self.data(mode=mode,selected_glycans=[],glycans_O='674.2382'))['results'][0]['rows'][0]
            self.assertEqual(r['glycan'],674.2382)
            self.assertEqual(r['glycan_info']['ids'],[])
            self.assertTrue(r['glycan_info']['custom'])

    def test_isobaric_reference_selections_keep_both_labels(self):
        r=analyze(self.data(selected_glycans=['o_tn','o_glcnac']))
        self.assertEqual(r['candidate_count'],1)
        self.assertEqual(set(r['results'][0]['rows'][0]['glycan_info']['ids']),{'o_tn','o_glcnac'})

    def test_unknown_selection_and_wrong_type_rejected(self):
        for selected in [['unknown'],'o_core1',[{}],['n_man5']]:
            with self.subTest(selected=selected),self.assertRaises(ValueError): analyze(self.data(selected_glycans=selected))

    def test_breakdown_reconciles_for_all_modes_and_multiple_AAs(self):
        seq='SASSTN'
        for mode in ['user','mono_free','mono_attached']:
            result=analyze(self.data(sequence=seq,min_length=len(seq),max_length=len(seq),mode=mode,observations='500,2'))['results'][0]
            r=result['rows'][0]; b=r['breakdown']
            self.assertEqual(sum(x['count'] for x in b['amino_acids']),len(seq))
            self.assertAlmostEqual(b['calculated_mass'],r['mass'])
            self.assertAlmostEqual(b['peptide_mass']+b['glycan_mass']-b['attachment_water'],r['mass'])
            self.assertAlmostEqual(r['corrected_mass']-result['target_mass'],r['error_da'])
            self.assertAlmostEqual(r['mass']/2-PROTON,r['predicted_mz'])
            if mode=='user':
                self.assertAlmostEqual(b['aa_input_total'],sum(FREE_AA[x] for x in seq))
                self.assertAlmostEqual(b['aa_water_loss'],len(seq)*WATER)
            else: self.assertAlmostEqual(b['aa_total'],sum(MONO[x] for x in seq))

if __name__=='__main__': unittest.main()
