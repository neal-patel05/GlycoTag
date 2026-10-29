import csv
import io
import json
import sqlite3
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

import glycosmos as g
from engine import analyze
from glycan_adducts import analyze_glycan_adducts
from scripts.build_glycosmos import build
from test_expanded_library import request
import test_vercel


class SpeciesTests(unittest.TestCase):
    def test_aliases_names_ids_and_partial_names(self):
        for names, count, taxon in [(('human', 'humans', '  HOMO   sapiens ', '9606'), 8016, '9606'),
                                     (('bovine', 'cow', 'cattle', 'Bos taurus', '9913'), 1755, '9913')]:
            for name in names:
                with self.subTest(name=name):
                    self.assertEqual(len(g.species_ids(name)), count)
                    self.assertEqual([t['id'] for t in g.matching_species(name)], [taxon])
        matches = g.matching_species('Escherichia coli')
        self.assertTrue(matches)
        self.assertTrue(all('Escherichia coli' in t['name'] for t in matches))
        for term in ('960', 'unknown organism', '%', "%' OR 1=1 --"):
            self.assertFalse(g.species_ids(term))
        self.assertTrue(all('_' in t['name'] for t in g.matching_species('_')))
        for term in (None, [], {}, 'a'*201):
            with self.assertRaises(ValueError):
                g.normalize_species(term)

    def test_species_combines_with_all_linkages_and_text_filters(self):
        expected = {'all':8016, 'N':2425, 'O':892, 'both':3317}
        for kind, count in expected.items():
            result = g.browse({'species':['human'], 'linkage':[kind]})
            self.assertEqual(result['total'], count)
            self.assertEqual(result['linkage_counts'], expected)
            self.assertTrue(all('9606' in {t['id'] for t in row['organisms']} for row in result['glycans']))
        # Non-N/O records are included without pretending they have a linkage.
        other = sorted(g.species_ids('human') - g.linkage_ids('both'))[0]
        self.assertEqual(g.browse({'species':['human'], 'q':[other]})['total'], 1)
        self.assertEqual(g.browse({'species':['human'], 'linkage':['both'], 'q':[other]})['total'], 0)
        filtered = g.browse({'species':['bovine'], 'q':['G00039MO'], 'supported':['1']})
        self.assertIn('G00039MO', [row['id'] for row in filtered['glycans']])
        self.assertTrue(all(row['free_mass'] is not None for row in filtered['glycans']))
        self.assertEqual(filtered['linkage_counts']['all'], 1755)
        self.assertEqual(g.browse({'species':['']})['total'], 257024)
        self.assertEqual(g.browse({'species':['not-an-organism']})['total'], 0)
        page1 = g.browse({'species':['human']})
        page2 = g.browse({'species':['human'], 'page':['2']})
        self.assertTrue({x['id'] for x in page1['glycans']}.isdisjoint(x['id'] for x in page2['glycans']))

    def test_source_organisms_and_metadata(self):
        organisms = g.record('G00039MO')['organisms']
        self.assertIn({'id':'9606', 'name':'Homo sapiens'}, organisms)
        self.assertIn({'id':'9913', 'name':'Bos taurus'}, organisms)
        self.assertIn({'id':'1446740', 'name':'Escherichia coli O45:H2 str. 01-3147'}, organisms)
        d = g.browse({})
        self.assertEqual(d['metadata']['taxonomy_annotated'], 21689)
        self.assertEqual(d['metadata']['taxonomy_unannotated'], 235335)
        self.assertEqual(d['metadata']['taxonomy_count'], 2027)
        self.assertTrue(any(t['name']=='Homo sapiens' for t in d['species_options']))

    def test_scoped_select_all_exclusions_and_explicit_ids(self):
        for kind in ('all', 'N', 'O', 'both'):
            allowed = g.species_ids('human') & g.linkage_ids(kind)
            groups, coverage = g.selection_groups({'all':True, 'species':'human', 'linkage':kind,
                                                   'excluded':['G00039MO','G98562EA']}, 'both', {'N','O'}, 'mono_free')
            self.assertEqual(coverage['selected'], len(allowed - {'G00039MO','G98562EA'}))
            self.assertEqual(coverage['selected'], coverage['supported']+coverage['skipped'])
            self.assertEqual(coverage['species_matches'][0]['id'], '9606')
            self.assertTrue(all(set(row['ids']) <= allowed - {'G00039MO'} for row in groups['O'].values()))
        _, coverage = g.selection_groups({'ids':['G00039MO','G98562EA'], 'species':'human'}, 'O', {'O'}, 'mono_free')
        self.assertEqual(coverage['selected'], 1)
        _, coverage = g.selection_groups({'all':True, 'species':'not-an-organism'}, 'O', {'O'}, 'mono_free')
        self.assertEqual(coverage['selected'], 0)

    def test_both_analysis_paths_keep_species_scope_and_provenance(self):
        data = request(glycans_O='', glycosmos_selection={'all':True, 'species':'bovine', 'linkage':'O'})
        for result in (analyze(data), analyze_glycan_adducts(data)):
            self.assertEqual(result['glycosmos_coverage']['selected'], 176)
            self.assertEqual(result['glycosmos_selection']['species'], 'bovine')
            self.assertEqual(result['glycosmos_coverage']['species_matches'][0]['name'], 'Bos taurus')
        rows = analyze(data)['results'][0]['rows']
        allowed = g.species_ids('bovine') & g.linkage_ids('O')
        self.assertTrue(all(set(row['glycan_info']['ids']) <= allowed for row in rows))
        winner = analyze_glycan_adducts(data)['results'][0]['best_assignment']
        self.assertTrue(set(winner['glycan_info']['ids']) <= allowed)

    def test_http_query_and_invalid_species(self):
        client = test_vercel.DeploymentTests()
        code, result = client.request('/api?route=glycosmos&species=bovine&linkage=O')
        self.assertEqual(code, 200)
        self.assertEqual(result['total'], 176)
        code, result = client.request('/api?route=glycosmos&species=Homo%20sapiens&q=G00039MO')
        self.assertEqual(code, 200)
        self.assertIn('G00039MO', [row['id'] for row in result['glycans']])
        code, _ = client.request('/api?route=glycosmos&species='+'a'*201)
        self.assertEqual(code, 400)

    def test_import_preserves_taxonomy_and_deduplicates_relationships(self):
        fields = ['Accession Number','WURCS','IUPAC Condensed','Monoisotopic mass','Motif name(s)','Taxonomy']
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            source, database = path/'source.csv', path/'library.sqlite3'
            with source.open('w', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=fields)
                writer.writeheader()
                writer.writerow({'Accession Number':'G00001AA', 'Taxonomy':'Homo sapiens:9606,Homo sapiens:9606,Escherichia coli O45:H2:1446740'})
                writer.writerow({'Accession Number':'G00002AA', 'Taxonomy':''})
            with redirect_stdout(io.StringIO()):
                build(source, database, retrieved='2026-10-26')
            with sqlite3.connect(database) as db:
                self.assertEqual(db.execute('SELECT count(*) FROM glycan_taxa').fetchone()[0], 2)
                self.assertEqual(db.execute('SELECT name FROM taxa WHERE id=?', ('1446740',)).fetchone()[0], 'Escherichia coli O45:H2')
                metadata = json.loads(db.execute('SELECT value FROM metadata').fetchone()[0])
                self.assertEqual(metadata['taxonomy_annotated'], 1)
                self.assertEqual(metadata['taxonomy_unannotated'], 1)
                self.assertEqual(metadata['retrieved'], '2026-10-26')
