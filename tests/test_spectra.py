import io
import json
import unittest
from app import Handler
from spectra import import_peak_table, parse_peak_table


class PeakTableTests(unittest.TestCase):
    def test_qual_style_table_preserves_peak_charge_pairing_and_unknowns(self):
        data = parse_peak_table('Sample: protein\nMass\tIntensity\tCharge\n1000.123456789\t900\t3\n701.2\t40\t0\n800.1\t50\t\n')
        self.assertEqual(data['peaks'][0]['mz'], 1000.123456789)
        self.assertEqual(data['peaks'][0]['charge'], 3)
        self.assertEqual(data['peaks'][0]['intensity'], 900)
        self.assertEqual(data['unknown_charge_count'], 2)
        self.assertIsNone(data['peaks'][2]['charge'])

    def test_signed_charge_and_scan(self):
        data = parse_peak_table('m/z,Z,Scan\n500.25,3+,9\n600.4,2−,10\n700.1,-4,11\n')
        self.assertEqual(data['polarities'], ['negative', 'positive'])
        self.assertEqual([p['charge'] for p in data['peaks']], [3, 2, 4])
        self.assertEqual(data['peaks'][2]['scan'], 11)

    def test_text_delimiters_encodings_and_quoted_csv(self):
        for body in ['m/z;Z\n500.2;2', 'm/z Z\n500.2 2', '"m/z","Charge"\n"500.2","2"']:
            for encoding in ['utf-8-sig', 'utf-16']:
                self.assertEqual(import_peak_table(body.encode(encoding))['peaks'][0]['mz'], 500.2)

    def test_intensity_is_never_treated_as_charge(self):
        for text in ['500.2,2', 'Mass,Intensity\n500.2,2', 'm/z,Charge\n']:
            with self.assertRaises(ValueError):
                parse_peak_table(text)

    def test_reject_invalid_values_without_partial_import(self):
        for row in ['nan,2', 'inf,2', '-50,2', '500,2.5', '500,101', '500,+2-', '500,2,extra']:
            with self.subTest(row=row), self.assertRaises(ValueError):
                parse_peak_table('m/z,Charge\n600,3\n' + row)
        for text in ['m/z,Charge,Intensity\n500,2,inf', 'Mass,Charge,Scan\n500,2,-1', 'm/z,Mass,Charge\n500,500,2']:
            with self.assertRaises(ValueError):
                parse_peak_table(text)

    def test_size_row_and_encoding_limits(self):
        for body in [b'', b'x' * 2_000_001, b'\x80\x81', b'm/z,Charge\n' + b'1' * 140_000 + b',2']:
            with self.assertRaises(ValueError):
                import_peak_table(body)
        with self.assertRaises(ValueError):
            parse_peak_table('m/z,Charge\n' + '500,2\n' * 10_001)

    def test_import_route_and_failures(self):
        for body, expected in [(b'Mass,Charge\n500.125,2', 200), (b'bad table', 400), (b'', 400)]:
            h = object.__new__(Handler)
            h.path = '/api/import-peaks'
            h.headers = {'Content-Length': str(len(body))}
            h.rfile = io.BytesIO(body)
            result = []
            h.send = lambda status, data: result.append((status, data))
            h.do_POST()
            self.assertEqual(result[0][0], expected)
            json.dumps(result[0][1], allow_nan=False)


if __name__ == '__main__':
    unittest.main()
