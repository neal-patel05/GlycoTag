import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch

import glycosmos
from app import Handler


class HostedImageTests(unittest.TestCase):
    def tearDown(self):
        glycosmos.image_ids.cache_clear()

    def test_manifest_preserves_availability_without_large_archive(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(glycosmos, 'IMAGE_DATABASE', Path(directory)/'absent.sqlite3'):
            glycosmos.image_ids.cache_clear()
            self.assertEqual(len(glycosmos.image_ids()), 192581)
            self.assertTrue(glycosmos.record('G00039MO')['has_image'])
            self.assertFalse(glycosmos.record('G00000CV')['has_image'])
            self.assertIsNone(glycosmos.structure_image('G00039MO'))
            self.assertEqual(glycosmos.hosted_image_url('G00039MO'), 'https://image.glycosmos.org/snfg/png/G00039MO')
            self.assertIsNone(glycosmos.hosted_image_url('G00000CV'))
            with self.assertRaises(ValueError):
                glycosmos.hosted_image_url('../anything')

    def test_hosted_route_redirects_only_known_image_ids(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(glycosmos, 'IMAGE_DATABASE', Path(directory)/'absent.sqlite3'):
            glycosmos.image_ids.cache_clear()
            handler = object.__new__(Handler)
            handler.path = '/api/glycosmos/image/G00039MO'
            with patch.object(handler, 'send_response') as status, patch.object(handler, 'send_header') as header, patch.object(handler, 'end_headers'):
                handler.do_GET()
                status.assert_called_once_with(302)
                self.assertIn(('Location', 'https://image.glycosmos.org/snfg/png/G00039MO'), [call.args for call in header.call_args_list])
            handler.path = '/api/glycosmos/image/G00000CV'
            with patch.object(handler, 'send') as response:
                handler.do_GET()
                self.assertEqual(response.call_args.args[0], 404)

    def test_local_archive_continues_to_serve_png(self):
        self.assertIsNone(glycosmos.hosted_image_url('G00039MO'))
        self.assertTrue(glycosmos.structure_image('G00039MO').startswith(b'\x89PNG'))
