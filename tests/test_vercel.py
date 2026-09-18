import io
import json
import unittest
from unittest.mock import patch
from api.index import handler
from engine import PROTON, MONO

class DeploymentTests(unittest.TestCase):
    def request(self,path,data=None):
        # Exercise the actual Vercel handler without needing a listening socket.
        h=object.__new__(handler)
        h.path=path
        body=json.dumps(data).encode() if data is not None else b''
        h.rfile=io.BytesIO(body)
        h.headers={'Content-Length':str(len(body))}
        response=[]
        h.send=lambda status,body,content_type='application/json':response.append((status,body))
        h.do_GET() if data is None else h.do_POST()
        return response[0]

    def test_rewritten_get_routes(self):
        code,body=self.request('/api?route=glycans')
        self.assertEqual(code,200)
        self.assertEqual(len(body['glycans']),21)
        code,body=self.request('/api?route=example')
        self.assertEqual(code,200)
        self.assertEqual(len(body['sequence']),698)

    def test_rewritten_uniprot_route(self):
        with patch('app.import_entry',return_value={'sequence':'S','sites_text':'1,O'}) as get:
            code,body=self.request('/api?route=uniprot%2FP02787')
            self.assertEqual(code,200)
            get.assert_called_once_with('P02787')

    def test_rewritten_post_keeps_precision_and_proton_convention(self):
        mz=(MONO['S']+674.2382)/2-PROTON
        code,body=self.request('/api?route=analyze',dict(sequence='S',sites='1,O',
            ion_mode='negative',glycans_O='674.2382',observations=f'{mz},2'))
        self.assertEqual(code,200)
        row=body['results'][0]['rows'][0]
        self.assertAlmostEqual(row['error_da'],0)
        self.assertEqual(row['breakdown']['amino_acids'][0]['residue_mass_text'],'87.03202840472')

    def test_invalid_requests_and_direct_public_path(self):
        self.assertEqual(self.request('/api?route=unknown')[0],404)
        self.assertEqual(self.request('/api?route=analyze',[])[0],400)
        self.assertEqual(self.request('/api/glycans')[0],200)

    def test_hosted_response_limit_returns_clear_error(self):
        h=object.__new__(handler)
        with patch('app.Handler.send') as send:
            h.send(200,{'large':'x'*4000001})
            self.assertEqual(send.call_args.args[0],413)
            self.assertIn('fewer measurements',send.call_args.args[1]['error'])

if __name__=='__main__':unittest.main()
