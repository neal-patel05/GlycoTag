#!/usr/bin/env python3
"""Run locally with python3 app.py. Python 3.9+, no packages required."""
import argparse
import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.error import URLError
from urllib.parse import parse_qs, urlsplit
from engine import analyze, EXAMPLE
from glycans import catalog_payload
from uniprot import import_entry
from adducts import catalog_payload as adduct_payload
from glycosmos import browse, structure_image, hosted_image_url, record
from glycan_adducts import analyze_glycan_adducts

ROOT = Path(__file__).resolve().parent
class Handler(BaseHTTPRequestHandler):
    def send(self, status, body, content_type='application/json'):
        raw = json.dumps(body,allow_nan=False).encode() if content_type == 'application/json' else body
        self.send_response(status)
        self.send_header('Content-Type',content_type + ('; charset=utf-8' if content_type.startswith('text/') or content_type == 'application/json' else ''))
        self.send_header('Content-Length',str(len(raw)))
        self.send_header('X-Content-Type-Options','nosniff')
        self.send_header('Cache-Control', 'public, max-age=86400' if content_type == 'image/png' else 'no-store')
        self.end_headers()
        self.wfile.write(raw)

    def do_GET(self):
        url = urlsplit(self.path)
        if url.path.startswith('/api/glycosmos/image/'):
            try:
                identifier = url.path.removeprefix('/api/glycosmos/image/')
                png = structure_image(identifier)
                remote = hosted_image_url(identifier) if png is None else None
                if remote:
                    self.send_response(302)
                    self.send_header('Location', remote)
                    self.send_header('Cache-Control', 'public, max-age=86400')
                    self.send_header('Content-Length', '0')
                    self.end_headers()
                elif png is None:
                    self.send(404,dict(error='No structure image supplied in the GlyCosmos image archive.'))
                else:
                    self.send(200,png,'image/png')
            except ValueError as exc:
                self.send(400,dict(error=str(exc)))
        elif url.path.startswith('/api/glycosmos/record/'):
            try:
                glycan = record(url.path.removeprefix('/api/glycosmos/record/'))
                self.send(200 if glycan else 404, glycan or dict(error='GlyCosmos record not found.'))
            except ValueError as exc:
                self.send(400, dict(error=str(exc)))
        elif url.path == '/api/glycosmos':
            try:
                self.send(200, browse(parse_qs(url.query)))
            except ValueError as exc:
                self.send(400, dict(error=str(exc)))
        elif url.path == '/api/adducts':
            self.send(200, adduct_payload())
        elif url.path == '/api/example':
            self.send(200,dict(sequence=EXAMPLE))
        elif url.path == '/api/glycans':
            self.send(200,catalog_payload())
        elif url.path.startswith('/api/uniprot/'):
            try:
                self.send(200,import_entry(url.path.removeprefix('/api/uniprot/')))
            except (URLError,TimeoutError,OSError):
                self.send(502,dict(error='UniProt could not be reached or accession was not found. Paste the FASTA sequence instead.'))
            except (ValueError,KeyError,TypeError):
                self.send(400,dict(error='Check the UniProt accession / entry ID, or paste the sequence and sites manually.'))
        elif url.path in ('/','/style.css','/main.js','/glypeptid-logo.png'):
            file,mime = {'/':('index.html','text/html'),'/style.css':('style.css','text/css'),'/main.js':('main.js','text/javascript'),'/glypeptid-logo.png':('glypeptid-logo.png','image/png')}[url.path]
            self.send(200,(ROOT/'static'/file).read_bytes(),mime)
        else:
            self.send(404,dict(error='Not found'))

    def do_POST(self):
        path = urlsplit(self.path).path
        if path not in ('/api/analyze', '/api/glycan-adducts'):
            return self.send(404,dict(error='Not found'))
        try:
            length = int(self.headers.get('Content-Length',0))
            if not 0 < length <= 200000:
                raise ValueError('Request size must be between 1 and 200,000 bytes.')
            data = json.loads(self.rfile.read(length))
            if not isinstance(data,dict):
                raise ValueError('Expected a JSON object.')
            analyze_request = analyze_glycan_adducts if path == '/api/glycan-adducts' else analyze
            self.send(200,analyze_request(data))
        except (ValueError,TypeError,KeyError) as exc:
            self.send(400,dict(error=str(exc)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port',type=int,default=8000)
    args = parser.parse_args()
    server = ThreadingHTTPServer(('127.0.0.1',args.port),Handler)
    print(f'GlyPeptID is running at http://localhost:{args.port}',flush=True)
    server.serve_forever()
