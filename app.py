#!/usr/bin/env python3
"""Run locally with python3 app.py. Python 3.9+, no packages required."""
import argparse
import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.error import URLError
from engine import analyze, EXAMPLE
from glycans import catalog_payload
from uniprot import import_entry
from spectra import import_peak_table, MAX_EXPORT_BYTES

ROOT = Path(__file__).resolve().parent
class Handler(BaseHTTPRequestHandler):
    def send(self, status, body, content_type='application/json'):
        raw = json.dumps(body,allow_nan=False).encode() if content_type == 'application/json' else body
        self.send_response(status)
        self.send_header('Content-Type',content_type + ('; charset=utf-8' if content_type.startswith('text/') or content_type == 'application/json' else ''))
        self.send_header('Content-Length',str(len(raw)))
        self.send_header('X-Content-Type-Options','nosniff')
        self.end_headers()
        self.wfile.write(raw)

    def do_GET(self):
        if self.path == '/api/example':
            self.send(200,dict(sequence=EXAMPLE))
        elif self.path == '/api/glycans':
            self.send(200,catalog_payload())
        elif self.path.startswith('/api/uniprot/'):
            try:
                self.send(200,import_entry(self.path.removeprefix('/api/uniprot/')))
            except (URLError,TimeoutError,OSError):
                self.send(502,dict(error='UniProt could not be reached or accession was not found. Paste the FASTA sequence instead.'))
            except (ValueError,KeyError,TypeError):
                self.send(400,dict(error='Check the UniProt accession / entry ID, or paste the sequence and sites manually.'))
        elif self.path in ('/','/style.css','/main.js','/glypeptid-logo.png','/spectra.js'):
            file,mime = {'/':('index.html','text/html'),'/style.css':('style.css','text/css'),'/main.js':('main.js','text/javascript'),'/glypeptid-logo.png':('glypeptid-logo.png','image/png'),'/spectra.js':('spectra.js','text/javascript')}[self.path]
            self.send(200,(ROOT/'static'/file).read_bytes(),mime)
        else:
            self.send(404,dict(error='Not found'))

    def do_POST(self):
        if self.path == '/api/import-peaks':
            try:
                length = int(self.headers.get('Content-Length', 0))
                if not 0 < length <= MAX_EXPORT_BYTES:
                    raise ValueError('Peak tables must be between 1 byte and 2 MB.')
                self.send(200, import_peak_table(self.rfile.read(length)))
            except ValueError as exc:
                self.send(400, dict(error=str(exc)))
            return
        if self.path != '/api/analyze':
            return self.send(404,dict(error='Not found'))
        try:
            length = int(self.headers.get('Content-Length',0))
            if not 0 < length <= 200000:
                raise ValueError('Request size must be between 1 and 200,000 bytes.')
            data = json.loads(self.rfile.read(length))
            if not isinstance(data,dict):
                raise ValueError('Expected a JSON object.')
            self.send(200,analyze(data))
        except (ValueError,TypeError,KeyError) as exc:
            self.send(400,dict(error=str(exc)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port',type=int,default=8000)
    args = parser.parse_args()
    server = ThreadingHTTPServer(('127.0.0.1',args.port),Handler)
    print(f'GlyPeptID is running at http://localhost:{args.port}',flush=True)
    server.serve_forever()
