"""Vercel entrypoint; shares the tested local application's API handlers."""
from urllib.parse import parse_qs, urlsplit
import json
from app import Handler


class handler(Handler):
    def send(self, status, body, content_type='application/json'):
        if content_type == 'application/json' and len(json.dumps(body,allow_nan=False).encode()) > 4000000:
            status,body = 413,dict(error='This result is too large for hosting. Submit fewer measurements per request.')
        return super().send(status,body,content_type)

    def route(self):
        url = urlsplit(self.path)
        # vercel.json passes the public API suffix explicitly through rewrites.
        route = parse_qs(url.query).get('route', [None])[0]
        self.path = '/api/' + route if route is not None else url.path

    def do_GET(self):
        self.route()
        return super().do_GET()

    def do_POST(self):
        self.route()
        return super().do_POST()
