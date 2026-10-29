#!/usr/bin/env python3
"""Import the public https://image.glycosmos.org/snfg-png.zip archive locally.

Usage: python3 scripts/build_glycosmos_images.py /path/to/snfg-png.zip
The large local image cache is not a deployable asset or a Git artifact.
"""
import argparse
import hashlib
import gzip
import json
from pathlib import Path
import re
import sqlite3
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def build(source):
    destination = ROOT / 'data/glycosmos-images.sqlite3'
    temporary = destination.with_suffix('.building')
    temporary.unlink(missing_ok=True)
    count = 0
    with sqlite3.connect(temporary) as db, zipfile.ZipFile(source) as archive:
        db.execute('CREATE TABLE images (id TEXT PRIMARY KEY, png BLOB NOT NULL)')
        for item in archive.infolist():
            match = re.fullmatch(r'(G\d{5}[A-Z]{2})(?:\.png)?', Path(item.filename).name)
            if not match:
                continue
            png = archive.read(item)
            if not png.startswith(b'\x89PNG\r\n\x1a\n'):
                raise ValueError('Invalid PNG: ' + item.filename)
            db.execute('INSERT OR IGNORE INTO images VALUES (?,?)', (match.group(1), png))
            count += 1
    if not count:
        raise ValueError('Archive contains no glycan PNG images.')
    temporary.replace(destination)
    with sqlite3.connect(destination) as db:
        ids = [row[0] for row in db.execute('SELECT id FROM images ORDER BY id')]
    # Small portable manifest preserves image availability on hosted installs.
    manifest = destination.with_name('glycosmos-image-ids.json.gz')
    manifest.write_bytes(gzip.compress(json.dumps(ids, separators=(',', ':')).encode(), mtime=0))
    digest = hashlib.sha256()
    with source.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024), b''):
            digest.update(chunk)
    metadata = dict(source='https://image.glycosmos.org/snfg-png.zip', sha256=digest.hexdigest(),
                    images=count, attribution='GlyCosmos Portal', license='CC BY 4.0')
    destination.with_suffix('.metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(json.dumps(metadata))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('zip', type=Path)
    build(parser.parse_args().zip)
