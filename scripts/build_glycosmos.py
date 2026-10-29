#!/usr/bin/env python3
"""Build the full local library from GlyCosmos's public glycan CSV.

python3 scripts/build_glycosmos.py /path/to/glycosmos_glycans_list.csv
Download source: https://download.glycosmos.org/glycosmos_glycans_list.csv
No network calls are made by this script. The completed database is replaced
atomically only after parsing and validating the entire input.
"""
import argparse
from collections import Counter
import csv
from datetime import date
import hashlib
import gzip
import json
from pathlib import Path
import re
import sqlite3
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from glycosmos import DATABASE, RESIDUES, composition_mass


def build(source, output=DATABASE, retrieved=None):
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix('.building')
    temporary.unlink(missing_ok=True)
    db = sqlite3.connect(temporary)
    reasons = Counter()
    total = supported = annotated = 0
    taxa = {}
    try:
        db.execute('CREATE TABLE glycans (id TEXT PRIMARY KEY, iupac TEXT, wurcs TEXT, motifs TEXT, composition TEXT, free_mass REAL, attached_mass REAL, reported_mass REAL, mass_note TEXT, taxonomy TEXT)')
        db.execute('CREATE TABLE glycan_taxa (glycan_id TEXT, taxon_id TEXT, PRIMARY KEY (glycan_id,taxon_id))')
        with source.open(encoding='utf-8-sig', newline='') as stream:
            reader = csv.DictReader(stream)
            if not {'Accession Number', 'WURCS', 'IUPAC Condensed', 'Monoisotopic mass', 'Motif name(s)', 'Taxonomy'}.issubset(reader.fieldnames or []):
                raise ValueError('Expected the complete GlyCosmos glycans CSV.')
            for row in reader:
                identifier = row['Accession Number']
                if not re.fullmatch(r'G\d{5}[A-Z]{2}', identifier) or None in row or any(v is None for v in row.values()):
                    raise ValueError('Malformed or incomplete CSV row.')
                reported = float(row['Monoisotopic mass']) if row['Monoisotopic mass'] else None
                try:
                    masses = composition_mass(row['WURCS'], reported)
                    supported += 1
                    note = 'ExPASy GlycanMass composition calculation; underivatized, monoisotopic, free reducing end.'
                except ValueError as error:
                    masses = dict(composition={}, free_mass=None, attached_mass=None)
                    note = str(error)
                    reasons[note] += 1
                db.execute('INSERT INTO glycans VALUES (?,?,?,?,?,?,?,?,?,?)', (
                    identifier, row['IUPAC Condensed'], row['WURCS'], row['Motif name(s)'],
                    json.dumps(masses['composition'], sort_keys=True), masses['free_mass'],
                    masses['attached_mass'], reported, note, row['Taxonomy']))
                if row['Taxonomy']:
                    annotated += 1
                for entry in filter(None, row['Taxonomy'].split(',')):
                    name, separator, taxon_id = entry.rpartition(':')
                    if not separator or not name.strip() or not taxon_id.isdigit():
                        raise ValueError(f'Malformed taxonomy annotation for {identifier}.')
                    if taxon_id in taxa and taxa[taxon_id] != name.strip():
                        raise ValueError(f'Conflicting taxonomy names for {taxon_id}.')
                    taxa[taxon_id] = name.strip()
                    db.execute('INSERT OR IGNORE INTO glycan_taxa VALUES (?,?)', (identifier, taxon_id))
                total += 1
        if not total:
            raise ValueError('Empty catalogue.')
        metadata = dict(source='https://download.glycosmos.org/glycosmos_glycans_list.csv',
                        license='CC BY 4.0', attribution='GlyCosmos Portal',
                        retrieved=retrieved or str(date.today()), sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                        total=total, supported=supported, unavailable=total-supported,
                        taxonomy_annotated=annotated, taxonomy_unannotated=total-annotated,
                        taxonomy_count=len(taxa), taxonomy_source_field='Taxonomy',
                        mass_source='https://web.expasy.org/glycanmass/',
                        residues={k:str(v) for k,v in RESIDUES.items()}, unavailable_reasons=dict(reasons))
        db.execute('CREATE TABLE taxa (id TEXT PRIMARY KEY, name TEXT NOT NULL)')
        db.executemany('INSERT INTO taxa VALUES (?,?)', sorted(taxa.items()))
        db.execute('CREATE INDEX taxa_to_glycans ON glycan_taxa(taxon_id,glycan_id)')
        db.execute('CREATE TABLE metadata (value TEXT)')
        db.execute('INSERT INTO metadata VALUES (?)', (json.dumps(metadata),))
        db.commit()
    finally:
        db.close()
    temporary.replace(output)
    # Derive display labels independently of whether ExPASy supports a mass.
    with sqlite3.connect(output) as source_db:
        monomers = source_db.execute("SELECT id,iupac,wurcs FROM glycans WHERE wurcs LIKE 'WURCS=2.0/1,1,0/%' AND iupac!='' ORDER BY id")
        labels = {}
        for accession, name, wurcs in monomers:
            unit = re.search(r'\[([^]]+)\]', wurcs).group(1)
            labels[unit] = dict(label=re.sub(r'\([ab?]\d.*$', '', name), source_id=accession)
    (output.parent/'residue-units.json').write_text(json.dumps(labels, separators=(',', ':'))+'\n')
    # A compressed snapshot keeps the full catalogue small enough to share.
    # Runtime extracts it to the OS temporary directory, never the static site.
    snapshot = output.with_suffix('.sqlite3.gz')
    packed = snapshot.with_suffix('.building')
    with packed.open('wb') as stream, gzip.GzipFile(filename='', mode='wb', fileobj=stream, mtime=0) as target:
        target.write(output.read_bytes())
    packed.replace(snapshot)
    output.with_suffix('.metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('csv', type=Path)
    parser.add_argument('--retrieved', help='Actual source download date (YYYY-MM-DD) when rebuilding an existing download.')
    args = parser.parse_args()
    build(args.csv, retrieved=args.retrieved)
