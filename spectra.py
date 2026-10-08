"""Validate exported peak tables without inferring unrecorded charge states."""
import csv
import math
import re

MAX_EXPORT_BYTES = 2_000_000
MAX_PEAKS = 10_000


def parse_peak_table(text):
    """Read an explicit m/z + charge table (CSV, TSV, or whitespace separated).

    Qual Browser spectrum lists may have metadata above the table. Only a header
    naming m/z (or Mass) and charge is accepted; a two-column mass/intensity list
    must never be mistaken for m/z/charge. Unknown charges remain unknown.
    """
    lines = text.lstrip('\ufeff').splitlines()
    aliases = {
        'mz': 'mz', 'mass': 'mz', 'masscharge': 'mz', 'masstocharge': 'mz',
        'charge': 'charge', 'z': 'charge', 'chargestate': 'charge',
        'intensity': 'intensity', 'abundance': 'intensity', 'relativeabundance': 'intensity',
        'scan': 'scan', 'scannumber': 'scan', 'scanno': 'scan',
    }

    def columns(line, separator):
        return next(csv.reader([line], delimiter=separator)) if separator else re.split(r'\s+', line.strip())

    header = None
    for index, line in enumerate(lines[:100]):
        for separator in ('\t', ',', ';', None):
            names = [aliases.get(re.sub(r'[^a-z0-9]', '', item.lower())) for item in columns(line, separator)]
            if 'mz' in names and 'charge' in names:
                known = [name for name in names if name]
                if len(set(known)) != len(known):
                    raise ValueError('The peak table contains duplicate column names.')
                header = (index, separator, names)
                break
        if header:
            break
    if header is None:
        raise ValueError('Include a header with m/z (or Mass) and Charge (or Z). A mass/intensity list alone does not supply charge values.')

    index, separator, names = header
    peaks = []
    polarities = set()
    for line_number, line in enumerate(lines[index + 1:], index + 2):
        if not line.strip():
            continue
        cells = columns(line, separator)
        if len(cells) != len(names):
            raise ValueError(f'Line {line_number}: expected {len(names)} columns. Export a single spectrum table.')
        row = {name: cell.strip() for name, cell in zip(names, cells) if name}
        try:
            mz = float(row['mz'])
            if not math.isfinite(mz) or mz <= 0:
                raise ValueError()
        except ValueError:
            raise ValueError(f'Line {line_number}: m/z must be a finite positive number.') from None
        value = row['charge'].replace('−', '-')
        charge = None
        polarity = None
        if value.lower() not in ('', '0', '-', '?', 'n/a', 'na', 'unknown', 'unassigned'):
            match = re.fullmatch(r'([+-]?)(\d+)([+-]?)', value)
            if not match or (match[1] and match[3]):
                raise ValueError(f'Line {line_number}: charge must be an integer, such as 3, 3+, or 3−.')
            charge = int(match[2])
            if not 1 <= charge <= 100:
                raise ValueError(f'Line {line_number}: charge magnitude must be between 1 and 100.')
            sign = match[1] or match[3]
            if sign:
                polarity = 'negative' if sign == '-' else 'positive'
                polarities.add(polarity)
        intensity = None
        if row.get('intensity'):
            try:
                intensity = float(row['intensity'])
                if not math.isfinite(intensity) or intensity < 0:
                    raise ValueError()
            except ValueError:
                raise ValueError(f'Line {line_number}: intensity must be a finite nonnegative number.') from None
        scan = row.get('scan', '')
        if scan and (not scan.isdigit() or int(scan) < 1):
            raise ValueError(f'Line {line_number}: scan must be a positive integer.')
        peaks.append(dict(mz=mz, charge=charge, intensity=intensity,
                          scan=int(scan) if scan else None, polarity=polarity))
        if len(peaks) > MAX_PEAKS:
            raise ValueError(f'Import at most {MAX_PEAKS:,} peaks at once. Export a smaller spectrum table.')
    if not peaks:
        raise ValueError('The peak table has a header but no peaks.')
    return dict(peaks=peaks, unknown_charge_count=sum(p['charge'] is None for p in peaks),
                polarities=sorted(polarities), source='Imported peak table; values supplied by the file')


def import_peak_table(body):
    if not body or len(body) > MAX_EXPORT_BYTES:
        raise ValueError('Peak tables must be between 1 byte and 2 MB.')
    try:
        # Excel and Windows exports may use UTF-16 with a byte-order mark.
        text = body.decode('utf-16' if body.startswith((b'\xff\xfe', b'\xfe\xff')) else 'utf-8-sig')
    except UnicodeError:
        raise ValueError('Save the peak table as UTF-8 or UTF-16 text (CSV, TSV, or TXT).') from None
    try:
        return parse_peak_table(text)
    except csv.Error:
        raise ValueError('The peak table contains an invalid or excessively long CSV field.') from None
