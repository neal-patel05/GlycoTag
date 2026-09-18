"""Import UniProt annotations without inferring occupancy from sequence motifs."""
import json
import re
from urllib.request import Request, urlopen


def parse_record(record):
    seq = record['sequence']['value']
    accession = record['primaryAccession']
    sites, skipped = {}, []
    for feature in record.get('features', []):
        if feature.get('type') != 'Glycosylation':
            continue
        description = feature.get('description','')
        match = re.match(r'([NO])-linked\b',description,flags=re.I)
        location = feature.get('location',{})
        start,end = location.get('start',{}),location.get('end',{})
        pos = start.get('value')
        reason = None
        if not match:
            reason = 'Not a supported N- or O-linked annotation'
        elif location.get('sequence') not in (None,accession):
            reason = 'Annotation uses a different isoform sequence'
        elif (not isinstance(pos,int) or pos != end.get('value') or
              start.get('modifier','EXACT') != 'EXACT' or end.get('modifier','EXACT') != 'EXACT'):
            reason = 'Position is uncertain or spans a range'
        elif not 1 <= pos <= len(seq):
            reason = 'Position is outside this sequence'
        elif seq[pos-1] not in ('N' if match[1].upper() == 'N' else 'ST'):
            reason = 'Residue chemistry is outside this tool’s N / S / T support'
        if reason:
            skipped.append(dict(description=description,reason=reason,position=pos))
            continue
        kind = match[1].upper()
        site = sites.setdefault((pos,kind),dict(position=pos,type=kind,residue=seq[pos-1],descriptions=[],evidences=[]))
        if description not in site['descriptions']:
            site['descriptions'].append(description)
        for evidence in feature.get('evidences',[]):
            if evidence not in site['evidences']:
                site['evidences'].append(evidence)
    sites = [sites[key] for key in sorted(sites)]
    names = record.get('proteinDescription',{})
    name = names.get('recommendedName',{}).get('fullName',{}).get('value','')
    if not name:
        name = next((n.get('fullName',{}).get('value','') for n in names.get('submissionNames',[])),accession)
    return dict(sequence=seq,accession=accession,name=name,sites=sites,skipped=skipped,
                sites_text='\n'.join(f"{s['position']},{s['type']}" for s in sites),
                source=f'https://www.uniprot.org/uniprotkb/{accession}/entry',
                annotation_date=record.get('entryAudit',{}).get('lastAnnotationUpdateDate'),
                warnings=[])


def download(path):
    with urlopen(Request('https://rest.uniprot.org/uniprotkb/'+path,
                         headers={'User-Agent':'GlycoTag/1.1'}),timeout=20) as response:
        body = response.read(5000001)
    if len(body) > 5000000:
        raise ValueError('UniProt response is too large. Paste the sequence and sites manually.')
    return body.decode('utf-8')


def import_entry(identifier):
    identifier = identifier.strip().upper()
    if not re.fullmatch(r'(?:[A-Z0-9]{6,10}(?:-\d+)?|[A-Z0-9]{1,16}_[A-Z0-9]{1,16})',identifier):
        raise ValueError('Enter a UniProt accession or entry ID, e.g. P02787 or TRFE_HUMAN.')
    base = identifier.split('-')[0]
    result = parse_record(json.loads(download(base+'.json')))
    if '-' in identifier:
        fasta = download(identifier+'.fasta')
        seq = ''.join(line.strip() for line in fasta.splitlines() if not line.startswith('>'))
        if not seq:
            raise ValueError('UniProt returned an empty isoform sequence.')
        if seq != result['sequence']:
            result.update(sequence=seq,sites=[],sites_text='',skipped=[])
            result['warnings'].append('Isoform sequence differs from the annotated canonical sequence. Sites were not transferred; enter isoform-specific positions manually.')
        result['accession'] = identifier
    return result
