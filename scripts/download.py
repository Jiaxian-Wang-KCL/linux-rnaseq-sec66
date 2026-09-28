#!/usr/bin/env python3
"""Fetch complete ENA FASTQs and matched Ensembl 115 reference files."""
import csv
import gzip
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def digest(path, algorithm):
    h = hashlib.new(algorithm)
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def download(url, path, md5=None, size=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and md5:
        if digest(path, 'md5') != md5:
            raise RuntimeError(f'Checksum mismatch in existing {path}; inspect before replacing')
    if not path.exists():
        temp = path.with_suffix(path.suffix + '.partial')
        subprocess.run(['curl', '-fL', '--retry', '4', '--connect-timeout', '30',
                        '--continue-at', '-', '--output', str(temp), url], check=True)
        if size is not None and temp.stat().st_size != size:
            raise RuntimeError(f'Unexpected size: {temp}')
        if md5 and digest(temp, 'md5') != md5:
            raise RuntimeError(f'Checksum mismatch: {temp}')
        temp.rename(path)
    if size is not None and path.stat().st_size != size:
        raise RuntimeError(f'Unexpected size: {path}')
    return dict(url=url, path=str(path.relative_to(ROOT)), bytes=path.stat().st_size,
                sha256=digest(path, 'sha256'), expected_ena_md5=md5)


records = []
with (ROOT / 'metadata/samples.tsv').open() as f:
    samples = list(csv.DictReader(f, delimiter='\t'))
for s in samples:
    assert s['library_layout'] == 'SINGLE' and ';' not in s['fastq_ftp']
    print('Downloading complete sample:', s['sample'], s['run_accession'], flush=True)
    records.append(download('https://' + s['fastq_ftp'],
                            ROOT / 'data/raw' / (s['run_accession'] + '.fastq.gz'),
                            s['fastq_md5'], int(s['fastq_bytes'])))

base = 'https://ftp.ensembl.org/pub/release-115'
references = [
    (f'{base}/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz', 'genome.fa'),
    (f'{base}/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.gtf.gz', 'genes.gtf'),
]
for url, name in references:
    gz = ROOT / 'reference' / (name + '.gz')
    records.append(download(url, gz))
    out = ROOT / 'reference' / name
    if not out.exists():
        with gzip.open(gz, 'rb') as src, out.open('wb') as dst:
            while chunk := src.read(1024 * 1024):
                dst.write(chunk)
(ROOT / 'metadata/downloads.json').write_text(json.dumps(records, indent=2) + '\n')
print('All downloads verified. FASTQs were not subsampled.', flush=True)
