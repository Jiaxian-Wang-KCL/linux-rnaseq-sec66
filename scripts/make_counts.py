#!/usr/bin/env python3
"""Check strand direction across all replicates and extract raw integer counts."""
import csv
import json
import re
import shutil
from pathlib import Path

p = Path(__file__).resolve().parents[1]
samples = list(csv.DictReader((p/'metadata/samples.tsv').open(), delimiter='\t'))
names = [s['sample'] for s in samples]
assignment = {}
for mode in (0, 1, 2):
    with (p/f'results/strand_check/strand{mode}.tsv.summary').open() as f:
        rows = list(csv.reader(f, delimiter='\t'))
    assert [Path(x).stem for x in rows[0][1:]] == names
    assignment[mode] = list(map(int, next(r[1:] for r in rows if r[0] == 'Assigned')))
best = 1 if sum(assignment[1]) > sum(assignment[2]) else 2
other = 3 - best
checks = []
for i, name in enumerate(names):
    a, b, u = assignment[best][i], assignment[other][i], assignment[0][i]
    record = dict(sample=name, unstranded=u, forward=assignment[1][i], reverse=assignment[2][i],
                  selected_fraction_of_unstranded=a/max(u, 1), selected_over_opposite=a/max(b, 1))
    checks.append(record)
    if not (a/max(u, 1) > 0.8 and a/max(b, 1) > 4):
        raise RuntimeError(f'Library direction not clearly consistent across all samples: {record}')
(p/'metadata/strandedness.json').write_text(json.dumps(dict(
    selected_featureCounts_s=best, rationale='All six samples: selected/unstranded >0.8; selected/opposite >4.',
    samples=checks), indent=2)+'\n')
src = p/f'results/strand_check/strand{best}.tsv'
shutil.copyfile(src, p/'results/counts/featureCounts.tsv')
shutil.copyfile(str(src)+'.summary', p/'results/counts/featureCounts.tsv.summary')
with src.open() as f:
    rows = list(csv.reader((line for line in f if not line.startswith('#')), delimiter='\t'))
assert [Path(x).stem for x in rows[0][6:]] == names
seen = set()
with (p/'results/counts/count_matrix.tsv').open('w') as f:
    w = csv.writer(f, delimiter='\t'); w.writerow(['gene_id']+names)
    for row in rows[1:]:
        assert row[0] not in seen
        seen.add(row[0])
        counts = [int(v) for v in row[6:]]
        assert len(counts)==6 and all(x>=0 for x in counts)
        w.writerow([row[0]]+counts)

genes = {}
for line in (p/'reference/genes.gtf').open():
    if line.startswith('#'): continue
    fields = line.rstrip('\n').split('\t')
    attrs = dict(re.findall(r'(\w+) "([^"]+)"', fields[8]))
    gid = attrs.get('gene_id')
    if gid: genes[gid] = (attrs.get('gene_name', gid), attrs.get('gene_biotype', ''))
with (p/'metadata/gene_annotation.tsv').open('w') as f:
    w = csv.writer(f, delimiter='\t'); w.writerow(['gene_id','gene_name','biotype'])
    w.writerows((g, *v) for g,v in sorted(genes.items()))

with (p/'results/sample_qc.tsv').open('w') as f:
    w=csv.writer(f, delimiter='\t')
    w.writerow(['sample','condition','raw_reads','retained_reads','retained_percent',
                'alignment_percent','assigned_reads','assigned_percent_of_retained','strand'])
    for i,s in enumerate(samples):
        q=json.loads((p/f'results/fastp/{s["sample"]}.json').read_text())['summary']
        raw=q['before_filtering']['total_reads']; retained=q['after_filtering']['total_reads']
        assert raw == int(s['read_count']), 'Input FASTQ count differs from ENA record'
        summary=(p/f'logs/{s["sample"]}.hisat2.summary.txt').read_text()
        rate=float(re.search(r'([\d.]+)% overall alignment rate',summary).group(1))
        if rate < 70: raise RuntimeError(f'Low alignment rate for {s["sample"]}: {rate}%')
        assigned=assignment[best][i]
        if assigned/retained < 0.5: raise RuntimeError('Low assigned fraction; inspect before DE analysis')
        w.writerow([s['sample'],s['condition'],raw,retained,100*retained/raw,rate,
                    assigned,100*assigned/retained,best])
print(f'Validated {len(seen)} genes × {len(names)} samples; selected strand={best}')
