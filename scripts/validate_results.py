#!/usr/bin/env python3
"""Cross-check real output tables and record every FastQC module status."""
import csv
import json
import zipfile
from pathlib import Path

p=Path(__file__).resolve().parents[1]
def read(rel):
    with (p/rel).open() as f:
        return list(csv.DictReader(f,delimiter='\t'))

meta=read('metadata/samples.tsv')
qc=read('results/sample_qc.tsv')
matrix=read('results/counts/count_matrix.tsv')
de=read('results/differential_expression/all_tested_genes.tsv')
sig=read('results/differential_expression/DEGs_FDR05_absLFC1.tsv')
stats={r['metric']:float(r['value']) for r in read('results/analysis_summary.tsv')}
assert [x['sample'] for x in meta]==[x['sample'] for x in qc]
assert len(meta)==6 and len(matrix)==stats['genes_in_matrix']
assert len(de)==stats['genes_after_low_count_filter']
assert len(sig)==stats['FDR05_absLFC1_genes']
assert len({r['gene_id'] for r in matrix})==len(matrix)
for s,q in zip(meta,qc):
    assert sum(int(r[s['sample']]) for r in matrix)==int(q['assigned_reads'])
    assert int(s['read_count'])==int(q['raw_reads'])
    assert int(q['retained_reads'])<=int(q['raw_reads'])
assert all(float(r['padj'])<0.05 and abs(float(r['log2FoldChange']))>=1 for r in sig)
assert sum(r['significant_FDR05']=='TRUE' for r in de)==stats['FDR05_genes']
assert sum(float(r['log2FoldChange'])>0 for r in sig)==stats['higher_in_sec66del_absLFC1']
assert sum(float(r['log2FoldChange'])<0 for r in sig)==stats['lower_in_sec66del_absLFC1']

rows=[]
archives=sorted((p/'results').glob('fastqc*/*.zip'))
assert len(archives)==12
for zpath in archives:
    with zipfile.ZipFile(zpath) as z:
        name=next(n for n in z.namelist() if n.endswith('/summary.txt'))
        for line in z.read(name).decode().splitlines():
            status,module,sample=line.split('\t')
            rows.append([zpath.parent.name,sample,module,status])
with (p/'results/fastqc_module_status.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['stage','sample','module','status']);w.writerows(rows)
for name in ['PCA','MA','volcano','top20_heatmap','sample_distances']:
    assert (p/f'results/figures/{name}.png').stat().st_size>1000
assert (p/'results/multiqc/multiqc_report.html').stat().st_size>1000
(p/'results/validation.json').write_text(json.dumps({
    'status':'passed', 'samples':len(meta), 'genes_in_matrix':len(matrix),
    'genes_tested':len(de), 'FDR05_absLFC1_genes':len(sig),
    'checks':['sample order','count column sums equal assigned reads','raw read counts match ENA',
              'differential expression thresholds and summary totals','12 FastQC archives','5 figures and MultiQC report'],
    'interpretation_note':'Technical checks passed; see RESULTS.md for biological variability and FastQC flags.'
},indent=2)+'\n')
print('Results validated against source metadata, count summaries and statistical thresholds.')
