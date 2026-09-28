#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export LC_ALL=C.UTF-8
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
THREADS=3
mkdir -p logs results/{fastqc_raw,fastqc_trimmed,fastp,counts,strand_check,figures} trimmed alignments reference
exec > >(tee -a logs/pipeline.log) 2>&1
trap 'echo "FAILED at line $LINENO, $(date -Iseconds)"' ERR
echo "Pipeline started: $(date -Iseconds)"
python scripts/download.py
{
  uname -a
  fastqc --version
  fastp --version
  hisat2 --version | head -1
  samtools --version | head -2
  featureCounts -v
  Rscript -e 'cat("R ",as.character(getRversion()),"; DESeq2 ",as.character(packageVersion("DESeq2")),"\n",sep="")'
} > metadata/tool_versions.txt 2>&1
conda list --explicit > config/conda-linux-aarch64.lock.txt
conda env export --no-builds | sed '/^prefix:/d' > config/environment-resolved.yml
if [[ ! -f reference/index.complete ]]; then
  hisat2_extract_splice_sites.py reference/genes.gtf > reference/splice_sites.txt
  hisat2_extract_exons.py reference/genes.gtf > reference/exons.txt
  hisat2-build -p "$THREADS" --ss reference/splice_sites.txt --exon reference/exons.txt \
    reference/genome.fa reference/genome > logs/hisat2-build.log 2>&1
  touch reference/index.complete
fi

while IFS=$'\t' read -r sample condition replicate gsm run remainder; do
  [[ "$sample" == sample ]] && continue
  echo "Processing $sample ($run): $(date -Iseconds)"
  raw="data/raw/$run.fastq.gz"
  clean="trimmed/$sample.fastq.gz"
  bam="alignments/$sample.bam"
  if [[ ! -f "results/fastqc_raw/${run}_fastqc.html" ]]; then
    fastqc -t 1 --memory 512 -o results/fastqc_raw "$raw" > "logs/$sample.fastqc_raw.log" 2>&1
  fi
  if [[ ! -f "trimmed/$sample.complete" ]]; then
    fastp -i "$raw" -o "$clean.partial.gz" --thread "$THREADS" \
      --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
      --qualified_quality_phred 20 --unqualified_percent_limit 40 \
      --length_required 25 --compression 4 \
      --json "results/fastp/$sample.json" --html "results/fastp/$sample.html" \
      > "logs/$sample.fastp.log" 2>&1
    gzip -t "$clean.partial.gz"
    mv "$clean.partial.gz" "$clean"
    touch "trimmed/$sample.complete"
  fi
  if [[ ! -f "results/fastqc_trimmed/${sample}_fastqc.html" ]]; then
    fastqc -t 1 --memory 512 -o results/fastqc_trimmed "$clean" > "logs/$sample.fastqc_trimmed.log" 2>&1
  fi
  if [[ ! -f "$bam" ]]; then
    # Align without constraining strand; count all orientations below to verify library direction.
    hisat2 -p 2 --min-intronlen 20 --max-intronlen 10000 --dta \
      --known-splicesite-infile reference/splice_sites.txt -x reference/genome \
      -U "$clean" --summary-file "logs/$sample.hisat2.summary.txt" \
      2> "logs/$sample.hisat2.log" |
      samtools sort -@ 1 -m 384M -o "$bam.partial" -
    samtools quickcheck -v "$bam.partial"
    mv "$bam.partial" "$bam"
  fi
  samtools index "$bam"
  samtools flagstat "$bam" > "logs/$sample.flagstat.txt"
done < metadata/samples.tsv

mapfile -t samples < <(tail -n +2 metadata/samples.tsv | cut -f1)
bams=()
for sample in "${samples[@]}"; do bams+=("alignments/$sample.bam"); done
for strand in 0 1 2; do
  # Single-end reads; exclude multi-mappers (default), ambiguous overlaps (default), MAPQ <10.
  featureCounts -T "$THREADS" -s "$strand" -Q 10 -t exon -g gene_id \
    -a reference/genes.gtf -o "results/strand_check/strand${strand}.tsv" \
    "${bams[@]}" > "logs/featureCounts.strand${strand}.log" 2>&1
done
python scripts/make_counts.py
Rscript scripts/differential_expression.R > logs/deseq2.log 2>&1
multiqc results/fastqc_raw results/fastqc_trimmed results/fastp logs results/counts \
  --outdir results/multiqc --force > logs/multiqc.log 2>&1
python scripts/validate_results.py
echo "Pipeline completed: $(date -Iseconds)"
date -Iseconds > results/COMPLETED.txt
