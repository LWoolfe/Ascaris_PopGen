# 03 — BWA-MEM alignment and BAM quality control

## Purpose

BWA-MEM alignment and BAM quality control. This file archives the commands and analysis notes from the thesis/paper workflow.

## Inputs

Host-depleted paired FASTQs, reference. Source filenames and absolute HPC paths are retained in the original script blocks below.

## Workflow

See the original annotated analysis sequence below. Distinct script blocks are independent jobs/R sessions, **not one executable script**.

## Annotated scripts

## Overview

Host-depleted paired-end reads from all 156 samples were aligned to the chromosome-level *Ascaris suum* reference:

```text
/parallel_scratch/lw01165/01_Reference/Ascaris_primary_scaffolds_plus_X.fasta
```

The mapping and mapping-level quality-control workflow consisted of:

1. alignment of paired reads using BWA-MEM;
2. coordinate sorting and indexing of BAM files using samtools;
3. generation of initial alignment summaries using `samtools flagstat`;
4. identification and marking of duplicate reads using GATK/Picard MarkDuplicates;
5. calculation of sample-level median mapping quality, breadth and median depth;
6. construction of non-overlapping 25 kb genomic windows;
7. calculation of breadth and true mean depth in each 25 kb window;
8. aggregation of window-level results across samples;
9. classification and visualisation of sample-level mapping QC in R;
10. chromosome-level and genome-position summaries of depth and breadth in R.

The host-depleted reads remained as separate samples throughout mapping. Mapping and QC jobs were run in batches where necessary, but the same analytical commands and thresholds were applied across the dataset.

---

## 1. Alignment to the chromosome-level reference

### Purpose

Host-depleted paired reads were aligned to the chromosome-level reference using BWA-MEM. Read-group information was added during alignment so that each BAM retained its sample identity for downstream GATK processing.

The alignments were then:

- coordinate-sorted using samtools;
- written as one BAM per sample;
- indexed using samtools;
- summarised using `samtools flagstat`.

### Inputs

```text
Reference:
Ascaris_primary_scaffolds_plus_X.fasta

Reads:
SAMPLE_decon_R1.fastq.gz
SAMPLE_decon_R2.fastq.gz

Sample list:
Sample_list.txt
```

### Mapping script

The following recovered SLURM submission script was used for the 128 publicly available SRA datasets. The newly generated samples were mapped using the same BWA-MEM, read-group, sorting, indexing and `flagstat` commands in a earlier batchs as they became available.

```bash
#!/bin/bash 
#SBATCH --job-name=bwa_map
#SBATCH --partition=shared
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=72:00:00
#SBATCH --array=1-128%6
#SBATCH --output=/parallel_scratch/lw01165/03_Mapping/logs/bwa_%A_%a.out
#SBATCH --error=/parallel_scratch/lw01165/03_Mapping/logs/bwa_%A_%a.err

set -euo pipefail

# ------------------------------------------------------------------------
# SOFTWARE ENVIRONMENT
# ------------------------------------------------------------------------

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate BWA

command -v bwa >/dev/null || {
  echo "bwa not found"
  exit 1
}

command -v samtools >/dev/null || {
  echo "samtools not found"
  exit 1
}

# ------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------

REF="/parallel_scratch/lw01165/01_Reference/Ascaris_primary_scaffolds_plus_X.fasta"

READS_DIR="/parallel_scratch/lw01165/sra_files/samples_deacon"
SAMPLE_LIST="/parallel_scratch/lw01165/sra_files/ids.clean.txt"

OUT_DIR="/parallel_scratch/lw01165/03_Mapping"
QC_DIR="${OUT_DIR}/Align_QC"
LOGDIR="${OUT_DIR}/logs"
TMPDIR="/parallel_scratch/lw01165/tmp_bwa_${SLURM_JOB_ID:-local}"

mkdir -p "$OUT_DIR" "$QC_DIR" "$LOGDIR" "$TMPDIR"

# ------------------------------------------------------------------------
# REFERENCE-INDEX CHECKS
# ------------------------------------------------------------------------

for ext in amb ann bwt pac sa; do
  [[ -s "${REF}.${ext}" ]] || {
    echo "❌ Missing BWA index: ${REF}.${ext}"
    exit 1
  }
done

[[ -s "${REF}.fai" ]] || {
  echo "❌ Missing samtools index: ${REF}.fai"
  exit 1
}

# ------------------------------------------------------------------------
# SELECT SAMPLE
# ------------------------------------------------------------------------

[[ -s "$SAMPLE_LIST" ]] || {
  echo "❌ Missing sample list: $SAMPLE_LIST"
  exit 1
}

SAMPLE="$(
  sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST" |
  tr -d '\r' |
  sed 's/[[:space:]]*$//'
)"

[[ -n "${SAMPLE:-}" ]] || {
  echo "No sample for task ${SLURM_ARRAY_TASK_ID}"
  exit 0
}

# ------------------------------------------------------------------------
# INPUTS AND OUTPUT
# ------------------------------------------------------------------------

R1="${READS_DIR}/${SAMPLE}_decon_R1.fastq.gz"
R2="${READS_DIR}/${SAMPLE}_decon_R2.fastq.gz"
BAM="${OUT_DIR}/${SAMPLE}.bam"

[[ -s "$R1" && -s "$R2" ]] || {
  echo "❌ Missing reads for $SAMPLE: $R1 or $R2"
  exit 1
}

# ------------------------------------------------------------------------
# SKIP COMPLETED SAMPLES
# ------------------------------------------------------------------------

if [[ -s "$BAM" && -s "${BAM}.bai" ]]; then
  echo "✅ $SAMPLE BAM and BAI already exist — skipping"
  exit 0
fi

THREADS="${SLURM_CPUS_PER_TASK:-8}"
SORT_MEM="2G"

# ------------------------------------------------------------------------
# READ GROUP AND PLATFORM UNIT
# ------------------------------------------------------------------------

HEADER="$(zcat "$R1" | head -n1 || true)"

FLOWCELL="$(
  awk -F':' '{if (NF >= 3) print $3}' <<< "$HEADER"
)"

LANE="$(
  awk -F':' '{if (NF >= 4) print $4}' <<< "$HEADER"
)"

PU="${FLOWCELL:-${SAMPLE}}.${LANE:-1}"
RGID="${PU}"

echo "==> $(date) Mapping $SAMPLE on $(hostname) (threads=$THREADS)"

# ------------------------------------------------------------------------
# ALIGN AND COORDINATE-SORT
# ------------------------------------------------------------------------

bwa mem \
  -t "$THREADS" \
  -R "@RG\tID:${RGID}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1\tPU:${PU}" \
  "$REF" \
  "$R1" \
  "$R2" |
samtools sort \
  -@ "$THREADS" \
  -m "$SORT_MEM" \
  -T "${TMPDIR}/${SAMPLE}" \
  -o "$BAM" \
  -

# ------------------------------------------------------------------------
# INDEX AND GENERATE INITIAL ALIGNMENT SUMMARY
# ------------------------------------------------------------------------

samtools index \
  -@ "$THREADS" \
  "$BAM"

samtools flagstat \
  "$BAM" \
  > "${QC_DIR}/${SAMPLE}.flagstat.txt"

echo "==> $(date) Done $SAMPLE"
```

### Main outputs

```text
/parallel_scratch/lw01165/03_Mapping/SAMPLE.bam
/parallel_scratch/lw01165/03_Mapping/SAMPLE.bam.bai
/parallel_scratch/lw01165/03_Mapping/Align_QC/SAMPLE.flagstat.txt
```

### Interpretation

The BAM files contain coordinate-sorted alignments. The read-group sample field (`SM`) allows sample identity to be carried into downstream GATK analyses. `samtools flagstat` records initial alignment metrics, including mapped and properly paired reads.

---

## 2. Duplicate marking

### Purpose

Duplicate reads were identified and marked using GATK MarkDuplicates. Duplicate reads were flagged in the BAM rather than removed. A metrics file was generated for each sample, and the duplicate-marked BAM was indexed using samtools.

### Duplicate-marking script

```bash
#!/bin/bash 
#SBATCH --job-name=markdup
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=30:00:00
#SBATCH --array=1-156%6
#SBATCH --output=/parallel_scratch/lw01165/03_Mapping/3a_markdup/logs/markdup_%A_%a.out
#SBATCH --error=/parallel_scratch/lw01165/03_Mapping/3a_markdup/logs/markdup_%A_%a.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env

command -v gatk >/dev/null || { echo "❌ gatk not found in tools_env"; exit 1; }
command -v samtools >/dev/null || { echo "❌ samtools not found in tools_env"; exit 1; }

MAP_DIR="/parallel_scratch/lw01165/03_Mapping"
MARKDUP_DIR="${MAP_DIR}/3a_markdup"
LOGDIR="${MARKDUP_DIR}/logs"
SAMPLE_LIST="${MAP_DIR}/Sample_list.txt"

mkdir -p "$MARKDUP_DIR" "$LOGDIR"

[[ -s "$SAMPLE_LIST" ]] || { echo "❌ Missing sample list: $SAMPLE_LIST"; exit 1; }

SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST" | tr -d '\r' | sed 's/[[:space:]]*$//')"
[[ -n "${SAMPLE:-}" ]] || { echo "No sample for task ${SLURM_ARRAY_TASK_ID}"; exit 0; }

exec >"${LOGDIR}/${SAMPLE}.markdup.out" 2>"${LOGDIR}/${SAMPLE}.markdup.err"

IN_BAM="${MAP_DIR}/${SAMPLE}.bam"
OUT_BAM="${MARKDUP_DIR}/${SAMPLE}.markdup.bam"
METRICS="${MARKDUP_DIR}/${SAMPLE}.metrics.txt"

[[ -s "$IN_BAM" ]] || { echo "❌ Input BAM not found for $SAMPLE: $IN_BAM"; exit 1; }

if [[ -s "$OUT_BAM" && -s "${OUT_BAM}.bai" ]]; then
  echo "✅ $SAMPLE markdup BAM+BAI already exist — skipping"
  exit 0
fi

THREADS="${SLURM_CPUS_PER_TASK:-4}"

echo "==> $(date) MarkDuplicates for $SAMPLE on $(hostname) (threads=$THREADS)"

gatk MarkDuplicates \
  --INPUT "$IN_BAM" \
  --OUTPUT "$OUT_BAM" \
  --METRICS_FILE "$METRICS" \
  --VALIDATION_STRINGENCY LENIENT

samtools index -@ "$THREADS" "$OUT_BAM"

echo "==> $(date) Done $SAMPLE"
```

### Main outputs

```text
/parallel_scratch/lw01165/03_Mapping/3a_markdup/SAMPLE.markdup.bam
/parallel_scratch/lw01165/03_Mapping/3a_markdup/SAMPLE.markdup.bam.bai
/parallel_scratch/lw01165/03_Mapping/3a_markdup/SAMPLE.metrics.txt
```

### Interpretation

`MarkDuplicates` identifies read pairs with matching alignment positions and flags probable PCR or optical duplicates. The reads remain present, but downstream tools can recognise their duplicate status. The metrics files retain duplicate-count information for each sample.

---

## 3. Sample-level mapping statistics

### Purpose

This script calculates three mapping-level statistics for each duplicate-marked BAM:

- median mapping quality;
- breadth of coverage, defined as the fraction of reference positions with depth greater than zero;
- median depth across reference positions.

It also calculates the median of each metric across all samples.

### Mapping-QC script

```bash
#!/bin/bash
#SBATCH --job-name=mapping_qc_stats
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/parallel_scratch/lw01165/03_Mapping/3b_QC/mapping_qc_stats_%j.out
#SBATCH --error=/parallel_scratch/lw01165/03_Mapping/3b_QC/mapping_qc_stats_%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env

command -v samtools >/dev/null || { echo "❌ samtools not found"; exit 1; }
command -v datamash >/dev/null || { echo "❌ datamash not found"; exit 1; }

MAP_DIR="/parallel_scratch/lw01165/03_Mapping"
MARKDUP_DIR="${MAP_DIR}/3a_markdup"
QC_DIR="${MAP_DIR}/3b_QC"
SAMPLE_LIST="${MAP_DIR}/Sample_list.txt"

mkdir -p "$QC_DIR"

[[ -s "$SAMPLE_LIST" ]] || { echo "❌ Sample list not found: $SAMPLE_LIST"; exit 1; }

THREADS="${SLURM_CPUS_PER_TASK:-4}"
SUMMARY="${QC_DIR}/sample_mapping_stats.txt"

echo -e "sample\tmedian_MAPQ\tbreadth\tmedian_depth" > "$SUMMARY"

echo "==> Computing mapping stats for samples in $SAMPLE_LIST"
echo "    BAMs from: $MARKDUP_DIR"
echo "    Output to: $QC_DIR"
echo

while read -r SAMPLE; do
  SAMPLE="$(echo "$SAMPLE" | tr -d '\r' | sed 's/[[:space:]]*$//')"
  [[ -z "$SAMPLE" ]] && continue

  BAM="${MARKDUP_DIR}/${SAMPLE}.markdup.bam"
  [[ -s "$BAM" ]] || { echo "❌ Missing BAM for $SAMPLE: $BAM" >&2; exit 1; }

  echo "==> Sample: $SAMPLE"

  STATS="${QC_DIR}/${SAMPLE}.stats"
  samtools stats -@ "$THREADS" "$BAM" > "$STATS"

  MED_MAPQ=$(grep ^MAPQ "$STATS" | awk '
  BEGIN { total=0 }
  { total += $3; count[NR]=$3; q[NR]=$2 }
  END {
    if (total==0) { print "NA"; exit }
    half = total/2;
    running = 0;
    for (i=1; i<=NR; i++) {
      running += count[i];
      if (running >= half) {
        print q[i];
        exit
      }
    }
  }')

  BREADTH=$(samtools depth -a -@ "$THREADS" "$BAM" | \
    awk '{ if($3 > 0) covered++; total++ } END { if(total>0) print covered/total; else print 0 }')

  MED_DEPTH=$(samtools depth -a -@ "$THREADS" "$BAM" | \
    awk '{print $3}' | datamash median 1)

  echo -e "${SAMPLE}\t${MED_MAPQ}\t${BREADTH}\t${MED_DEPTH}" >> "$SUMMARY"

done < "$SAMPLE_LIST"

echo "==> Per-sample stats written to: $SUMMARY"

OVERALL="${QC_DIR}/overall_mapping_stats.txt"

MED_MAPQ_ALL=$(awk 'NR>1 {print $2}' "$SUMMARY" | datamash median 1)
MED_BREADTH_ALL=$(awk 'NR>1 {print $3}' "$SUMMARY" | datamash median 1)
MED_DEPTH_ALL=$(awk 'NR>1 {print $4}' "$SUMMARY" | datamash median 1)

{
  echo -e "metric\tmedian_MAPQ\tmedian_breadth\tmedian_depth"
  echo -e "across_samples\t${MED_MAPQ_ALL}\t${MED_BREADTH_ALL}\t${MED_DEPTH_ALL}"
} > "$OVERALL"

echo "==> Overall medians written to: $OVERALL"
echo "==> Done."
```

### Main outputs

```text
/parallel_scratch/lw01165/03_Mapping/3b_QC/sample_mapping_stats.txt
/parallel_scratch/lw01165/03_Mapping/3b_QC/overall_mapping_stats.txt
/parallel_scratch/lw01165/03_Mapping/3b_QC/SAMPLE.stats
```

---

## 4. Construction of non-overlapping 25 kb windows

### Purpose

The reference FASTA index was converted into a two-column chromosome-size file. `bedtools makewindows` was then used to divide all 24 chromosome-level sequences into consecutive, non-overlapping 25 kb windows.

NB:The final window on each chromosome may be shorter than 25 kb where chromosome length is not exactly divisible by 25,000.

### Window-generation script

```bash
#!/bin/bash
#SBATCH --job-name=make_25kb_windows
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/parallel_scratch/lw01165/01_Reference/make_25kb_%j.out
#SBATCH --error=/parallel_scratch/lw01165/01_Reference/make_25kb_%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env

command -v bedtools >/dev/null || { echo "❌ bedtools not found in tools_env"; exit 1; }
command -v samtools >/dev/null || { echo "❌ samtools not found in tools_env"; exit 1; }

REF_DIR="/parallel_scratch/lw01165/01_Reference"
REF="${REF_DIR}/Ascaris_primary_scaffolds_plus_X.fasta"
FAI="${REF}.fai"
GENOME="${REF_DIR}/Ascaris_primary_scaffolds_plus_X.genome.txt"
WINDOWS="${REF_DIR}/Ascaris_primary_scaffolds_plus_X.25kb.bed"

echo "==> Using reference: $REF"
echo

if [[ ! -s "$FAI" ]]; then
    echo "==> No .fai file found, running samtools faidx..."
    samtools faidx "$REF"
else
    echo "==> .fai index already exists."
fi

echo "==> Creating genome file: $GENOME"
cut -f1,2 "$FAI" > "$GENOME"

echo "==> Making 25 kb windows..."
bedtools makewindows \
    -g "$GENOME" \
    -w 25000 \
    > "$WINDOWS"

echo
echo "==> Done!"
echo "Generated:"
echo " - $GENOME"
echo " - $WINDOWS"
```

### Main outputs

```text
Ascaris_primary_scaffolds_plus_X.genome.txt
Ascaris_primary_scaffolds_plus_X.25kb.bed
```

---

## 5. Calculation of breadth and mean depth in 25 kb windows

### Purpose

This script calculates window-level coverage statistics separately for each duplicate-marked BAM.

Two `bedtools coverage` commands are used:

1. standard `bedtools coverage`, which provides the number of overlapping alignments, number of bases covered, interval length and fraction of bases covered;
2. `bedtools coverage -mean`, which provides the true mean read depth across each interval.

> **Important:** The fourth output column from standard `bedtools coverage` is the number of overlapping reads or alignments. It is not read depth. True mean depth is calculated separately using `bedtools coverage -mean`.

### Final per-sample 25 kb coverage script

```bash
#!/bin/bash
#SBATCH --job-name=cov25kb
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=5-00:00:00
#SBATCH --array=1-156%15
#SBATCH --output=/parallel_scratch/lw01165/03_Mapping/3c_25kb_QC/logs/cov25kb_%A_%a.out
#SBATCH --error=/parallel_scratch/lw01165/03_Mapping/3c_25kb_QC/logs/cov25kb_%A_%a.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env

command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found"; exit 1; }
command -v samtools  >/dev/null 2>&1 || { echo "ERROR: samtools not found"; exit 1; }

REF_BED="/parallel_scratch/lw01165/01_Reference/Ascaris_primary_scaffolds_plus_X.25kb.bed"
MAP_DIR="/parallel_scratch/lw01165/03_Mapping"
MARKDUP_DIR="${MAP_DIR}/3a_markdup"
QC_DIR="${MAP_DIR}/3c_25kb_QC"
LOGDIR="${QC_DIR}/logs"
SAMPLE_LIST="${MAP_DIR}/Sample_list.txt"
TMP_DIR="${QC_DIR}/tmp_sample_25kb_cov_stats"

mkdir -p "$QC_DIR" "$LOGDIR" "$TMP_DIR"

[[ -s "$REF_BED" ]]     || { echo "ERROR: 25kb BED not found: $REF_BED"; exit 1; }
[[ -s "$SAMPLE_LIST" ]] || { echo "ERROR: Sample list not found: $SAMPLE_LIST"; exit 1; }

SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST" | tr -d $'\r' | sed 's/[[:space:]]*$//')"
[[ -n "${SAMPLE:-}" ]] || { echo "No sample for task ${SLURM_ARRAY_TASK_ID}"; exit 0; }

IN_BAM="${MARKDUP_DIR}/${SAMPLE}.markdup.bam"
[[ -s "$IN_BAM" ]] || { echo "ERROR: Missing BAM for $SAMPLE: $IN_BAM"; exit 1; }

OUT_BREADTH="${TMP_DIR}/${SAMPLE}.25kb.breadth.tmp.tsv"
OUT_MEANDEPTH="${TMP_DIR}/${SAMPLE}.25kb.meandepth.tmp.tsv"
OUT_SAMPLE="${TMP_DIR}/${SAMPLE}.25kb.cov_stats.tsv"

if [[ -s "$OUT_SAMPLE" ]]; then
  echo "Output exists for ${SAMPLE}, skipping: $OUT_SAMPLE"
  exit 0
fi

echo "==> $(date) Coverage stats for $SAMPLE"
echo "    BAM: $IN_BAM"
echo "    BED: $REF_BED"
echo "    OUT: $OUT_SAMPLE"

SORTED_BED="${QC_DIR}/Ascaris_primary_scaffolds_plus_X.25kb.sorted.bed"

if [[ ! -s "$SORTED_BED" ]]; then
  echo "Creating sorted BED: $SORTED_BED"
  tmp="${SORTED_BED}.${SLURM_ARRAY_TASK_ID}.tmp"
  sort -k1,1 -k2,2n "$REF_BED" > "$tmp"
  mv "$tmp" "$SORTED_BED"
fi

if [[ ! -s "${IN_BAM}.bai" ]]; then
  echo "Indexing BAM..."
  samtools index "$IN_BAM"
fi

bedtools coverage \
  -sorted \
  -a "$SORTED_BED" \
  -b "$IN_BAM" \
  > "$OUT_BREADTH"

bedtools coverage \
  -sorted \
  -mean \
  -a "$SORTED_BED" \
  -b "$IN_BAM" \
  > "$OUT_MEANDEPTH"

paste "$OUT_BREADTH" <(awk 'BEGIN{OFS="\t"} {print $NF}' "$OUT_MEANDEPTH") | \
awk -v S="$SAMPLE" 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,S}' \
> "$OUT_SAMPLE"

rm -f "$OUT_BREADTH" "$OUT_MEANDEPTH"

echo "==> $(date) Done $SAMPLE"
```

### Per-sample output columns

```text
chrom
start
end
n_overlapping_reads
bases_covered
window_length
fraction_covered
mean_depth
sample
```

---

## 6. Combination of per-sample 25 kb coverage tables

### Purpose

This script combines all per-sample window tables into one master file in `Sample_list.txt` order. It checks that every expected sample file is present before the output is used.

### Combination script

```bash
#!/bin/bash
#SBATCH --job-name=combine25kb
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --output=/parallel_scratch/lw01165/03_Mapping/3c_25kb_QC/logs/combine25kb_%j.out
#SBATCH --error=/parallel_scratch/lw01165/03_Mapping/3c_25kb_QC/logs/combine25kb_%j.err

set -eo pipefail

MAP_DIR="/parallel_scratch/lw01165/03_Mapping"
QC_DIR="${MAP_DIR}/3c_25kb_QC"
TMP_DIR="${QC_DIR}/tmp_sample_25kb_cov_stats"
SAMPLE_LIST="${MAP_DIR}/Sample_list.txt"
MASTER_OUT="${QC_DIR}/New_25kb_cov_stats.tsv"
CLEAN_TMP=0

mkdir -p "$QC_DIR"

[[ -d "$TMP_DIR" ]] || { echo "ERROR: Temp directory not found: $TMP_DIR"; exit 1; }
[[ -s "$SAMPLE_LIST" ]] || { echo "ERROR: Sample list not found: $SAMPLE_LIST"; exit 1; }

echo "==> $(date) Combining 25 kb coverage stats"
echo "    TMP_DIR: $TMP_DIR"
echo "    OUT:     $MASTER_OUT"

echo -e "chrom\tstart\tend\tn_overlapping_reads\tbases_covered\twindow_length\tfraction_covered\tmean_depth\tsample" > "$MASTER_OUT"

missing=0

while read -r SAMPLE; do
  SAMPLE="$(echo "$SAMPLE" | tr -d $'\r' | sed 's/[[:space:]]*$//')"
  [[ -z "$SAMPLE" ]] && continue

  FILE="${TMP_DIR}/${SAMPLE}.25kb.cov_stats.tsv"

  if [[ ! -s "$FILE" ]]; then
    echo "WARNING: Missing per-sample file: $FILE"
    missing=$((missing + 1))
    continue
  fi

  cat "$FILE" >> "$MASTER_OUT"
done < "$SAMPLE_LIST"

echo "==> Combined file written:"
echo "    $MASTER_OUT"

echo "==> Number of data rows:"
tail -n +2 "$MASTER_OUT" | wc -l

echo "==> Missing sample files:"
echo "    $missing"

if [[ "$missing" -gt 0 ]]; then
  echo "ERROR: Some sample files were missing. Check array logs before using master file."
  exit 1
fi

if [[ "$CLEAN_TMP" -eq 1 ]]; then
  echo "Cleaning temporary per-sample files..."
  rm -f "${TMP_DIR}"/*.25kb.cov_stats.tsv
fi

echo "==> $(date) Done combining"
```

### Master output

```text
/parallel_scratch/lw01165/03_Mapping/3c_25kb_QC/New_25kb_cov_stats.tsv
```

---

## 7. Optional shell-based chromosome breadth summary

### Purpose

This provides an additional shell-based summary of breadth across 25 kb observations for each chromosome. The final R analysis also calculates chromosome-level summaries and is the main source of the figures.

The input filename below has been aligned with the final combined coverage output.

```bash
#!/bin/bash
#SBATCH --job-name=chr25kb_summary
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=chr25kb_summary_%j.out
#SBATCH --error=chr25kb_summary_%j.err

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env

set -euo pipefail

QC_DIR="/parallel_scratch/lw01165/03_Mapping/3c_25kb_QC"
IN="${QC_DIR}/New_25kb_cov_stats.tsv"
OUT="${QC_DIR}/per_chromosome_25kb_breadth_summary.tsv"
CHR_ORDER="${QC_DIR}/chrom_order.txt"
CHR_QUANT="${QC_DIR}/.chr_quant.tsv"
CHR_UNSORTED="${QC_DIR}/.chr_summary_unsorted.tsv"

[[ -s "$IN" ]] || { echo "ERROR: missing input $IN" >&2; exit 1; }

awk -F'\t' 'NR>1 {print $1}' "$IN" | awk '!seen[$1]++' > "$CHR_ORDER"

awk -F'\t' 'NR>1 && $1!="" && $7!="" {print $1 "\t" $7}' "$IN" \
| sort -t $'\t' -k1,1 -k2,2n \
| datamash -t $'\t' -g 1 count 2 q1 2 median 2 q3 2 perc:10 2 \
> "$CHR_QUANT"

{
  echo -e "chrom\tn_windows\tq25_breadth\tmedian_breadth\tq75_breadth\tiqr_breadth\tp10_breadth"
  awk -F'\t' 'BEGIN{OFS="\t"}
  {
    chrom=$1
    n=$2
    q1=$3
    med=$4
    q3=$5
    p10=$6
    iqr=q3-q1
    print chrom,n,q1,med,q3,iqr,p10
  }' "$CHR_QUANT"
} > "$CHR_UNSORTED"

{
  echo -e "chrom\tn_windows\tq25_breadth\tmedian_breadth\tq75_breadth\tiqr_breadth\tp10_breadth"
  awk -F'\t' 'NR==FNR {ord[$1]=++i; next}
              FNR>1 {print ord[$1], $0}' "$CHR_ORDER" "$CHR_UNSORTED" \
  | sort -t $'\t' -k1,1n \
  | cut -f2-
} > "$OUT"

rm -f "$CHR_QUANT" "$CHR_UNSORTED"

echo "DONE"
echo "Output written to: $OUT"
```

---

## 8. Sample-level mapping-QC classification and figure

### Purpose

The sample-level mapping-QC table was imported into R and used to classify samples based on depth and breadth.

Final thresholds were:

```text
Median depth:         5×
Breadth of coverage:  0.85
```

A sample was classified as a formal mapping-QC failure only when it failed both thresholds:

```text
median depth < 5×
AND
breadth < 0.85
```

The three lowest-scoring retained samples were displayed as `Watch` samples for review. The Watch category did not itself determine exclusion.

### R script

```r
# ============================================================
# Mapping QC: classify samples, write tables, publication plot
# ============================================================

setwd("S:/Research/Lauren_PopGen")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(ggrepel)
})

qc_file <- "sample_mapping_stats.txt"

if (!file.exists(qc_file)) {
  stop("Cannot find input file: ", qc_file)
}

qc <- fread(qc_file)

required_cols <- c(
  "sample",
  "median_MAPQ",
  "breadth",
  "median_depth"
)

missing_cols <- setdiff(required_cols, names(qc))

if (length(missing_cols) > 0) {
  stop(
    "The mapping-QC file is missing these columns: ",
    paste(missing_cols, collapse = ", "),
    "\nAvailable columns are: ",
    paste(names(qc), collapse = ", ")
  )
}

qc[, `:=`(
  sample       = as.character(sample),
  median_depth = as.numeric(median_depth),
  breadth      = as.numeric(breadth),
  median_MAPQ  = as.numeric(median_MAPQ)
)]

if (anyNA(qc$sample)) {
  stop("Missing sample identifiers were found.")
}

if (anyNA(qc$median_depth) || anyNA(qc$breadth)) {
  stop("Missing or non-numeric depth/breadth values were found.")
}

depth_cutoff <- 5
breadth_cutoff <- 0.85

qc[, fail_depth := median_depth < depth_cutoff]
qc[, fail_breadth := breadth < breadth_cutoff]

qc[, qc_class := fifelse(
  fail_depth & fail_breadth,
  "Failure",
  "Keep"
)]

qc[, qc_class := factor(
  qc_class,
  levels = c("Keep", "Failure")
)]

qc[, depth_log10 := log10(pmax(median_depth, 1e-6))]
qc[, score := breadth + depth_log10]

n_watch_seed <- 3

watch_seed <- qc[qc_class == "Keep"][order(score)][
  1:min(.N, n_watch_seed), sample
]

qc[, qc_class_plot := "Keep"]
qc[sample %in% watch_seed, qc_class_plot := "Watch"]
qc[qc_class == "Failure", qc_class_plot := "Failure"]

qc[
  breadth == breadth_cutoff | median_depth == depth_cutoff,
  qc_class_plot := "Watch"
]

qc[, qc_class_plot := factor(
  qc_class_plot,
  levels = c("Keep", "Watch", "Failure")
)]

fwrite(qc, "Table_S1_mapping_QC.csv")

fwrite(
  qc[qc_class_plot == "Failure",
     .(sample, breadth, median_depth, median_MAPQ,
       fail_breadth, fail_depth, score)],
  "failure_samples.csv"
)

fwrite(
  qc[qc_class_plot == "Watch",
     .(sample, breadth, median_depth, median_MAPQ,
       fail_breadth, fail_depth, score)],
  "watch_samples.csv"
)

qc_kept <- qc[qc_class_plot != "Failure"]
fwrite(qc_kept, "mapping_QC_excluding_failures.csv")

legend_labels_dt <- qc[, .(
  samples = paste(sample, collapse = ", "),
  n = .N
), by = qc_class_plot]

legend_labels <- setNames(
  paste0(
    legend_labels_dt$qc_class_plot,
    " (n=", legend_labels_dt$n, ")",
    ifelse(
      legend_labels_dt$qc_class_plot == "Keep",
      "",
      paste0(": ", legend_labels_dt$samples)
    )
  ),
  as.character(legend_labels_dt$qc_class_plot)
)

qc_cols <- c(
  "Keep" = "#1B9E77",
  "Watch" = "#7570B3",
  "Failure" = "darkred"
)

p1 <- ggplot(qc, aes(x = breadth, y = median_depth, colour = qc_class_plot)) +
  geom_point(data = qc[qc_class_plot == "Keep"], alpha = 0.85, size = 2) +
  geom_point(data = qc[qc_class_plot == "Watch"], alpha = 1, size = 2.6) +
  geom_point(data = qc[qc_class_plot == "Failure"], alpha = 0.9, size = 2.6) +
  geom_vline(
    xintercept = breadth_cutoff,
    linetype = "dashed", linewidth = 0.3, colour = "grey70"
  ) +
  geom_hline(
    yintercept = depth_cutoff,
    linetype = "dashed", linewidth = 0.3, colour = "grey70"
  ) +
  scale_colour_manual(
    values = qc_cols,
    labels = legend_labels,
    drop = FALSE
  ) +
  scale_y_log10(
    limits = c(0.6, NA),
    breaks = c(1, 2, 5, 10, 20, 50, 100, 200),
    labels = scales::label_number()
  ) +
  scale_x_continuous(
    breaks = seq(0.5, 1, by = 0.1),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  coord_cartesian(xlim = c(0.45, 1)) +
  labs(
    x = "Breadth of coverage",
    y = "Median depth (×, log10 scale)",
    colour = NULL,
    title = "Sample-level sequencing depth and breadth"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    plot.margin = margin(t = 5.5, r = 45, b = 5.5, l = 5.5)
  )

print(p1)

ggsave("Fig_QC_depth_breadth.pdf", p1, width = 7, height = 4.5, useDingbats = FALSE)
ggsave("Fig_QC_depth_breadth.png", p1, width = 7, height = 4.5, dpi = 600)

ggsave(
  "Fig_QC_depth_breadth.tiff",
  plot = p1,
  device = "tiff",
  width = 11,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  type = "cairo"
)

summary_stats <- qc[, .(
  n_samples = .N,
  n_keep = sum(qc_class_plot == "Keep", na.rm = TRUE),
  n_watch = sum(qc_class_plot == "Watch", na.rm = TRUE),
  n_failure = sum(qc_class_plot == "Failure", na.rm = TRUE),
  median_MAPQ = median(median_MAPQ, na.rm = TRUE),
  median_breadth = median(breadth, na.rm = TRUE),
  median_depth = median(median_depth, na.rm = TRUE),
  depth_range = paste0(
    min(median_depth, na.rm = TRUE),
    "–",
    max(median_depth, na.rm = TRUE)
  )
)]

fwrite(summary_stats, "mapping_QC_summary.csv")
print(summary_stats)
```

### Main outputs

```text
Table_S1_mapping_QC.csv
failure_samples.csv
watch_samples.csv
mapping_QC_excluding_failures.csv
mapping_QC_summary.csv
Fig_QC_depth_breadth.pdf
Fig_QC_depth_breadth.png
Fig_QC_depth_breadth.tiff
```

---

## 9. Per-chromosome summaries of 25 kb breadth and depth

### Purpose

The final combined 25 kb coverage table was imported into R. Breadth and true mean depth were summarised across all sample-window observations for each chromosome.

For each chromosome, the script calculated the number of observations, 25th percentile, median, 75th percentile and interquartile range. Chromosomes were displayed as 1–19 followed by X1–X5.

###  R script

```r
# ============================================================
# Genome-wide coverage heterogeneity in 25 kb windows
# ============================================================

setwd("S:/Research/Lauren_PopGen/Mapping_Assessment")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

win_file <- "New_25kb_cov_stats.tsv"

BREADTH_LIM <- c(0.90, 1.00)
BREADTH_BREAKS <- seq(BREADTH_LIM[1], BREADTH_LIM[2], by = 0.02)
DEPTH_LIM <- NULL
IQR_BAR_WIDTH <- 0.15

chrom_map <- c(
  "CM024168.1"="1",  "CM024169.1"="2",  "CM024170.1"="3",
  "CM024171.1"="4",  "CM024172.1"="5",  "CM024173.1"="6",
  "CM024174.1"="7",  "CM024175.1"="8",  "CM024176.1"="9",
  "CM024177.1"="10", "CM024178.1"="11", "CM024179.1"="12",
  "CM024180.1"="13", "CM024181.1"="14", "CM024182.1"="15",
  "CM024183.1"="16", "CM024184.1"="17", "CM024185.1"="18",
  "CM024186.1"="19",
  "CM024187.1"="X1", "CM024188.1"="X2", "CM024189.1"="X3",
  "CM024190.1"="X4", "CM024191.1"="X5"
)

chrom_levels_plot <- c(as.character(1:19), "X1", "X2", "X3", "X4", "X5")

if (!file.exists(win_file)) {
  stop("Cannot find coverage file: ", win_file)
}

win <- fread(win_file)

required_cols <- c(
  "chrom", "start", "end", "n_overlapping_reads",
  "bases_covered", "window_length", "fraction_covered",
  "mean_depth", "sample"
)

missing_cols <- setdiff(required_cols, names(win))

if (length(missing_cols) > 0) {
  stop(
    "The coverage file is missing these columns: ",
    paste(missing_cols, collapse = ", "),
    "\nAvailable columns are: ",
    paste(names(win), collapse = ", ")
  )
}

setnames(win, "fraction_covered", "breadth")

win[, `:=`(
  chrom = as.character(chrom),
  sample = as.character(sample),
  start = as.integer(start),
  end = as.integer(end),
  n_overlapping_reads = as.numeric(n_overlapping_reads),
  bases_covered = as.numeric(bases_covered),
  window_length = as.numeric(window_length),
  breadth = as.numeric(breadth),
  mean_depth = as.numeric(mean_depth)
)]

win[, chrom_label := chrom_map[chrom]]
win[is.na(chrom_label), chrom_label := chrom]

extra_levels <- setdiff(unique(win$chrom_label), chrom_levels_plot)
final_levels <- c(chrom_levels_plot, extra_levels)
win[, chrom_label := factor(chrom_label, levels = final_levels)]

chr_breadth <- win[, .(
  n_windows = .N,
  q25_breadth = as.numeric(quantile(breadth, 0.25, na.rm = TRUE)),
  median_breadth = median(breadth, na.rm = TRUE),
  q75_breadth = as.numeric(quantile(breadth, 0.75, na.rm = TRUE)),
  iqr_breadth = IQR(breadth, na.rm = TRUE)
), by = chrom_label]

chr_depth <- win[, .(
  n_windows = .N,
  q25_depth = as.numeric(quantile(mean_depth, 0.25, na.rm = TRUE)),
  median_depth = median(mean_depth, na.rm = TRUE),
  q75_depth = as.numeric(quantile(mean_depth, 0.75, na.rm = TRUE)),
  iqr_depth = IQR(mean_depth, na.rm = TRUE)
), by = chrom_label]

chr_breadth[, chrom_label := factor(as.character(chrom_label), levels = final_levels)]
chr_depth[, chrom_label := factor(as.character(chrom_label), levels = final_levels)]
setorder(chr_breadth, chrom_label)
setorder(chr_depth, chrom_label)

fwrite(chr_breadth, "per_chromosome_25kb_breadth_summary.tsv", sep = "\t")
fwrite(chr_depth, "per_chromosome_25kb_depth_summary.tsv", sep = "\t")

p2_breadth <- ggplot(chr_breadth, aes(x = chrom_label, y = median_breadth)) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = q25_breadth, ymax = q75_breadth),
    width = IQR_BAR_WIDTH
  ) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(chr_breadth$chrom_label))) +
  scale_y_continuous(
    limits = BREADTH_LIM,
    breaks = BREADTH_BREAKS,
    oob = squish
  ) +
  labs(
    x = NULL,
    y = "Breadth",
    title = "A. Per-chromosome breadth (25 kb windows)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

p2_depth <- ggplot(chr_depth, aes(x = chrom_label, y = median_depth)) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = q25_depth, ymax = q75_depth),
    width = IQR_BAR_WIDTH
  ) +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(chr_depth$chrom_label))) +
  {
    if (is.null(DEPTH_LIM)) scale_y_continuous()
    else scale_y_continuous(limits = DEPTH_LIM)
  } +
  labs(
    x = "Chromosome",
    y = "Mean depth per 25 kb window",
    title = "B. Per-chromosome depth (25 kb windows)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("Figure2A_breadth_median_IQR.pdf", p2_breadth, width = 7, height = 8)
ggsave("Figure2A_breadth_median_IQR.png", p2_breadth, width = 7, height = 8, dpi = 300)
ggsave("Figure2B_depth_median_IQR.pdf", p2_depth, width = 7, height = 8)
ggsave("Figure2B_depth_median_IQR.png", p2_depth, width = 7, height = 8, dpi = 300)

p2_main <- p2_breadth / p2_depth +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Genome-wide coverage heterogeneity (25 kb windows)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0)
    )
  )

ggsave("Figure2_main_breadth_and_depth.pdf", p2_main, width = 8.5, height = 11)
ggsave("Figure2_main_breadth_and_depth.png", p2_main, width = 8.5, height = 11, dpi = 300)

print(p2_main)
```

---

## 10. Coverage across genomic position

### Purpose

The supplementary R analysis summarised each autosomal 25 kb interval across samples and plotted median breadth and median depth across cumulative genomic position.

For each genomic window, it calculated sample number, median and mean breadth, breadth interquartile range, median and mean depth, and depth interquartile range.

It also retained tables identifying:

- windows with median breadth below 0.95;
- windows in the lowest 1% of median depth;
- high-depth windows above the plotted depth-axis limit.

### R script

```r
# ============================================================
# Supplementary figure: coverage across genomic position
# ============================================================

setwd("S:/Research/Lauren_PopGen/Mapping_Assessment")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

base_dir <- "S:/Research/Lauren_PopGen/Mapping_Assessment"
win_file <- file.path(base_dir, "New_25kb_cov_stats.tsv")
outdir <- file.path(base_dir, "coverage_QC_outputs")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(win_file)) {
  stop("Cannot find file: ", win_file)
}

chrom_map <- c(
  "CM024168.1"="1",  "CM024169.1"="2",  "CM024170.1"="3",
  "CM024171.1"="4",  "CM024172.1"="5",  "CM024173.1"="6",
  "CM024174.1"="7",  "CM024175.1"="8",  "CM024176.1"="9",
  "CM024177.1"="10", "CM024178.1"="11", "CM024179.1"="12",
  "CM024180.1"="13", "CM024181.1"="14", "CM024182.1"="15",
  "CM024183.1"="16", "CM024184.1"="17", "CM024185.1"="18",
  "CM024186.1"="19",
  "CM024187.1"="X1", "CM024188.1"="X2", "CM024189.1"="X3",
  "CM024190.1"="X4", "CM024191.1"="X5"
)

autosome_levels <- as.character(1:19)

win <- fread(win_file)

required_cols <- c(
  "chrom", "start", "end", "n_overlapping_reads",
  "bases_covered", "window_length", "fraction_covered",
  "mean_depth", "sample"
)

missing_cols <- setdiff(required_cols, names(win))

if (length(missing_cols) > 0) {
  stop(
    "The coverage file is missing these columns: ",
    paste(missing_cols, collapse = ", "),
    "\nAvailable columns are: ",
    paste(names(win), collapse = ", ")
  )
}

setnames(win, "fraction_covered", "breadth")

win[, `:=`(
  chrom = as.character(chrom),
  sample = as.character(sample),
  start = as.integer(start),
  end = as.integer(end),
  n_overlapping_reads = as.numeric(n_overlapping_reads),
  bases_covered = as.numeric(bases_covered),
  window_length = as.numeric(window_length),
  breadth = as.numeric(breadth),
  mean_depth = as.numeric(mean_depth)
)]

win[, chrom_label := chrom_map[chrom]]
win[is.na(chrom_label), chrom_label := chrom]
win[, chrom_label := as.character(chrom_label)]

win_auto <- win[chrom_label %in% autosome_levels]

if (nrow(win_auto) == 0) {
  stop("No autosomal windows found after chromosome filtering. Check chromosome names.")
}

win_pos <- win_auto[, .(
  n_samples = uniqueN(sample),
  median_breadth = median(breadth, na.rm = TRUE),
  mean_breadth = mean(breadth, na.rm = TRUE),
  q25_breadth = as.numeric(quantile(breadth, 0.25, na.rm = TRUE)),
  q75_breadth = as.numeric(quantile(breadth, 0.75, na.rm = TRUE)),
  median_depth = median(mean_depth, na.rm = TRUE),
  mean_depth = mean(mean_depth, na.rm = TRUE),
  q25_depth = as.numeric(quantile(mean_depth, 0.25, na.rm = TRUE)),
  q75_depth = as.numeric(quantile(mean_depth, 0.75, na.rm = TRUE))
), by = .(chrom_label, start, end)]

win_pos[, chrom_label := factor(as.character(chrom_label), levels = autosome_levels)]
setorder(win_pos, chrom_label, start)

chr_lengths <- win_pos[, .(
  chr_length = max(end, na.rm = TRUE)
), by = chrom_label]

chr_lengths[, chrom_label := factor(as.character(chrom_label), levels = autosome_levels)]
setorder(chr_lengths, chrom_label)
chr_lengths[, chr_start := cumsum(shift(chr_length, fill = 0))]
chr_lengths[, chr_mid := chr_start + chr_length / 2]

win_pos <- merge(win_pos, chr_lengths, by = "chrom_label", all.x = TRUE)
win_pos[, genome_pos := chr_start + ((start + end) / 2)]
setorder(win_pos, chrom_label, start)

fwrite(
  win_pos,
  file.path(outdir, "per_25kb_window_coverage_across_samples.tsv"),
  sep = "\t"
)

low_breadth_windows <- win_pos[median_breadth < 0.95]
lowest_1_percent_depth_cutoff <- quantile(win_pos$median_depth, 0.01, na.rm = TRUE)
low_depth_windows <- win_pos[median_depth <= lowest_1_percent_depth_cutoff]

fwrite(
  low_breadth_windows,
  file.path(outdir, "low_breadth_windows_median_breadth_below_0.95.tsv"),
  sep = "\t"
)

fwrite(
  low_depth_windows,
  file.path(outdir, "lowest_1_percent_depth_windows.tsv"),
  sep = "\t"
)

p_pos_breadth <- ggplot(win_pos, aes(x = genome_pos, y = median_breadth)) +
  geom_point(size = 0.35, alpha = 0.7) +
  scale_x_continuous(
    breaks = chr_lengths$chr_mid,
    labels = as.character(chr_lengths$chrom_label),
    expand = c(0.005, 0.005)
  ) +
  scale_y_continuous(
    limits = c(0.90, 1.00),
    breaks = seq(0.90, 1.00, by = 0.02),
    oob = squish
  ) +
  labs(
    x = NULL,
    y = "Median breadth per 25 kb window",
    title = "A. Breadth across genomic position"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 9)
  )

depth_upper <- as.numeric(quantile(win_pos$median_depth, 0.995, na.rm = TRUE))
depth_upper <- ceiling(depth_upper * 10) / 10

high_depth_windows <- win_pos[median_depth > depth_upper]

fwrite(
  high_depth_windows,
  file.path(outdir, "high_depth_windows_above_axis_limit.tsv"),
  sep = "\t"
)

p_pos_depth <- ggplot(win_pos, aes(x = genome_pos, y = median_depth)) +
  geom_point(size = 0.35, alpha = 0.7) +
  scale_x_continuous(
    breaks = chr_lengths$chr_mid,
    labels = as.character(chr_lengths$chrom_label),
    expand = c(0.005, 0.005)
  ) +
  scale_y_continuous(
    limits = c(0, depth_upper),
    oob = squish
  ) +
  labs(
    x = "Chromosome",
    y = "Median depth per 25 kb window",
    title = "B. Depth across genomic position"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 9)
  )

p_pos_main <- p_pos_breadth / p_pos_depth +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Genome-wide coverage across genomic position (25 kb windows)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0)
    )
  )

ggsave(
  file.path(outdir, "SuppFig_25kb_coverage_across_genome.pdf"),
  p_pos_main,
  width = 11,
  height = 7
)

ggsave(
  file.path(outdir, "SuppFig_25kb_coverage_across_genome.png"),
  p_pos_main,
  width = 11,
  height = 7,
  dpi = 600
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    file.path(outdir, "SuppFig_25kb_coverage_across_genome.svg"),
    p_pos_main,
    width = 11,
    height = 7
  )
} else {
  message("Package 'svglite' is not installed, so SVG was not saved. PDF and PNG were saved.")
}

print(p_pos_main)
```

---

## 11. Main inputs and outputs

### Inputs

```text
SAMPLE_decon_R1.fastq.gz
SAMPLE_decon_R2.fastq.gz
Ascaris_primary_scaffolds_plus_X.fasta
Ascaris_primary_scaffolds_plus_X.fasta.fai
Ascaris_primary_scaffolds_plus_X.25kb.bed
Sample_list.txt
```

### Mapping outputs

```text
SAMPLE.bam
SAMPLE.bam.bai
SAMPLE.flagstat.txt
```

### Duplicate-marking outputs

```text
SAMPLE.markdup.bam
SAMPLE.markdup.bam.bai
SAMPLE.metrics.txt
```

### Sample-level mapping-QC outputs

```text
sample_mapping_stats.txt
overall_mapping_stats.txt
SAMPLE.stats
Table_S1_mapping_QC.csv
failure_samples.csv
watch_samples.csv
mapping_QC_excluding_failures.csv
mapping_QC_summary.csv
Fig_QC_depth_breadth.pdf
Fig_QC_depth_breadth.png
Fig_QC_depth_breadth.tiff
```

### Window-level outputs

```text
Ascaris_primary_scaffolds_plus_X.genome.txt
Ascaris_primary_scaffolds_plus_X.25kb.bed
SAMPLE.25kb.cov_stats.tsv
New_25kb_cov_stats.tsv
per_chromosome_25kb_breadth_summary.tsv
per_chromosome_25kb_depth_summary.tsv
per_25kb_window_coverage_across_samples.tsv
low_breadth_windows_median_breadth_below_0.95.tsv
lowest_1_percent_depth_windows.tsv
high_depth_windows_above_axis_limit.tsv
```

### Figure outputs

```text
Figure2A_breadth_median_IQR.pdf/png
Figure2B_depth_median_IQR.pdf/png
Figure2_main_breadth_and_depth.pdf/png
SuppFig_25kb_coverage_across_genome.pdf/png
SuppFig_25kb_coverage_across_genome.svg
```

---

### 12. Software versions

```text
BWA-MEM: v0.7.19-r1273
GATK: v4.6.2.0
Picard: v3.4.0
samtools: v1.22
```

---

## Outputs and QC checks

Indexed duplicate-marked BAMs; sample-level and 25-kb QC. Consult recorded logs and intermediate tables to distinguish expected settings from actual execution success.
