# 01 — Raw-read acquisition, trimming and host-read depletion

## Purpose

Raw-read acquisition, trimming and host-read depletion. This file archives the commands and analysis notes from the thesis/paper workflow.

## Inputs

Paired-end FASTQ, SRA run accessions and manually curated sample identifiers. Source filenames and absolute HPC paths are retained in the original script blocks below.

## Workflow

See the original annotated analysis sequence below. Distinct script blocks are independent jobs/R sessions, **not one executable script**.

## Annotated scripts

## Overview

The study included 28 newly generated paired-end whole-genome sequencing datasets and 128 publicly available datasets downloaded from the NCBI Sequence Read Archive, giving 156 samples in total.

The preprocessing workflow was initially developed and tested using the newly generated datasets. The same procedures and parameters were subsequently applied to the publicly available datasets after download. Samples were processed in separate batches as they became available, but all datasets were taken through the same main stages:

1. acquisition or incorporation of paired-end FASTQ files;
2. confirmation of paired-read availability and read counts;
3. adapter and quality trimming using fastp;
4. removal of possible human- and pig-derived reads using Deacon;
5. quality assessment using FastQC;
6. aggregation of quality-control reports using MultiQC.

The 28 newly generated datasets and the 128 downloaded datasets were not concatenated into single FASTQ files. They remained separate samples and were combined at the dataset and sample-list level for subsequent mapping and variant calling.

---

## 1. Download of publicly available SRA datasets

Publicly available whole-genome sequencing datasets were downloaded using the SRA Toolkit. Run accessions were supplied in a text file containing one accession per line.

The initial download procedure used `prefetch` to retrieve each SRA record, followed by `fasterq-dump` to produce paired-end FASTQ files. The FASTQ files were subsequently compressed using `pigz`.

```bash
#!/bin/bash

# Text file containing one SRA run accession per line.
IDFILE="accessions.txt"

# Download each SRA run, convert it into paired FASTQ files and
# compress the resulting read files.
while read ACC; do

  # Download the SRA record.
  prefetch "$ACC"

  # Convert the SRA record into separate paired-end FASTQ files.
  fasterq-dump --split-files "${ACC}.sra"

  # Compress the paired FASTQ outputs.
  pigz "${ACC}_1.fastq"
  pigz "${ACC}_2.fastq"

done < "$IDFILE"
```

Expected outputs for each paired-end SRA run were:

```text
ACCESSION_1.fastq.gz
ACCESSION_2.fastq.gz
```

Additional download and validation jobs were used where individual accessions were incomplete or absent after the initial download.

---

## 2. Recovery and validation of missing SRA datasets

Thirteen accessions that were missing or incomplete following the initial download were processed using a more detailed recovery script.

This script:

* selected one accession per SLURM array task;
* retried failed downloads;
* directed the SRA cache and temporary files to the parallel filesystem;
* retrieved the SRA record using `prefetch`;
* located the downloaded `.sra` file;
* generated paired FASTQ files using `fasterq-dump`;
* compressed the files using `pigz`;
* checked that both files existed;
* tested gzip integrity;
* confirmed that R1 and R2 contained equal numbers of reads.

```bash
#!/bin/bash
#SBATCH --job-name=sra_missing
#SBATCH --partition=shared
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=90:00:00
#SBATCH --array=1-13%4
#SBATCH --output=/parallel_scratch/lw01165/sra_files/scripts/sra_missing_%A_%a.out
#SBATCH --error=/parallel_scratch/lw01165/sra_files/scripts/sra_missing_%A_%a.err

set -euo pipefail

# ------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------

# List of SRA accessions requiring recovery.
IDFILE="/parallel_scratch/lw01165/missing_ids.txt"

# Final output directory for paired FASTQ files.
OUTDIR="/parallel_scratch/lw01165/sra_files"

# SRA download cache on the parallel filesystem.
CACHE="/parallel_scratch/lw01165/sra_cache"

# Temporary working directory.
TMPBASE="/parallel_scratch/lw01165/sra_tmp_${SLURM_JOB_ID:-local}"

# Number of threads and permitted retry attempts.
THREADS="${SLURM_CPUS_PER_TASK:-8}"
RETRIES=4

mkdir -p "$OUTDIR" "$CACHE" "$TMPBASE"


# ------------------------------------------------------------------------
# SOFTWARE ENVIRONMENT
# ------------------------------------------------------------------------

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate SRA_tools

# Confirm that all required programs are available.
command -v prefetch >/dev/null || {
  echo "prefetch not found"
  exit 1
}

command -v fasterq-dump >/dev/null || {
  echo "fasterq-dump not found"
  exit 1
}

command -v pigz >/dev/null || {
  echo "pigz not found"
  exit 1
}

command -v vdb-config >/dev/null || {
  echo "vdb-config not found"
  exit 1
}


# ------------------------------------------------------------------------
# CONFIGURE SRA CACHE
# ------------------------------------------------------------------------

# Direct the SRA Toolkit cache to the parallel filesystem rather
# than the smaller default user storage area.
vdb-config -s "/repository/user/main/public/root=${CACHE}"


# ------------------------------------------------------------------------
# SELECT ACCESSION
# ------------------------------------------------------------------------

# Confirm that the accession list exists.
[[ -s "$IDFILE" ]] || {
  echo "ID list missing: $IDFILE"
  exit 1
}

# Select the accession corresponding to the current SLURM array task.
ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$IDFILE" \
  | tr -d '\r' \
  | sed 's/[[:space:]]*$//')

# Exit without error if no accession is assigned to the task.
[[ -n "${ACC:-}" ]] || {
  echo "No accession for task ${SLURM_ARRAY_TASK_ID}"
  exit 0
}

echo "==> $(date) START ${ACC} (task ${SLURM_ARRAY_TASK_ID})"

# Create a sample-specific temporary directory.
TMP="${TMPBASE}/${ACC}"
mkdir -p "$TMP"

# Expected paired FASTQ outputs.
R1="${OUTDIR}/${ACC}_1.fastq.gz"
R2="${OUTDIR}/${ACC}_2.fastq.gz"


# ------------------------------------------------------------------------
# HELPER FUNCTIONS
# ------------------------------------------------------------------------

# Retry a command using exponential backoff.
backoff() {
  local max=$1
  shift
  local n=1

  until "$@"; do
    ((n >= max)) && return 1
    sleep $((2**(n-1)))
    ((n++))
  done
}


# Confirm that both paired files:
#   - exist and are non-empty;
#   - pass gzip integrity checks;
#   - contain equal numbers of FASTQ records.
validate_pair() {
  local a1="$1"
  local a2="$2"

  [[ -s "$a1" && -s "$a2" ]] || return 13

  gzip -t "$a1" "$a2" || return 10

  local c1
  local c2

  c1=$(pigz -cd "$a1" | wc -l)
  c1=$((c1 / 4))

  c2=$(pigz -cd "$a2" | wc -l)
  c2=$((c2 / 4))

  [[ "$c1" -eq "$c2" ]] || return 12

  echo "R1=$c1 R2=$c2"
}


# ------------------------------------------------------------------------
# SKIP ACCESSIONS ALREADY DOWNLOADED AND VALIDATED
# ------------------------------------------------------------------------

if validate_pair "$R1" "$R2" >/dev/null 2>&1; then
  echo "✅ ${ACC} already present and valid — skipping"
  exit 0
fi


# ------------------------------------------------------------------------
# DOWNLOAD SRA RECORD
# ------------------------------------------------------------------------

# Use Aspera transport where available.
PREF_OPTS=(-p -O "$CACHE")

command -v ascp >/dev/null 2>&1 && PREF_OPTS+=(--transport ascp)

echo "📥 prefetch ${ACC}"

if ! backoff "$RETRIES" prefetch "${PREF_OPTS[@]}" "$ACC"; then
  echo "❌ prefetch failed for ${ACC}"
  exit 1
fi


# Locate the downloaded SRA file in the configured cache.
SRA_PATH="$(find "$CACHE" -type f -name "${ACC}.sra" -print -quit)"

[[ -s "$SRA_PATH" ]] || {
  echo "❌ ${ACC}.sra missing after prefetch"
  exit 1
}


# ------------------------------------------------------------------------
# CONVERT SRA TO PAIRED FASTQ
# ------------------------------------------------------------------------

export TMPDIR="$TMP"

echo "🔄 fasterq-dump ${ACC} (threads=$THREADS)"

if ! backoff "$RETRIES" \
  fasterq-dump \
    -e "$THREADS" \
    --split-files \
    -p \
    -t "$TMP" \
    -O "$OUTDIR" \
    "$SRA_PATH"; then

  echo "❌ fasterq-dump failed for ${ACC}"
  exit 1
fi


# ------------------------------------------------------------------------
# COMPRESS FASTQ FILES
# ------------------------------------------------------------------------

echo "🗜️ pigz compress"

[[ -s "${OUTDIR}/${ACC}_1.fastq" ]] && \
  pigz -p "$THREADS" -f "${OUTDIR}/${ACC}_1.fastq"

[[ -s "${OUTDIR}/${ACC}_2.fastq" ]] && \
  pigz -p "$THREADS" -f "${OUTDIR}/${ACC}_2.fastq"


# ------------------------------------------------------------------------
# VALIDATE FINAL PAIRED FILES
# ------------------------------------------------------------------------

echo "🔎 validate ${ACC}"

if out=$(validate_pair "$R1" "$R2"); then

  echo "✅ OK ${ACC} ${out}"

  # Remove temporary files after successful validation.
  rm -rf "$TMP" 2>/dev/null || true

  exit 0

else

  echo "⚠️ validation failed for ${ACC} — keeping files for inspection"
  exit 2
fi
```

The recovery script provided a stronger validation step than the initial download loop by locating the full cached SRA path and checking the integrity and paired-read counts of the compressed outputs.

---

## 3. Incorporation and checking of newly generated reads

The 28 newly generated datasets were received as paired compressed FASTQ files. These datasets were initially used to test and establish the preprocessing workflow before the same approach was applied to the public datasets.

A read-count script was used to record the number of reads in each R1 and R2 file.

A FASTQ record contains four lines, so the number of reads was calculated as:

```text
number of decompressed lines / 4
```

```bash
#!/bin/bash
#SBATCH --job-name=check_counts_simple
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=check_counts_simple_%j.out
#SBATCH --error=check_counts_simple_%j.err

# ------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------

# Directory containing the newly generated FASTQ files.
RAW_DIR="/parallel_scratch/lw01165/02_Samples"

# List containing one sample identifier per line.
ID_LIST="${RAW_DIR}/sample.ids.txt"

# Output table.
OUT="${RAW_DIR}/raw_read_counts.tsv"


# ------------------------------------------------------------------------
# OUTPUT HEADER
# ------------------------------------------------------------------------

echo -e "Sample\tR1_reads\tR2_reads" > "$OUT"

cd "$RAW_DIR" || exit 1


# ------------------------------------------------------------------------
# SELECT DECOMPRESSION PROGRAM
# ------------------------------------------------------------------------

# Use gzip by default, but use pigz where available.
CAT="gzip -cd"

command -v pigz >/dev/null && CAT="pigz -cd"


# ------------------------------------------------------------------------
# COUNT PAIRED READS
# ------------------------------------------------------------------------

while IFS= read -r sample; do

  # Skip blank lines and comments.
  [[ -z "$sample" || "$sample" =~ ^# ]] && continue

  # Expected paired FASTQ files for each sample.
  R1="${sample}_1.fq.gz"
  R2="${sample}_2.fq.gz"

  # Count R1 records or report NA if the file is absent.
  if [[ -f "$R1" ]]; then

    c1_lines=$($CAT "$R1" | wc -l)
    c1=$((c1_lines / 4))

  else

    c1="NA"
  fi

  # Count R2 records or report NA if the file is absent.
  if [[ -f "$R2" ]]; then

    c2_lines=$($CAT "$R2" | wc -l)
    c2=$((c2_lines / 4))

  else

    c2="NA"
  fi

  # Write the sample identifier and paired-read counts.
  echo -e "${sample}\t${c1}\t${c2}" >> "$OUT"

done < "$ID_LIST"

echo "Done. Results written to $OUT"
```

This generated:

```text
raw_read_counts.tsv
```

The newly generated files used the suffixes:

```text
_1.fq.gz
_2.fq.gz
```

whereas the downloaded SRA files used:

```text
_1.fastq.gz
_2.fastq.gz
```

The difference reflects their original file naming rather than a difference in the sequence data.

---

## 4. Trimming and quality filtering using fastp

The paired reads were trimmed and quality-filtered using fastp.

The main parameters were:

```text
paired-end adapter detection: enabled
5′ trimming from R1:          20 bp
5′ trimming from R2:          20 bp
minimum retained length:      50 bp
```

The workflow was initially tested using the newly generated samples and was then applied to the downloaded SRA datasets in batches.

The relevant sample list was changed between batches, but the same fastp commands and filtering parameters were used.

```bash
#!/bin/bash
#SBATCH --job-name=fastp_batched
#SBATCH --cpus-per-task=24
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/parallel_scratch/lw01165/sra_files/logs/fastp_%A_%a.out
#SBATCH --error=/parallel_scratch/lw01165/sra_files/logs/fastp_%A_%a.err
#SBATCH --requeue

set -euo pipefail


# ------------------------------------------------------------------------
# PATHS AND INPUTS
# ------------------------------------------------------------------------

BASEDIR="/parallel_scratch/lw01165/sra_files"

# One sample identifier per line.
#
# The sample list was changed as required when processing the newly
# generated and publicly downloaded sample batches.
SAMPLE_LIST="${BASEDIR}/ids.clean.txt"

# Confirmed output directory containing the fastp-processed reads.
OUTDIR="${BASEDIR}/samples_fastp"

LOGDIR="${BASEDIR}/logs"
STATUSDIR="${OUTDIR}/status"

mkdir -p "$OUTDIR" "$LOGDIR" "$STATUSDIR"


# Input file naming pattern.
R1_SUFFIX="_1.fastq.gz"
R2_SUFFIX="_2.fastq.gz"


# ------------------------------------------------------------------------
# SOFTWARE ENVIRONMENT
# ------------------------------------------------------------------------

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate fast_env


# Use the parallel filesystem for temporary outputs.
export TMPDIR="/parallel_scratch/lw01165/tmp_fastp/${SLURM_JOB_ID:-local}"
mkdir -p "$TMPDIR"


# ------------------------------------------------------------------------
# FASTP PARAMETERS
# ------------------------------------------------------------------------

THREADS=${SLURM_CPUS_PER_TASK:-24}

# Limit the number of threads used by pigz during output compression.
PIGZ_THREADS=6

# Discard reads shorter than 50 bp after trimming.
MIN_LEN=50

# Remove the first 20 bp from the 5′ end of both reads.
TRIM_FRONT=20

# Process five samples per SLURM array task.
SAMPLES_PER_TASK=5


# ------------------------------------------------------------------------
# ARRAY CHUNK CALCULATION
# ------------------------------------------------------------------------

# Require the script to be submitted as a SLURM array job.
: "${SLURM_ARRAY_TASK_ID:?Run this script as a SLURM array job}"

TASK_ID="$SLURM_ARRAY_TASK_ID"

start=$(( (TASK_ID - 1) * SAMPLES_PER_TASK + 1 ))
end=$(( TASK_ID * SAMPLES_PER_TASK ))

echo "[$(date)] Job ${SLURM_JOB_ID:-noslurm} Task $TASK_ID processing lines $start..$end from $SAMPLE_LIST"
echo "Threads: fastp=$THREADS pigz=$PIGZ_THREADS"


# ------------------------------------------------------------------------
# PER-SAMPLE PROCESSING FUNCTION
# ------------------------------------------------------------------------

run_one_sample () {

  local sample="$1"

  # Paired input FASTQ files.
  local r1="${BASEDIR}/${sample}${R1_SUFFIX}"
  local r2="${BASEDIR}/${sample}${R2_SUFFIX}"

  # Final paired fastp outputs.
  local out_r1="${OUTDIR}/${sample}_trimmed_R1.fastq.gz"
  local out_r2="${OUTDIR}/${sample}_trimmed_R2.fastq.gz"

  # Per-sample fastp HTML report.
  local out_html="${OUTDIR}/${sample}_fastp.html"

  # Completion and failure flags.
  local doneflag="${STATUSDIR}/${sample}.done"
  local failflag="${STATUSDIR}/${sample}.fail"


  # ----------------------------------------------------------------------
  # RESUME CHECK
  # ----------------------------------------------------------------------

  # Skip samples already marked as completed.
  if [[ -f "$doneflag" ]]; then
    echo "✅ [$sample] .done exists — skipping."
    return 0
  fi

  # Accept existing paired outputs as evidence of completion.
  if [[ -f "$out_r1" && -f "$out_r2" ]]; then
    echo "✅ [$sample] Outputs exist — marking done and skipping."
    touch "$doneflag"
    return 0
  fi


  # ----------------------------------------------------------------------
  # INPUT CHECK
  # ----------------------------------------------------------------------

  if [[ ! -f "$r1" || ! -f "$r2" ]]; then
    echo "❌ [$sample] Missing input: $r1 or $r2"
    touch "$failflag"
    return 0
  fi


  # ----------------------------------------------------------------------
  # SAMPLE-SPECIFIC TEMPORARY DIRECTORY
  # ----------------------------------------------------------------------

  local SCRATCH="${TMPDIR}/${sample}"
  mkdir -p "$SCRATCH"

  echo "🚀 [$sample] Starting fastp on $(hostname) at $(date)"


  # ----------------------------------------------------------------------
  # RUN FASTP
  # ----------------------------------------------------------------------

  set +e

  fastp \
    -i "$r1" \
    -I "$r2" \
    -o "${SCRATCH}/${sample}_trimmed_R1.fastq" \
    -O "${SCRATCH}/${sample}_trimmed_R2.fastq" \
    --thread "$THREADS" \
    --html "${SCRATCH}/${sample}_fastp.html" \
    --detect_adapter_for_pe \
    --trim_front1 "$TRIM_FRONT" \
    --trim_front2 "$TRIM_FRONT" \
    --length_required "$MIN_LEN"

  local status=$?

  set -e


  # Record failed samples without terminating the entire array task.
  if [[ $status -ne 0 ]]; then
    echo "❌ [$sample] fastp failed with exit code $status"
    touch "$failflag"
    rm -rf "$SCRATCH"
    return 0
  fi


  # ----------------------------------------------------------------------
  # COMPRESS OUTPUT FILES
  # ----------------------------------------------------------------------

  echo "📦 [$sample] Compressing to OUTDIR (pigz threads=$PIGZ_THREADS)"

  # Write temporary compressed outputs before renaming them into place.
  # This prevents incomplete files from being mistaken for completed outputs.
  local tmp_out_r1="${out_r1}.tmp.${SLURM_JOB_ID:-local}"
  local tmp_out_r2="${out_r2}.tmp.${SLURM_JOB_ID:-local}"

  pigz -p "$PIGZ_THREADS" \
    -c "${SCRATCH}/${sample}_trimmed_R1.fastq" \
    > "$tmp_out_r1"

  pigz -p "$PIGZ_THREADS" \
    -c "${SCRATCH}/${sample}_trimmed_R2.fastq" \
    > "$tmp_out_r2"

  mv -f "$tmp_out_r1" "$out_r1"
  mv -f "$tmp_out_r2" "$out_r2"


  # ----------------------------------------------------------------------
  # RETAIN FASTP REPORT
  # ----------------------------------------------------------------------

  cp -f "${SCRATCH}/${sample}_fastp.html" "$out_html"


  # ----------------------------------------------------------------------
  # MARK SAMPLE AS COMPLETE
  # ----------------------------------------------------------------------

  touch "$doneflag"
  rm -rf "$SCRATCH"

  echo "🏁 [$sample] Finished at $(date)"
}


# ------------------------------------------------------------------------
# MAIN BATCH LOOP
# ------------------------------------------------------------------------

if [[ ! -f "$SAMPLE_LIST" ]]; then
  echo "❌ Sample list not found: $SAMPLE_LIST"
  exit 1
fi

for i in $(seq "$start" "$end"); do

  sample="$(sed -n "${i}p" "$SAMPLE_LIST" \
    | tr -d '\r' || true)"

  [[ -z "${sample:-}" ]] && continue

  # Produce a separate searchable log for each sample.
  sample_log="${LOGDIR}/${sample}_fastp_${SLURM_JOB_ID:-noslurm}.log"

  {
    echo "=============================="
    echo "Sample: $sample"
    echo "Job: ${SLURM_JOB_ID:-noslurm} Task: ${SLURM_ARRAY_TASK_ID:-na}"
    echo "Start: $(date)"
    echo "=============================="

    run_one_sample "$sample"

    echo "End: $(date)"

  } 2>&1 | tee -a "$sample_log"

done

echo "✅ Task $TASK_ID complete at $(date)"
```

The main fastp outputs were:

```text
SAMPLE_trimmed_R1.fastq.gz
SAMPLE_trimmed_R2.fastq.gz
SAMPLE_fastp.html
```

The confirmed fastp output directory was:

```text
/parallel_scratch/lw01165/sra_files/samples_fastp
```

---

## 5. Construction of the combined human and pig host reference

Potential host-derived reads were removed using a combined human and pig reference.

The reference assemblies were:

```text
Homo sapiens: GCF_000001405.40
Sus scrofa:   GCF_000003025.6
```

The downloaded genomic FASTA files were concatenated into a single file:

```text
human_pig_genome.fna
```

The concatenation was equivalent to:

```bash
# Human reference:
#   Homo sapiens
#   RefSeq accession GCF_000001405.40
HUMAN_FASTA="/path/to/human_reference_genomic.fna"

# Pig reference:
#   Sus scrofa
#   RefSeq accession GCF_000003025.6
PIG_FASTA="/path/to/pig_reference_genomic.fna"

# Combine both reference FASTA files into one host reference.
cat "$HUMAN_FASTA" "$PIG_FASTA" \
  > human_pig_genome.fna
```

The precise local filenames of the downloaded human and pig FASTA files depended on the filenames supplied by NCBI. The analytical operation was concatenation of the two complete genomic FASTA files.

The large source and intermediate host-reference files were subsequently deleted when additional storage space was required. The host-depletion outputs and relevant analysis records were retained.

---

## 6. Construction of the Deacon host-reference index

A Deacon minimiser index was built from the combined human and pig reference.

```bash
#!/bin/bash
#SBATCH --job-name=deacon_index
#SBATCH --output=/parallel_scratch/lw01165/01_Reference/deacon_human_pig/deacon_index_%j.out
#SBATCH --error=/parallel_scratch/lw01165/01_Reference/deacon_human_pig/deacon_index_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

set -euo pipefail


# ------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------

WORKDIR="/parallel_scratch/lw01165/01_Reference/deacon_human_pig"

cd "$WORKDIR"

# Combined human and pig genomic FASTA.
FASTA="${WORKDIR}/human_pig_genome.fna"

# Deacon index output.
IDX_FILE="${WORKDIR}/human_pig_genome.idx"

# Temporary working directory.
TMPDIR="${WORKDIR}/tmp_${SLURM_JOB_ID}"

mkdir -p "$TMPDIR"
export TMPDIR

THREADS=${SLURM_CPUS_PER_TASK:-8}


# ------------------------------------------------------------------------
# SOFTWARE ENVIRONMENT
# ------------------------------------------------------------------------

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate deacon


# ------------------------------------------------------------------------
# INPUT CHECK
# ------------------------------------------------------------------------

if [[ ! -s "$FASTA" ]]; then
  echo "❌ ERROR: Combined FASTA not found: $FASTA"
  exit 1
fi


# ------------------------------------------------------------------------
# BUILD DEACON INDEX
# ------------------------------------------------------------------------

echo "==> Building Deacon index"
echo "FASTA: $FASTA"
echo "Output idx: $IDX_FILE"
echo "Threads: $THREADS"

deacon index build "$FASTA" > "$IDX_FILE"

echo "✅ Finished building index: $(date)"
```

The main index output was:

```text
human_pig_genome.idx
```

---

## 7. Removal of host-derived reads using Deacon

The fastp-processed paired reads were filtered against the combined human and pig Deacon index.

The workflow was initially applied to the newly generated samples and subsequently to the publicly available datasets. The same file structure and filtering command were used, with the relevant sample list supplied for each processing batch.

```bash
#!/bin/bash
#SBATCH --job-name=deacon_deplete
#SBATCH --output=/parallel_scratch/lw01165/sra_files/logs/deacon_deplete_%j.out
#SBATCH --error=/parallel_scratch/lw01165/sra_files/logs/deacon_deplete_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=36:00:00

set -euo pipefail


# ------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------

# Combined human and pig Deacon index.
REFDIR="/parallel_scratch/lw01165/01_Reference/deacon_human_pig"
IDX="${REFDIR}/human_pig_genome.idx"

BASEDIR="/parallel_scratch/lw01165/sra_files"

# fastp-processed input reads.
IN_DIR="${BASEDIR}/samples_fastp"

# Host-depleted output reads.
OUT_DIR="${BASEDIR}/samples_deacon"

# Per-sample logs and status files.
LOGDIR="${BASEDIR}/decon_logs"
STATUSDIR="${OUT_DIR}/status"

mkdir -p "$OUT_DIR" "$LOGDIR" "$STATUSDIR"

# Sample identifiers for the relevant processing batch.
SAMPLE_LIST="${BASEDIR}/ids.clean.txt"

THREADS=${SLURM_CPUS_PER_TASK:-16}


# ------------------------------------------------------------------------
# SOFTWARE ENVIRONMENT
# ------------------------------------------------------------------------

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate deacon


# ------------------------------------------------------------------------
# PROGRAM AND INPUT CHECKS
# ------------------------------------------------------------------------

command -v deacon >/dev/null || {
  echo "❌ deacon not found in PATH"
  exit 1
}

[[ -s "$IDX" ]] || {
  echo "❌ Index not found: $IDX"
  exit 1
}

[[ -s "$SAMPLE_LIST" ]] || {
  echo "❌ Sample list not found: $SAMPLE_LIST"
  exit 1
}


# ------------------------------------------------------------------------
# RUN INFORMATION
# ------------------------------------------------------------------------

echo "==> Host:     $(hostname)"
echo "==> Index:    $IDX"
echo "==> Samples:  $SAMPLE_LIST (n=$(wc -l < "$SAMPLE_LIST"))"
echo "==> IN_DIR:   $IN_DIR"
echo "==> OUT_DIR:  $OUT_DIR"
echo "==> Threads:  $THREADS"
echo "==> Start:    $(date)"


# ------------------------------------------------------------------------
# PER-SAMPLE DEACON FUNCTION
# ------------------------------------------------------------------------

run_one() {

  local sample="$1"

  # fastp-processed paired reads.
  local r1="${IN_DIR}/${sample}_trimmed_R1.fastq.gz"
  local r2="${IN_DIR}/${sample}_trimmed_R2.fastq.gz"

  # Host-depleted paired outputs.
  local o1="${OUT_DIR}/${sample}_decon_R1.fastq.gz"
  local o2="${OUT_DIR}/${sample}_decon_R2.fastq.gz"

  # Per-sample Deacon JSON summary.
  local js="${OUT_DIR}/${sample}_deacon_summary.json"

  # Completion and failure flags.
  local doneflag="${STATUSDIR}/${sample}.done"
  local failflag="${STATUSDIR}/${sample}.fail"


  # ----------------------------------------------------------------------
  # RESUME CHECKS
  # ----------------------------------------------------------------------

  if [[ -f "$doneflag" ]]; then
    echo "✅ [$sample] done flag exists — skipping"
    return 0
  fi


  # ----------------------------------------------------------------------
  # INPUT CHECK
  # ----------------------------------------------------------------------

  if [[ ! -s "$r1" || ! -s "$r2" ]]; then
    echo "❌ [$sample] missing input(s): $r1 or $r2"
    touch "$failflag"
    return 0
  fi


  # ----------------------------------------------------------------------
  # EXISTING OUTPUT CHECK
  # ----------------------------------------------------------------------

  if [[ -s "$o1" && -s "$o2" ]]; then
    echo "✅ [$sample] outputs exist — marking done"
    touch "$doneflag"
    return 0
  fi


  # ----------------------------------------------------------------------
  # RUN DEACON
  # ----------------------------------------------------------------------

  echo "🚀 [$sample] running deacon filter..."

  set +e

  deacon filter \
    -d "$IDX" \
    "$r1" \
    "$r2" \
    -o "$o1" \
    -O "$o2" \
    -s "$js" \
    --threads "$THREADS"

  rc=$?

  set -e


  # Record failed samples without terminating the full batch.
  if [[ $rc -ne 0 ]]; then
    echo "❌ [$sample] deacon failed (exit $rc)"
    touch "$failflag"
    return 0
  fi

  touch "$doneflag"

  echo "🏁 [$sample] finished"
}


# ------------------------------------------------------------------------
# MAIN SAMPLE LOOP
# ------------------------------------------------------------------------

while IFS= read -r raw; do

  # Remove trailing whitespace.
  sample="$(echo "$raw" | sed 's/[[:space:]]*$//')"

  # Skip blank lines and comments.
  [[ -z "$sample" ]] && continue
  [[ "$sample" =~ ^# ]] && continue

  # Produce a separate log file for each sample.
  {
    run_one "$sample"
  } 2>&1 | tee -a \
    "${LOGDIR}/${sample}_deacon_${SLURM_JOB_ID:-noslurm}.log"

done < "$SAMPLE_LIST"


# ------------------------------------------------------------------------
# COMPLETION SUMMARY
# ------------------------------------------------------------------------

echo "✅ All samples processed at $(date)"

echo "Done: $(ls -1 "$STATUSDIR"/*.done 2>/dev/null | wc -l) \
Fail: $(ls -1 "$STATUSDIR"/*.fail 2>/dev/null | wc -l)"
```

The main Deacon outputs were:

```text
SAMPLE_decon_R1.fastq.gz
SAMPLE_decon_R2.fastq.gz
SAMPLE_deacon_summary.json
```

The JSON summaries recorded the outcome of host filtering for each sample.

---

## 8. FastQC assessment

FastQC was run on different batches and processing stages as the reads became available.

The example below assessed all host-depleted FASTQ files in:

```text
/parallel_scratch/lw01165/sra_files/samples_deacon
```

Equivalent FastQC jobs were run for the newly generated samples and other relevant read batches.

```bash
#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --partition=shared
#SBATCH --output=fastqc_%j.out
#SBATCH --error=fastqc_%j.err

set -eo pipefail


# ------------------------------------------------------------------------
# SOFTWARE ENVIRONMENT
# ------------------------------------------------------------------------

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate fast_env

command -v fastqc >/dev/null || {
  echo "fastqc not found in PATH"
  exit 1
}


# ------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------

# Host-depleted paired FASTQ files.
INDIR="/parallel_scratch/lw01165/sra_files/samples_deacon/"

# FastQC output directory.
OUTDIR="/parallel_scratch/lw01165/sra_files/fastqc_reports"

mkdir -p "$OUTDIR"


# ------------------------------------------------------------------------
# RUN FASTQC
# ------------------------------------------------------------------------

cd "$INDIR"

# Include compressed or uncompressed FASTQ files using either the
# .fq or .fastq filename extension.
shopt -s nullglob

FILES=(
  *.fq.gz
  *.fastq.gz
  *.fq
  *.fastq
)

# Stop without error if no FASTQ files are found.
if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "No FASTQ files found in $INDIR"
  exit 0
fi

echo "Found ${#FILES[@]} files. Running FastQC..."

fastqc \
  -t "${SLURM_CPUS_PER_TASK:-8}" \
  -o "$OUTDIR" \
  --noextract \
  "${FILES[@]}"

echo "Done. Reports in: $OUTDIR"
```

FastQC produced an HTML report and compressed results file for each FASTQ file.

---

## 9. MultiQC aggregation

MultiQC was run after FastQC to combine the individual reports into batch-level summaries.

Multiple MultiQC jobs were run because the 28 newly generated samples and 128 public datasets were processed at different times and were stored in different directories during preprocessing.

The example below summarised the FastQC reports for one of the newly generated sequencing batches.

```bash
#!/bin/bash
#SBATCH --job-name=multiqc_fastqc
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/parallel_scratch/lw01165/sra_files/logs/multiqc_fastqc_%j.out
#SBATCH --error=/parallel_scratch/lw01165/sra_files/logs/multiqc_fastqc_%j.err

set -euo pipefail


# ------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------

# Directory containing the FastQC reports for this sample batch.
FASTQC_DIR="/parallel_scratch/lw01165/02_Samples/X204SC24112216-Z01-F003/04.fastqc"

# Directory in which the MultiQC summary was written.
OUT_DIR="${FASTQC_DIR}/05.multiqc"

mkdir -p "$OUT_DIR"


# ------------------------------------------------------------------------
# SOFTWARE ENVIRONMENT
# ------------------------------------------------------------------------

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate fast_env

command -v multiqc >/dev/null || {
  echo "❌ multiqc not found in PATH"
  exit 1
}


# ------------------------------------------------------------------------
# RUN INFORMATION
# ------------------------------------------------------------------------

echo "==> Host:    $(hostname)"
echo "==> Input:   $FASTQC_DIR"
echo "==> Output:  $OUT_DIR"
echo "==> Start:   $(date)"


# ------------------------------------------------------------------------
# RUN MULTIQC
# ------------------------------------------------------------------------

multiqc "$FASTQC_DIR" -o "$OUT_DIR"


# ------------------------------------------------------------------------
# COMPLETION MESSAGE
# ------------------------------------------------------------------------

echo "✅ MultiQC finished at $(date)"
echo "Report available:"
echo "   ${OUT_DIR}/multiqc_report.html"
```

Equivalent commands were run with the relevant FastQC input directory for each additional batch.

---

## 10. Main preprocessing outputs

The main files passed into the next stage were the host-depleted paired FASTQ files:

```text
SAMPLE_decon_R1.fastq.gz
SAMPLE_decon_R2.fastq.gz
```

Supporting QC outputs included:

```text
raw_read_counts.tsv
SAMPLE_fastp.html
SAMPLE_deacon_summary.json
FastQC HTML reports
FastQC ZIP archives
multiqc_report.html
multiqc_data/
```

---

## Outputs and QC checks

FASTQ and Deacon/FastQC/MultiQC quality reports. Consult recorded logs and intermediate tables to distinguish expected settings from actual execution success.

## Provenance and limitations

- 128 downloaded SRA + 28 newly generated WGS samples; original code includes recovery steps.
- Sample identifiers and host/country metadata were manually curated in `sample_QC_list.csv`; it is not generated by an algorithm in these files.
- Commands, versions, seeds and paths are retained where present; annotations are not evidence that all jobs have been freshly rerun.
- For the discovery and correction of the PIXY all-sites input issue, see [`99_PIXY_allsites_correction.md`](99_PIXY_allsites_correction.md).
