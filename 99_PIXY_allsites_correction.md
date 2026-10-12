# 99 — PIXY all-sites mappability correction (separate audit)

## Purpose

PIXY all-sites mappability correction (separate audit). This file archives the commands and analysis notes from the thesis/paper workflow.

## Inputs

Historical 146-sample all-sites VCF, GenMap BED and curated SNP VCF. Source filenames and absolute HPC paths are retained in the original script blocks below.

## Workflow

See the original annotated analysis sequence below. Distinct script blocks are independent jobs/R sessions, **not one executable script**.

## Annotated scripts

> **Important:** this amendment records construction and QC of corrected input only. It did **not** generate the historical thesis PIXY figures/statistics. Do not imply any post-correction reanalysis took place.


## Purpose

This document records the correction made to the PIXY input preparation pipeline used for the *Ascaris* whole-genome analysis.

During final checking of the analysis pipeline, it was identified that the original all-sites VCF used for PIXY analyses had **not been restricted to uniquely mappable regions**, although the filtered SNP dataset had been restricted to these regions.

The purpose of this correction was therefore to reconstruct a PIXY-compatible VCF in which:

1. invariant sites were restricted to uniquely mappable regions;
2. invariant sites had no more than 10% missing genotype data;
3. the existing final filtered biallelic SNP dataset was retained;
4. invariant and SNP records were combined into a single sorted and indexed VCF.

No remapping, variant calling, genotype calling, or complete population-genomic reanalysis was performed.

---

# Input files

## Original all-sites VCF

```text
/parallel_scratch/lw01165/05_analysis/05c_allsites_146/Ascaris.autosomes19.ALLSITES.146.vcf.gz
```

This contained:

```text
210,854,819 records
146 samples
19 autosomes
```

## Existing final SNP VCF

```text
/parallel_scratch/lw01165/05_analysis/05b_146_data/Ascaris.autosomes.SNPs.146.maxmiss0.90.vcf.gz
```

Initial audit:

```text
13,484,735 records
13,484,735 biallelic SNPs
0 multiallelic records
0 non-SNP/non-biallelic records
146 samples
```

## Uniquely mappable regions

```text
/parallel_scratch/lw01165/01_Reference/mappability_unique.merged.bed
```

BED summary:

```text
192,116 intervals
246,625,982 uniquely mappable bp
```

---

# Working directory

```text
/parallel_scratch/lw01165/05_analysis/050_PIXY_allsites_redo/
```

Main directories generated:

```text
01_inputs/
02_uniquemap/
03_qc/
04_pixy_input/
tmp/
```

Current QC files:

```text
00_input_audit.txt
01_uniquemap_contig_counts.tsv
01_uniquemap_retention.tsv
02_record_type_check.tsv
03_final_pixy_contig_counts.tsv
03_final_pixy_record_types.tsv
03_final_pixy_summary.tsv
allsites_contig_counts.tsv
allsites_ref_alt_summary.tsv
mappability_bed_contigs.txt
samples_extra_in_allsites.txt
samples_missing_from_allsites.txt
```

---

# Step 00 — Input audit

## Purpose

The first step checked that the all-sites and SNP VCFs:

- contained the same 146 samples;
- contained samples in the same order;
- were indexed;
- contained the expected autosomal records;
- contained the expected final SNP dataset;
- were compatible with the uniquely mappable BED file.

The final SNP VCF was used to define the canonical sample order.

## Key checks obtained

```text
Final SNP VCF samples: 146
All-sites VCF samples: 146

Missing from all-sites VCF: 0
Extra in all-sites VCF: 0

Sample identities and sample order already match
```

Record counts:

```text
All-sites VCF: 210,854,819
Final SNP VCF: 13,484,735
```

Final SNP VCF content:

```text
Total records:          13,484,735
Biallelic SNP records:  13,484,735
Non-biallelic/non-SNP:  0
Multiallelic records:   0
```

## Audit script

```bash
#!/bin/bash
#SBATCH --job-name=pixy00_audit
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/parallel_scratch/lw01165/05_analysis/050_PIXY_allsites_redo/00_audit_inputs.%j.out
#SBATCH --error=/parallel_scratch/lw01165/05_analysis/050_PIXY_allsites_redo/00_audit_inputs.%j.err

set +u

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env

set -eo pipefail

BASE=/parallel_scratch/lw01165/05_analysis/050_PIXY_allsites_redo

ALLSITES=/parallel_scratch/lw01165/05_analysis/05c_allsites_146/Ascaris.autosomes19.ALLSITES.146.vcf.gz

FINAL_SNP=/parallel_scratch/lw01165/05_analysis/05b_146_data/Ascaris.autosomes.SNPs.146.maxmiss0.90.vcf.gz

MAPPABILITY_BED=/parallel_scratch/lw01165/01_Reference/mappability_unique.merged.bed

INPUT_DIR=${BASE}/01_inputs
QC_DIR=${BASE}/03_qc

mkdir -p "${INPUT_DIR}" "${QC_DIR}"

SAMPLE_LIST=${INPUT_DIR}/samples_146.canonical.txt

echo "============================================================"
echo "PIXY redo: input audit"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "============================================================"

echo
echo "Software versions:"
bcftools --version | head -n 1
vcftools --version 2>&1 | head -n 1 || true

for file in \
    "${ALLSITES}" \
    "${FINAL_SNP}" \
    "${MAPPABILITY_BED}"
do
    if [[ ! -s "${file}" ]]; then
        echo "ERROR: Missing or empty input: ${file}"
        exit 1
    fi
done

echo
echo "Inputs:"
ls -lh "${MAPPABILITY_BED}" "${FINAL_SNP}" "${ALLSITES}"

bcftools query -l "${FINAL_SNP}" \
    > "${INPUT_DIR}/samples_final_snp.txt"

bcftools query -l "${ALLSITES}" \
    > "${INPUT_DIR}/samples_allsites.txt"

N_SNP_SAMPLES=$(wc -l < "${INPUT_DIR}/samples_final_snp.txt")
N_ALL_SAMPLES=$(wc -l < "${INPUT_DIR}/samples_allsites.txt")

echo
echo "Sample counts:"
echo "Final SNP VCF: ${N_SNP_SAMPLES}"
echo "All-sites VCF: ${N_ALL_SAMPLES}"

if [[ "${N_SNP_SAMPLES}" -ne 146 ]] || \
   [[ "${N_ALL_SAMPLES}" -ne 146 ]]; then
    echo "ERROR: Expected 146 samples"
    exit 1
fi

cp "${INPUT_DIR}/samples_final_snp.txt" "${SAMPLE_LIST}"

echo
echo "Sample identity comparison:"

MISSING=$(
    comm -23 \
        <(sort "${INPUT_DIR}/samples_final_snp.txt") \
        <(sort "${INPUT_DIR}/samples_allsites.txt") |
    wc -l
)

EXTRA=$(
    comm -13 \
        <(sort "${INPUT_DIR}/samples_final_snp.txt") \
        <(sort "${INPUT_DIR}/samples_allsites.txt") |
    wc -l
)

echo "Missing from all-sites VCF: ${MISSING}"
echo "Extra in all-sites VCF:     ${EXTRA}"

if diff -q \
    "${INPUT_DIR}/samples_final_snp.txt" \
    "${INPUT_DIR}/samples_allsites.txt" \
    >/dev/null
then
    echo "Sample identities and sample order already match"
else
    echo "NOTE: identities match but sample order differs"
fi

ALL_RECORDS=$(bcftools index -n "${ALLSITES}")
SNP_RECORDS=$(bcftools index -n "${FINAL_SNP}")

echo
echo "Indexed record counts:"
echo "All-sites VCF: ${ALL_RECORDS}"
echo "Final SNP VCF: ${SNP_RECORDS}"

BIALLELIC_SNPS=$(
    bcftools view \
        -m2 \
        -M2 \
        -v snps \
        -H \
        "${FINAL_SNP}" |
    wc -l
)

MULTI=$(
    bcftools view \
        -m3 \
        -H \
        "${FINAL_SNP}" |
    wc -l
)

NON_BIALLELIC=$((SNP_RECORDS - BIALLELIC_SNPS))

echo
echo "Final SNP VCF content:"
echo "Total records:          ${SNP_RECORDS}"
echo "Biallelic SNP records:  ${BIALLELIC_SNPS}"
echo "Non-biallelic/non-SNP:  ${NON_BIALLELIC}"
echo "Multiallelic records:   ${MULTI}"

BED_INTERVALS=$(wc -l < "${MAPPABILITY_BED}")

BED_BP=$(
    awk '
      {
        total += ($3 - $2)
      }
      END {
        print total
      }
    ' "${MAPPABILITY_BED}"
)

echo
echo "Uniquely mappable BED summary:"
echo "Intervals: ${BED_INTERVALS}"
echo "Total uniquely mappable bp: ${BED_BP}"

echo
echo "BED contigs:"
cut -f1 "${MAPPABILITY_BED}" |
sort -u

echo
echo "All-sites VCF contigs and record counts:"
bcftools index -s "${ALLSITES}"

echo "============================================================"
echo "Audit completed: $(date)"
echo "============================================================"
```

---

# Step 01 — Restrict the all-sites VCF to uniquely mappable regions

## Purpose

This step corrects the main preprocessing inconsistency.

The original all-sites VCF was intersected with:

```text
/parallel_scratch/lw01165/01_Reference/mappability_unique.merged.bed
```

The sample set and sample order were also explicitly retained.

No SNP hard filtering was carried out at this stage.

## Script

```bash
#!/bin/bash
#SBATCH --job-name=pixy01_umap
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --time=2-00:00:00
#SBATCH --output=/parallel_scratch/lw01165/05_analysis/050_PIXY_allsites_redo/01_mask_allsites_uniquemap.%j.out
#SBATCH --error=/parallel_scratch/lw01165/05_analysis/050_PIXY_allsites_redo/01_mask_allsites_uniquemap.%j.err

set +u

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env

set -eo pipefail

BASE=/parallel_scratch/lw01165/05_analysis/050_PIXY_allsites_redo

IN_VCF=/parallel_scratch/lw01165/05_analysis/05c_allsites_146/Ascaris.autosomes19.ALLSITES.146.vcf.gz

MAPPABILITY_BED=/parallel_scratch/lw01165/01_Reference/mappability_unique.merged.bed

SAMPLE_LIST=${BASE}/01_inputs/samples_146.canonical.txt

OUTDIR=${BASE}/02_uniquemap
QC_DIR=${BASE}/03_qc
TMPDIR=${BASE}/tmp/mask_${SLURM_JOB_ID}

OUT_VCF=${OUTDIR}/Ascaris.autosomes19.ALLSITES.146.uniquemap.vcf.gz
TMP_VCF=${TMPDIR}/Ascaris.autosomes19.ALLSITES.146.uniquemap.tmp.vcf.gz

mkdir -p \
  "${OUTDIR}" \
  "${QC_DIR}" \
  "${TMPDIR}"

trap 'rm -rf "${TMPDIR}"' EXIT

echo "============================================================"
echo "Restricting all-sites VCF to uniquely mappable regions"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "============================================================"

echo "Input VCF:       ${IN_VCF}"
echo "Mappability BED: ${MAPPABILITY_BED}"
echo "Sample order:    ${SAMPLE_LIST}"
echo "Output VCF:      ${OUT_VCF}"

for file in \
  "${IN_VCF}" \
  "${MAPPABILITY_BED}" \
  "${SAMPLE_LIST}"
do

  if [[ ! -s "${file}" ]]; then
    echo "ERROR: Missing or empty input: ${file}"
    exit 1
  fi

done

if [[ ! -s "${IN_VCF}.tbi" && ! -s "${IN_VCF}.csi" ]]; then

  echo "Input VCF index missing; creating tabix index"

  bcftools index \
    --threads "${SLURM_CPUS_PER_TASK:-16}" \
    -t \
    "${IN_VCF}"

fi

echo
echo "Applying sample-order and uniquely mappable-region restrictions..."

bcftools view \
  --threads "${SLURM_CPUS_PER_TASK:-16}" \
  -S "${SAMPLE_LIST}" \
  -T "${MAPPABILITY_BED}" \
  --targets-overlap 1 \
  -Oz \
  -o "${TMP_VCF}" \
  "${IN_VCF}"

bcftools index \
  --threads "${SLURM_CPUS_PER_TASK:-16}" \
  -t \
  "${TMP_VCF}"

bcftools query -l "${TMP_VCF}" \
  > "${TMPDIR}/masked.samples.txt"

if ! diff -q \
  "${SAMPLE_LIST}" \
  "${TMPDIR}/masked.samples.txt" \
  >/dev/null
then

  echo "ERROR: Output sample names or sample order do not match"
  diff "${SAMPLE_LIST}" "${TMPDIR}/masked.samples.txt" || true
  exit 1

fi

N_SAMPLES=$(bcftools query -l "${TMP_VCF}" | wc -l)

if [[ "${N_SAMPLES}" -ne 146 ]]; then
  echo "ERROR: Masked VCF contains ${N_SAMPLES} samples rather than 146"
  exit 1
fi

mv -f "${TMP_VCF}" "${OUT_VCF}"
mv -f "${TMP_VCF}.tbi" "${OUT_VCF}.tbi"

INPUT_RECORDS=$(bcftools index -n "${IN_VCF}")
OUTPUT_RECORDS=$(bcftools index -n "${OUT_VCF}")

PROP_RETAINED=$(
  awk \
    -v retained="${OUTPUT_RECORDS}" \
    -v original="${INPUT_RECORDS}" \
    'BEGIN {
       if (original > 0) {
         printf "%.8f", retained / original
       }
       else {
         print "NA"
       }
     }'
)

SUMMARY=${QC_DIR}/01_uniquemap_retention.tsv

{
  echo -e "metric\tvalue"
  echo -e "input_records\t${INPUT_RECORDS}"
  echo -e "uniquemap_records\t${OUTPUT_RECORDS}"
  echo -e "proportion_retained\t${PROP_RETAINED}"
  echo -e "n_samples\t${N_SAMPLES}"
} > "${SUMMARY}"

bcftools index -s "${OUT_VCF}" \
  > "${QC_DIR}/01_uniquemap_contig_counts.tsv"

echo
echo "Retention summary:"
column -t "${SUMMARY}" || cat "${SUMMARY}"

echo
echo "Per-contig record counts:"
cat "${QC_DIR}/01_uniquemap_contig_counts.tsv"

echo
echo "Final output:"
ls -lh "${OUT_VCF}" "${OUT_VCF}.tbi"

echo "============================================================"
echo "Completed: $(date)"
echo "============================================================"
```

## Result

```text
input_records        210854819
uniquemap_records    189429567
proportion_retained  0.89838861
n_samples            146
```

Therefore:

```text
21,425,252 total all-sites records
```

were outside the uniquely mappable regions.

This represents approximately:

```text
10.16%
```

of the original all-sites records.

**Important:** this is the number of total all-sites records removed, not specifically the number of invariant sites removed.

---

# Step 02 — Determine how invariant sites are represented

## Purpose

Before extracting invariant records, the representation of non-variable sites in the all-sites VCF was checked.

A representative 1 Mb region was inspected:

```text
CM024168.1:1-1000000
```

## Script

```bash
#!/bin/bash
#SBATCH --job-name=pixy02_types
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/parallel_scratch/lw01165/05_analysis/050_PIXY_allsites_redo/02_check_record_types.%j.out
#SBATCH --error=/parallel_scratch/lw01165/05_analysis/050_PIXY_allsites_redo/02_check_record_types.%j.err

set +u
source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env
set -eo pipefail

BASE=/parallel_scratch/lw01165/05_analysis/050_PIXY_allsites_redo

VCF=${BASE}/02_uniquemap/Ascaris.autosomes19.ALLSITES.146.uniquemap.vcf.gz
OUT=${BASE}/03_qc/02_record_type_check.tsv

mkdir -p "${BASE}/03_qc"

[[ -s "${VCF}" ]] || {
  echo "ERROR: Missing VCF: ${VCF}"
  exit 1
}

echo -e "category\tcount" > "${OUT}"

bcftools query \
  -r CM024168.1:1-1000000 \
  -f '%REF\t%ALT\n' \
  "${VCF}" |
awk '
  {
    ref=$1
    alt=$2

    if (alt == ".") {
      dot++
    }
    else if (alt ~ /<NON_REF>/) {
      nonref++
    }
    else if (alt ~ /,/) {
      multi++
    }
    else if (length(ref) == 1 && length(alt) == 1) {
      snp++
    }
    else {
      other++
    }

    total++
  }

  END {
    print "total\t" total
    print "ALT_dot\t" dot+0
    print "ALT_NON_REF\t" nonref+0
    print "single_ALT_SNP\t" snp+0
    print "multiallelic_ALT\t" multi+0
    print "other\t" other+0
  }
' >> "${OUT}"

column -t "${OUT}" || cat "${OUT}"
```

## Result

```text
category          count
total             417098
ALT_dot           384880
ALT_NON_REF       0
single_ALT_SNP    27417
multiallelic_ALT  2143
other             2658
```

Therefore invariant/reference-only records are represented as:

```text
ALT=.
```

The values above refer only to the 1 Mb test region and are not genome-wide counts.

---

# Step 03 — Construct the corrected PIXY VCF

## Purpose

The final corrected PIXY input was constructed from two independent branches.

### Invariant branch

From the uniquely mappable all-sites VCF:

```text
retain ALT=.
retain F_MISSING <= 0.10
```

### SNP branch

Use the existing final SNP VCF and explicitly confirm:

```text
canonical 146 samples
uniquely mappable regions
biallelic SNPs only
F_MISSING <= 0.10
```

The two branches were checked for overlapping positions and then concatenated and sorted.

## Script

```bash
#!/bin/bash
#SBATCH --job-name=pixy03_build
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --time=3-00:00:00
#SBATCH --output=/parallel_scratch/lw01165/05_analysis/050_PIXY_allsites_redo/03_build_pixy_input.%j.out
#SBATCH --error=/parallel_scratch/lw01165/05_analysis/050_PIXY_allsites_redo/03_build_pixy_input.%j.err

set +u

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env

set -eo pipefail

BASE=/parallel_scratch/lw01165/05_analysis/050_PIXY_allsites_redo

MASKED_ALLSITES=${BASE}/02_uniquemap/Ascaris.autosomes19.ALLSITES.146.uniquemap.vcf.gz

FINAL_SNP=/parallel_scratch/lw01165/05_analysis/05b_146_data/Ascaris.autosomes.SNPs.146.maxmiss0.90.vcf.gz

MAPPABILITY_BED=/parallel_scratch/lw01165/01_Reference/mappability_unique.merged.bed
SAMPLE_LIST=${BASE}/01_inputs/samples_146.canonical.txt

OUTDIR=${BASE}/04_pixy_input
QC_DIR=${BASE}/03_qc
TMPDIR=${BASE}/tmp/build_pixy_${SLURM_JOB_ID}

INVARIANT_BCF=${TMPDIR}/Ascaris.invariant.uniquemap.maxmiss0.90.bcf
SNP_BCF=${TMPDIR}/Ascaris.SNPs.uniquemap.maxmiss0.90.bcf
COMBINED_UNSORTED=${TMPDIR}/Ascaris.PIXY.unsorted.bcf

FINAL_VCF=${OUTDIR}/Ascaris.autosomes19.PIXY.146.uniquemap.maxmiss0.90.vcf.gz

mkdir -p \
  "${OUTDIR}" \
  "${QC_DIR}" \
  "${TMPDIR}"

trap 'rm -rf "${TMPDIR}"' EXIT

echo "============================================================"
echo "Constructing corrected PIXY input"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "============================================================"

echo "Masked all-sites VCF: ${MASKED_ALLSITES}"
echo "Final SNP VCF:        ${FINAL_SNP}"
echo "Final PIXY VCF:       ${FINAL_VCF}"

for file in \
  "${MASKED_ALLSITES}" \
  "${FINAL_SNP}" \
  "${MAPPABILITY_BED}" \
  "${SAMPLE_LIST}"
do

  if [[ ! -s "${file}" ]]; then
    echo "ERROR: Missing or empty input: ${file}"
    exit 1
  fi

done

command -v bcftools >/dev/null 2>&1 || {
  echo "ERROR: bcftools not available"
  exit 1
}

bcftools --version | head -n 1

# ---------------------------------------------------------------------
# Invariant branch
#
# ALT="." = invariant/reference-only sites.
# F_MISSING<=0.10 = no more than 10% missing genotype calls.
# ---------------------------------------------------------------------

echo
echo "Extracting invariant records with no more than 10% missing data..."

bcftools view \
  --threads "${SLURM_CPUS_PER_TASK:-16}" \
  -i 'ALT="." && F_MISSING<=0.10' \
  -Ob \
  -o "${INVARIANT_BCF}" \
  "${MASKED_ALLSITES}"

bcftools index \
  --threads "${SLURM_CPUS_PER_TASK:-16}" \
  "${INVARIANT_BCF}"

N_INVARIANT=$(bcftools index -n "${INVARIANT_BCF}")

echo "Retained invariant records: ${N_INVARIANT}"

# ---------------------------------------------------------------------
# SNP branch
#
# Recheck the restrictions already expected in the final SNP VCF.
# ---------------------------------------------------------------------

echo
echo "Preparing filtered SNP branch..."

bcftools view \
  --threads "${SLURM_CPUS_PER_TASK:-16}" \
  -S "${SAMPLE_LIST}" \
  -T "${MAPPABILITY_BED}" \
  --targets-overlap 1 \
  -m2 \
  -M2 \
  -v snps \
  -i 'F_MISSING<=0.10' \
  -Ob \
  -o "${SNP_BCF}" \
  "${FINAL_SNP}"

bcftools index \
  --threads "${SLURM_CPUS_PER_TASK:-16}" \
  "${SNP_BCF}"

N_SNP=$(bcftools index -n "${SNP_BCF}")

echo "Retained SNP records: ${N_SNP}"

# ---------------------------------------------------------------------
# Confirm sample identity/order.
# ---------------------------------------------------------------------

bcftools query -l "${INVARIANT_BCF}" \
  > "${TMPDIR}/invariant.samples.txt"

bcftools query -l "${SNP_BCF}" \
  > "${TMPDIR}/snp.samples.txt"

if ! diff -q \
  "${SAMPLE_LIST}" \
  "${TMPDIR}/invariant.samples.txt" \
  >/dev/null
then

  echo "ERROR: Invariant branch sample order does not match"
  exit 1

fi

if ! diff -q \
  "${SAMPLE_LIST}" \
  "${TMPDIR}/snp.samples.txt" \
  >/dev/null
then

  echo "ERROR: SNP branch sample order does not match"
  exit 1

fi

echo "Sample identities and order match in both branches"

# ---------------------------------------------------------------------
# Confirm invariant branch contains only ALT=.
# ---------------------------------------------------------------------

N_BAD_INVARIANT=$(
  bcftools view \
    -H \
    -e 'ALT="."' \
    "${INVARIANT_BCF}" |
  wc -l
)

if [[ "${N_BAD_INVARIANT}" -ne 0 ]]; then
  echo "ERROR: Invariant branch contains ${N_BAD_INVARIANT} non-invariant records"
  exit 1
fi

# ---------------------------------------------------------------------
# Ensure invariant and SNP branches do not share positions.
# ---------------------------------------------------------------------

echo
echo "Checking for overlapping positions between branches..."

bcftools query \
  -f '%CHROM\t%POS\n' \
  "${INVARIANT_BCF}" |
sort -T "${TMPDIR}" \
  > "${TMPDIR}/invariant.positions.tsv"

bcftools query \
  -f '%CHROM\t%POS\n' \
  "${SNP_BCF}" |
sort -T "${TMPDIR}" \
  > "${TMPDIR}/snp.positions.tsv"

N_OVERLAP=$(
  comm -12 \
    "${TMPDIR}/invariant.positions.tsv" \
    "${TMPDIR}/snp.positions.tsv" |
  wc -l
)

echo "Shared positions between branches: ${N_OVERLAP}"

if [[ "${N_OVERLAP}" -ne 0 ]]; then
  echo "ERROR: Invariant and SNP branches overlap at ${N_OVERLAP} positions"
  exit 1
fi

rm -f \
  "${TMPDIR}/invariant.positions.tsv" \
  "${TMPDIR}/snp.positions.tsv"

# ---------------------------------------------------------------------
# Combine and sort the two branches.
# ---------------------------------------------------------------------

echo
echo "Combining invariant and SNP branches..."

bcftools concat \
  --allow-overlaps \
  -Ob \
  -o "${COMBINED_UNSORTED}" \
  "${INVARIANT_BCF}" \
  "${SNP_BCF}"

echo
echo "Sorting the combined PIXY input..."

bcftools sort \
  --max-mem 12G \
  --temp-dir "${TMPDIR}/sort" \
  -Oz \
  -o "${FINAL_VCF}" \
  "${COMBINED_UNSORTED}"

bcftools index \
  --threads "${SLURM_CPUS_PER_TASK:-16}" \
  -t \
  "${FINAL_VCF}"

# ---------------------------------------------------------------------
# Final validation.
# ---------------------------------------------------------------------

N_FINAL=$(bcftools index -n "${FINAL_VCF}")
EXPECTED=$((N_INVARIANT + N_SNP))

echo
echo "Final record totals:"
echo "Invariant records: ${N_INVARIANT}"
echo "SNP records:       ${N_SNP}"
echo "Expected total:    ${EXPECTED}"
echo "Observed total:    ${N_FINAL}"

if [[ "${N_FINAL}" -ne "${EXPECTED}" ]]; then
  echo "ERROR: Final record count does not equal invariant + SNP totals"
  exit 1
fi

N_FINAL_SAMPLES=$(bcftools query -l "${FINAL_VCF}" | wc -l)

if [[ "${N_FINAL_SAMPLES}" -ne 146 ]]; then
  echo "ERROR: Final PIXY VCF contains ${N_FINAL_SAMPLES} samples"
  exit 1
fi

bcftools query -l "${FINAL_VCF}" \
  > "${TMPDIR}/final.samples.txt"

if ! diff -q \
  "${SAMPLE_LIST}" \
  "${TMPDIR}/final.samples.txt" \
  >/dev/null
then

  echo "ERROR: Final PIXY VCF sample order does not match"
  exit 1
fi

bcftools query \
  -f '%REF\t%ALT\n' \
  "${FINAL_VCF}" |
awk '
  {
    total++

    if ($2 == ".") {
      invariant++
    }
    else if (length($1) == 1 && length($2) == 1 && $2 !~ /,/) {
      snp++
    }
    else {
      unexpected++
    }
  }

  END {
    print "category\tcount"
    print "total\t" total
    print "invariant_ALT_dot\t" invariant+0
    print "biallelic_SNP\t" snp+0
    print "unexpected\t" unexpected+0
  }
' > "${QC_DIR}/03_final_pixy_record_types.tsv"

bcftools index -s "${FINAL_VCF}" \
  > "${QC_DIR}/03_final_pixy_contig_counts.tsv"

{
  echo -e "metric\tvalue"
  echo -e "masked_allsites_records\t$(bcftools index -n "${MASKED_ALLSITES}")"
  echo -e "retained_invariant_records\t${N_INVARIANT}"
  echo -e "retained_snp_records\t${N_SNP}"
  echo -e "final_pixy_records\t${N_FINAL}"
  echo -e "shared_positions_between_branches\t${N_OVERLAP}"
  echo -e "n_samples\t${N_FINAL_SAMPLES}"
} > "${QC_DIR}/03_final_pixy_summary.tsv"

echo
echo "Construction summary:"
column -t "${QC_DIR}/03_final_pixy_summary.tsv" ||
cat "${QC_DIR}/03_final_pixy_summary.tsv"

echo
echo "Final record types:"
column -t "${QC_DIR}/03_final_pixy_record_types.tsv" ||
cat "${QC_DIR}/03_final_pixy_record_types.tsv"

echo
echo "Final output:"
ls -lh "${FINAL_VCF}" "${FINAL_VCF}.tbi"

echo "============================================================"
echo "Corrected PIXY input completed: $(date)"
echo "============================================================"
```

---

# Final corrected PIXY input

Final VCF:

```text
/parallel_scratch/lw01165/05_analysis/050_PIXY_allsites_redo/04_pixy_input/Ascaris.autosomes19.PIXY.146.uniquemap.maxmiss0.90.vcf.gz
```

Final QC:

```text
masked_allsites_records            189429567
retained_invariant_records         156583802
retained_snp_records               13484410
final_pixy_records                 170068212
shared_positions_between_branches  0
n_samples                          146
```

Record types:

```text
total              170068212
invariant_ALT_dot  156583802
biallelic_SNP      13484410
unexpected         0
```

Thus the corrected PIXY input contains:

```text
156,583,802 invariant sites
13,484,410 biallelic SNPs
170,068,212 total records
146 samples
```

---

# SNP difference relative to the existing final SNP dataset

Original final SNP VCF:

```text
13,484,735 SNPs
```

Corrected PIXY SNP branch:

```text
13,484,410 SNPs
```

Difference:

```text
325 SNPs
```

This represents approximately:

```text
0.0024%
```

of the original final SNP dataset.

Therefore, the SNP component changed very little during reconstruction of the PIXY input.

---

# Nature of the identified issue

The error concerned the **all-sites/invariant-site component** of the PIXY input.

The original workflow had applied the uniquely mappable-region restriction to the final SNP dataset, but the same restriction had not been carried through to the all-sites dataset used for PIXY.

Consequently, the original PIXY input could include invariant/reference positions in genomic regions that had already been excluded from the SNP dataset because they were not uniquely mappable.

This meant the callable/invariant component and SNP component were not being defined over exactly the same genomic space.

The corrected workflow ensures that both components are restricted to the same uniquely mappable regions.

---

# Scale of the issue

Applying the uniquely mappable-region restriction changed the all-sites VCF from:

```text
210,854,819 records
```

to:

```text
189,429,567 records
```

Therefore:

```text
21,425,252 total records
```

were outside the uniquely mappable regions.

This is approximately:

```text
10.16% of all records
```

However, this value includes all record types present in the original all-sites VCF. It should **not** be described as 21.4 million invariant sites.

The corrected invariant branch ultimately contained:

```text
156,583,802 invariant sites
```

after both:

1. restriction to uniquely mappable regions; and
2. removal of invariant sites with more than 10% missing genotype data.

The current QC outputs do not separately quantify how many invariant records were removed **solely** by the mappability restriction. No additional analysis was performed to derive this value.

---

# Interpretation for the PhD thesis

The correction was identified during final audit of the analysis pipeline.

The thesis results were generated using the original PIXY workflow. The complete PIXY-derived analysis was not repeated during final thesis preparation.

The corrected VCF has therefore been generated and retained to:

- document the identified preprocessing inconsistency;
- demonstrate the corrected workflow;
- provide the input required for complete reanalysis before publication.

Analyses based only on the independently filtered SNP dataset are not affected by the invariant-site preprocessing issue.

The complete PIXY analysis should be repeated using the corrected input before publication.

---

# Suggested code annotation for the original thesis pipeline

```text
NOTE — FINAL PIPELINE AUDIT

During final auditing of the thesis analysis, an inconsistency was
identified in the preparation of the all-sites VCF used for PIXY.

The final SNP dataset had been restricted to uniquely mappable genomic
regions, whereas the all-sites/invariant-site component had not been
subjected to the same mappability restriction.

A corrected PIXY input was therefore constructed separately.

Corrected dataset:
    invariant sites = 156,583,802
    SNPs            = 13,484,410
    total records   = 170,068,212
    samples         = 146

The original filtered SNP VCF contained 13,484,735 SNPs. Reapplication
of the final PIXY restrictions therefore altered the SNP component by
325 SNPs (0.0024%).

The complete PIXY-derived thesis analysis was not repeated following
identification of this issue. The corrected input is retained for
reproducibility and the full PIXY analysis will be rerun using the
corrected VCF prior to publication.
```

---

# Summary

The corrected pipeline does **not** represent a new variant-calling or mapping analysis.

It is a correction to the preparation of the all-sites VCF used for PIXY.

The final corrected dataset is:

```text
156,583,802 invariant sites
13,484,410 SNPs
170,068,212 total sites
146 samples
```

The independently filtered SNP dataset changed by only:

```text
325 SNPs (0.0024%)
```

The larger preprocessing discrepancy relates to the invariant/all-sites genomic space. Approximately 10.16% of all original all-sites records lay outside the uniquely mappable intervals, although the exact number of invariant sites specifically removed by this restriction was not separately calculated.

Complete PIXY-derived analyses will be repeated using the corrected input prior to publication.

---

## Outputs and QC checks

Corrected PIXY input: 156,583,802 invariant + 13,484,410 SNP records. Consult recorded logs and intermediate tables to distinguish expected settings from actual execution success.

## Provenance and limitations

- Corrected VCF prepared and audited; downstream PIXY analysis NOT rerun for thesis.
- Sample identifiers and host/country metadata were manually curated in `sample_QC_list.csv`; it is not generated by an algorithm in these files.
- Commands, versions, seeds and paths are retained where present; annotations are not evidence that all jobs have been freshly rerun.
- For the discovery and correction of the PIXY all-sites input issue, see [`99_PIXY_allsites_correction.md`](99_PIXY_allsites_correction.md).
