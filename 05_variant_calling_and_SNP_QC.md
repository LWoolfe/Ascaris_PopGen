# 05 — Variant calling and SNP QC

## Purpose

Document the **variant-calling and SNP quality-control workflow**, beginning with duplicate-marked BAMs produced by [`03_mapping_and_mapping_QC.md`](03_mapping_and_mapping_QC.md). 

## Inputs

| Input | Description |
|---|---|
| `03_Mapping/3a_markdup/SAMPLE.markdup.bam` | Indexed, duplicate-marked BAMs; prepared in stage 03 |
| `01_Reference/Ascaris_primary_scaffolds_plus_X.fasta` | Reference assembly prepared in stage 02 |
| `01_Reference/mappability_unique.merged.bed` | Unique-mappability mask from stage 02 |
| `04d_sample_QC/qcset.autosomes_146.txt` | Retained-sample list used for the **final**, not exploratory, VCF |

## Workflow

```text
Marked-duplicate BAMs
  → GATK HaplotypeCaller per-sample GVCFs
  → CombineGVCFs → GenotypeGVCFs
  → separate SNP / INDEL / MIXED records
  → site-level hard filters and PASS selection
  → unique-mappability restriction
  → biallelic SNP set → autosomal subset
  → autosomal missingness, depth and F statistics
  → sample QC review (7 additional hard failures; 11 flagged retained)
  → final 146-sample keep list
  → 90% per-site genotype call-rate filter
  → 13,484,735 biallelic autosomal SNPs
```

**Important chronology:** Three further samples were excluded based on prior mapping QC (documented in stage 03)

## Annotated scripts

### Part A — GVCF generation, joint genotyping and variant-site filtering

The blocks below are recovered command excerpts from the earlier variant-calling scripts; full QC and final VCF-construction jobs appear in Part B. They should not be pasted into one SLURM job.

### 1. Per-sample GVCF calling with GATK HaplotypeCaller

Each duplicate-marked BAM was called independently with GATK HaplotypeCaller in GVCF mode.

##### Input

```text
/parallel_scratch/lw01165/03_Mapping/3a_markdup/SAMPLE.markdup.bam
```

##### Reference

```text
/parallel_scratch/lw01165/01_Reference/Ascaris_primary_scaffolds_plus_X.fasta
```

##### Core command

```bash
gatk HaplotypeCaller \
    --emit-ref-confidence GVCF \
    --native-pair-hmm-threads 8 \
    -I "$IN_BAM" \
    -R "$REF_FA" \
    -O "$OUT_GVCF"
```

Some rerun scripts used the allocated CPU count for:

```bash
--native-pair-hmm-threads "${SLURM_CPUS_PER_TASK:-16}"
```

##### Output

```text
SAMPLE.g.vcf.gz
```

and corresponding index.

##### Note

Multiple HaplotypeCaller submission scripts are preserved because some samples were rerun after failed/incomplete jobs.

The underlying calling parameters are consistent in the recovered scripts:

```text
GVCF mode
same reference
same duplicate-marked BAM inputs
```

---

### 2. Combine per-sample GVCFs

The individual GVCFs were supplied to GATK `CombineGVCFs`.

An arguments file containing one input per line was generated:

```bash
ls *.g.vcf.gz |
awk '{print "-V "$0}' \
    > variantcall_sample_list.txt
```

The cohort GVCF was then created:

```bash
gatk CombineGVCFs \
    -R "$REF" \
    --arguments_file variantcall_sample_list.txt \
    -O merged_all_samples.g.vcf.gz \
    --tmp-dir "$TMPDIR"
```

The file was indexed:

```bash
gatk IndexFeatureFile \
    -I merged_all_samples.g.vcf.gz
```

##### Output

```text
merged_all_samples.g.vcf.gz
```

---

### 3. Joint genotyping

The combined GVCF was jointly genotyped:

```bash
gatk GenotypeGVCFs \
    -R "$REF" \
    -V merged_all_samples.g.vcf.gz \
    -O merged_all_samples.vcf.gz \
    --tmp-dir "$TMPDIR"
```

The cohort VCF was indexed:

```bash
gatk IndexFeatureFile \
    -I merged_all_samples.vcf.gz
```

##### Output

```text
merged_all_samples.vcf.gz
```

This was the main jointly genotyped cohort callset used for subsequent variant-type filtering.

---

### 4. Split joint callset by variant class

GATK `SelectVariants` was used to create separate datasets for:

```text
SNP
INDEL
MIXED
```

##### Commands

```bash
gatk SelectVariants \
    -R "$REF" \
    -V merged_all_samples.vcf.gz \
    --select-type-to-include SNP \
    -O merged_all_samples.SNPs.vcf.gz
```

```bash
gatk SelectVariants \
    -R "$REF" \
    -V merged_all_samples.vcf.gz \
    --select-type-to-include INDEL \
    -O merged_all_samples.INDELs.vcf.gz
```

```bash
gatk SelectVariants \
    -R "$REF" \
    -V merged_all_samples.vcf.gz \
    --select-type-to-include MIXED \
    -O merged_all_samples.MIXED.vcf.gz
```

All outputs were indexed.

---

### 5. Cohort variant-quality inspection

Before/alongside hard filtering, GATK annotations were exported for inspection.

Metrics included:

```text
QUAL
QD
FS
SOR
MQ
MQRankSum
ReadPosRankSum
```

##### Example

```bash
gatk VariantsToTable \
    -R "$REF" \
    -V merged_all_samples.vcf.gz \
    -O cohort.variant_qc.tsv \
    -F CHROM \
    -F POS \
    -F TYPE \
    -F QUAL \
    -F FILTER \
    -F QD \
    -F FS \
    -F SOR \
    -F MQ \
    -F MQRankSum \
    -F ReadPosRankSum \
    --show-filtered true
```

Percentile distributions were then calculated separately for SNPs and INDELs.

##### Annotation

The recovered supplementary material states that hard-filtering thresholds were selected with reference to the empirical distributions of these quality metrics.

---

### 6. SNP hard filtering

SNPs were tagged using GATK `VariantFiltration`.

##### Thresholds

```text
QD < 2.0   -> QD2
FS > 30.0  -> FS30
SOR > 3.0  -> SOR3
MQ < 40.0  -> MQ40
```

##### Command

```bash
gatk VariantFiltration \
    -R "$REF" \
    -V merged_all_samples.SNPs.vcf.gz \
    -O merged_all_samples.SNPs.tagged.vcf.gz \
    --filter-name "QD2" \
    --filter-expression "QD < 2.0" \
    --filter-name "FS30" \
    --filter-expression "FS > 30.0" \
    --filter-name "SOR3" \
    --filter-expression "SOR > 3.0" \
    --filter-name "MQ40" \
    --filter-expression "MQ < 40.0"
```

Only PASS SNPs were then retained:

```bash
gatk SelectVariants \
    -R "$REF" \
    -V merged_all_samples.SNPs.tagged.vcf.gz \
    --exclude-filtered \
    -O merged_all_samples.SNPs.filtered.vcf.gz
```


### 7. INDEL and MIXED hard filtering

INDELs used:

```text
QD < 2.0
FS > 50.0
SOR > 5.0
```

MIXED records were treated using the same threshold set as INDELs.

##### Example

```bash
gatk VariantFiltration \
    -R "$REF" \
    -V merged_all_samples.INDELs.vcf.gz \
    -O merged_all_samples.INDELs.tagged.vcf.gz \
    --filter-name "QD2" \
    --filter-expression "QD < 2.0" \
    --filter-name "FS50" \
    --filter-expression "FS > 50.0" \
    --filter-name "SOR5" \
    --filter-expression "SOR > 5.0"
```

PASS-only files were then produced with `SelectVariants --exclude-filtered`.

##### Note

The principal population-genomic analyses subsequently used the SNP dataset rather than the INDEL/MIXED datasets.

---

### 8. Post-filter QC

The tagged and PASS-only files were checked with `bcftools`.

Checks included:

- counts before and after filtering;
- FILTER-tag frequencies;
- confirmation that filtered outputs contained only PASS records;
- annotation presence;
- `bcftools stats`;
- Ts/Tv for relevant variant classes;
- inspection of example rejected variants.

This provides an explicit QC step between hard filtering and population-genomic processing.

---

### 9. Restriction to uniquely mappable regions

**Evidence: recovered historical code / documented outputs**

The PASS SNP, INDEL and MIXED datasets were intersected with:

```text
mappability_unique.merged.bed
```

For SNPs:

```bash
gatk SelectVariants \
    -R "$REF" \
    -V merged_all_samples.SNPs.tagged.vcf.gz \
    --exclude-filtered \
    -L "$BED" \
    -O merged_all_samples.SNPs.filtered.uniquemap.vcf.gz
```

Equivalent operations were performed for INDEL and MIXED records.

##### Main SNP output

```text
merged_all_samples.SNPs.filtered.uniquemap.vcf.gz
```

##### Annotation

This step is important: the principal SNP dataset **was** restricted to uniquely mappable regions.

The later all-sites discrepancy relates to a separate branch and is deliberately not folded into this document.

---
### Part B — Original sample-level SNP QC jobs

#### 10. `sampleQC`: retain biallelic SNPs and compute initial metrics

The SNP input has already passed the site-level filters and uniquely mappable-region restriction above. The first `.imiss`, `.idepth`, `.het` metrics span the chromosomes present at this point; the autosome-only metrics are recalculated in the following job.

```bash
#!/bin/bash
#SBATCH --job-name=sampleQC
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --output=sampleQC.%j.out
#SBATCH --error=sampleQC.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env

WORKDIR=/parallel_scratch/lw01165/04_Variants
INDIR=$WORKDIR/04c_uniquly_mappable_var
OUTDIR=$WORKDIR/04d_sample_QC

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# 1) Create biallelic SNP QC dataset
vcftools --gzvcf "$INDIR/merged_all_samples.SNPs.filtered.uniquemap.vcf.gz" \
  --min-alleles 2 --max-alleles 2 \
  --recode --recode-INFO-all \
  --out SNPs.QCset

# 2) Per-individual missingness
vcftools --vcf SNPs.QCset.recode.vcf \
  --missing-indv \
  --out SNPs.QCset

# 3) Per-individual mean depth (low-cov sanity check)
vcftools --vcf SNPs.QCset.recode.vcf \
  --depth \
  --out SNPs.QCset

# 4) Heterozygosity / inbreeding coefficient
vcftools --vcf SNPs.QCset.recode.vcf \
  --het \
  --out SNPs.QCset
```

#### 11. `split_auto_sex`: restrict to the 19 autosomes and recompute QC

Five X chromosomes were retained as a separate output; only autosomes proceeded to the main population-genomic analysis.

```bash
#!/bin/bash
#SBATCH --job-name=split_auto_sex
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --output=split_auto_sex.%j.out
#SBATCH --error=split_auto_sex.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env

WORKDIR=/parallel_scratch/lw01165/04_Variants/04d_sample_QC
cd $WORKDIR

INPUT_VCF=SNPs.QCset.recode.vcf
INPUT_GZ=SNPs.QCset.recode.vcf.gz

echo "Checking compression..."

# If gzipped file doesn't exist, bgzip it
if [ ! -f "$INPUT_GZ" ]; then
    echo "Compressing with bgzip..."
    bgzip -@ 8 $INPUT_VCF
fi

# Index if index does not exist
if [ ! -f "$INPUT_GZ.tbi" ]; then
    echo "Indexing VCF..."
    bcftools index -t $INPUT_GZ
fi

echo "Splitting autosomes..."

bcftools view -Oz \
  -r CM024168.1,CM024169.1,CM024170.1,CM024171.1,CM024172.1,CM024173.1,CM024174.1,CM024175.1,CM024176.1,CM024177.1,CM024178.1,CM024179.1,CM024180.1,CM024181.1,CM024182.1,CM024183.1,CM024184.1,CM024185.1,CM024186.1 \
  $INPUT_GZ > SNPs.QCset.autosomes.vcf.gz

bcftools index -t SNPs.QCset.autosomes.vcf.gz


echo "Splitting sex chromosomes..."

bcftools view -Oz \
  -r CM024187.1,CM024188.1,CM024189.1,CM024190.1,CM024191.1 \
  $INPUT_GZ > SNPs.QCset.sexchr.vcf.gz

bcftools index -t SNPs.QCset.sexchr.vcf.gz


echo "Recomputing QC metrics on autosomes..."

vcftools --gzvcf SNPs.QCset.autosomes.vcf.gz --missing-indv --out SNPs.autosomes
vcftools --gzvcf SNPs.QCset.autosomes.vcf.gz --depth       --out SNPs.autosomes
vcftools --gzvcf SNPs.QCset.autosomes.vcf.gz --het         --out SNPs.autosomes

echo "All done."
```

#### 12. Define sample-level hard failures and borderline flags

Thresholds are applied to the **autosomal** `.het`, `.idepth` and `.imiss` files. This is sample-level filtering

```bash
# remove: F < -0.4
awk 'NR>1 && $6 < -0.4 {print $1}' SNPs.autosomes.het > remove.Flt-0.4.txt

# remove: depth < 5
awk 'NR>1 && $3 < 5 {print $1}' SNPs.autosomes.idepth > remove.depthlt5.txt

# remove: missingness > 0.05
awk 'NR>1 && $5 > 0.05 {print $1}' SNPs.autosomes.imiss > remove.missgt0.05.txt

cat remove.Flt-0.4.txt remove.depthlt5.txt remove.missgt0.05.txt \
  | sort -u > remove.samples.txt

# flag: -0.4 to -0.2
awk 'NR>1 && $6 >= -0.4 && $6 < -0.2 {print $1}' SNPs.autosomes.het > flag.F-0.4to-0.2.txt

# flag: depth 5-8 (optional)
awk 'NR>1 && $3 >= 5 && $3 < 8 {print $1}' SNPs.autosomes.idepth > flag.depth5to8.txt

cat flag.F-0.4to-0.2.txt flag.depth5to8.txt \
  | sort -u | comm -23 - remove.samples.txt > flag.samples.txt

cut -f1 SNPs.autosomes.het | tail -n +2 | sort -u > all.samples.txt
comm -23 all.samples.txt <(cat remove.samples.txt flag.samples.txt | sort -u) > keep.samples.txt
```


The recovered source creates two different lists with different intended purposes:

- `remove.samples.txt`: hard-fail samples only (criteria F < −0.4, mean depth < 5× or missingness > 5%).
- `flag.samples.txt`: moderate excess heterozygosity or depth 5–8×, with hard failures excluded.
- `keep.samples.txt`: an **earlier strict-keep list that removes both** the hard failures **and the flagged samples**. This list generated an exploratory 135-sample branch and was **not** the final thesis dataset.

The final dataset **retained the 11 flagged samples**. The three mapping-QC exclusions plus seven additional variant/sample-QC exclusions account for 156 → 146. The original command that created `qcset.autosomes_146.txt` has not been recovered in this archive; the final downstream script references this existing list directly. Preserve that list as a provenance file if it is still available.

#### 13. Samples documented as excluded from the final analyses

| QC stage | Samples | Reason |
|---|---|---|
| Mapping QC (stage 03) | `AL12SLB_2`, `SRR31757526`, `SRR31757525` | Median depth < 5× **and** breadth < 0.85 |
| SNP QC | `DE_H1a`, `SRR31757540`, `SRR31757546` | Low mean depth and elevated missingness |
| SNP QC | `DE_WS11A`, `DE_WS13A`, `SRR31692059`, `UK16F13605_b` | Excess heterozygosity (F < −0.4) |

These exclusions are documented by the QC records and supplementary methods. The `qcset.autosomes_146.txt` file is the **authoritative input** to the final VCF-building command; the code does not recompute that list itself.

### Part C — Historical strict-keep branch (NOT final thesis VCF)

#### 14. `mk_cleanSNPvcf`: historical `keep.samples.txt` branch

This outputs `Ascaris.autosomes.SNPs.clean.vcf.gz` from `keep.samples.txt` (excluding hard-fail and borderline samples). It should **not** be confused with the 146-sample final VCF.

```bash
#!/bin/bash
#SBATCH --job-name=mk_cleanSNPvcf
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=06:00:00
#SBATCH --output=mk_cleanSNPvcf.%j.out
#SBATCH --error=mk_cleanSNPvcf.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env

QC_DIR=/parallel_scratch/lw01165/04_Variants/04d_sample_QC
ANALYSIS_DIR=/parallel_scratch/lw01165/05_analysis

IN_VCF=${QC_DIR}/SNPs.QCset.autosomes.vcf.gz
KEEP_LIST=${QC_DIR}/keep.samples.txt

OUT_VCF=${ANALYSIS_DIR}/Ascaris.autosomes.SNPs.clean.vcf.gz
OUT_PREFIX=${ANALYSIS_DIR}/Ascaris.autosomes.SNPs.clean

mkdir -p "${ANALYSIS_DIR}"

echo "Input VCF:  ${IN_VCF}"
echo "Keep list:  ${KEEP_LIST}"
echo "Output VCF: ${OUT_VCF}"

# Ensure input index exists
if [ ! -f "${IN_VCF}.tbi" ]; then
  echo "Index missing for ${IN_VCF} - creating..."
  bcftools index -t "${IN_VCF}"
fi

# Subset samples and keep biallelic SNPs
bcftools view \
  -S "${KEEP_LIST}" \
  -m2 -M2 -v snps \
  -Oz -o "${OUT_VCF}" \
  "${IN_VCF}"

# Index output
bcftools index -t "${OUT_VCF}"

# Simple stats for logs / methods
echo -e "metric\tvalue" > "${OUT_PREFIX}.stats.tsv"
echo -e "n_keep_samples\t$(wc -l < ${KEEP_LIST})" >> "${OUT_PREFIX}.stats.tsv"
echo -e "n_samples_in_vcf\t$(bcftools query -l ${OUT_VCF} | wc -l)" >> "${OUT_PREFIX}.stats.tsv"
echo -e "n_variants\t$(bcftools view -H ${OUT_VCF} | wc -l)" >> "${OUT_PREFIX}.stats.tsv"

echo "Done. Outputs:"
ls -lh "${OUT_VCF}" "${OUT_VCF}.tbi" "${OUT_PREFIX}.stats.tsv"
```

### Part D — Canonical final 146-sample VCF

#### 15. `mk_146_maxmiss90`: final biallelic, autosomal SNP VCF

It starts from the autosomal biallelic SNP QC VCF, uses `qcset.autosomes_146.txt` rather than the earlier `keep.samples.txt`, and retains flagged samples represented in that final keep list. It then filters SNP sites by ≥90% sample genotype call rate.

```bash
#!/bin/bash
#SBATCH --job-name=mk_146_maxmiss90
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=16
#SBATCH --mem=150G
#SBATCH --time=12:00:00
#SBATCH --output=mk_146_maxmiss90.%j.out
#SBATCH --error=mk_146_maxmiss90.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env

INDIR=/parallel_scratch/lw01165
OUTDIR=/parallel_scratch/lw01165/05_analysis/05b_146_data
mkdir -p "${OUTDIR}"

IN_VCF=${INDIR}/04_Variants/04d_sample_QC/SNPs.QCset.autosomes.vcf.gz
KEEP_LIST=${INDIR}/04_Variants/04d_sample_QC/qcset.autosomes_146.txt

MAX_MISSING=0.90

OUT_PREFIX=${OUTDIR}/Ascaris.autosomes.SNPs.146.maxmiss${MAX_MISSING}
OUT_VCF=${OUT_PREFIX}.vcf.gz

echo "Input VCF:  ${IN_VCF}"
echo "Keep list:  ${KEEP_LIST}"
echo "Output dir: ${OUTDIR}"
echo "Max-missing: ${MAX_MISSING}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-16}"

if [ ! -f "${IN_VCF}.tbi" ]; then
  echo "Index missing for ${IN_VCF} - creating..."
  bcftools index -t "${IN_VCF}"
fi

SUBSET_VCF=${OUT_PREFIX}.subset.biallelicSNPs.vcf.gz

bcftools view \
  -S "${KEEP_LIST}" \
  -m2 -M2 -v snps \
  -Oz -o "${SUBSET_VCF}" \
  "${IN_VCF}"

bcftools index -t "${SUBSET_VCF}"

vcftools --gzvcf "${SUBSET_VCF}" \
  --max-missing "${MAX_MISSING}" \
  --recode --recode-INFO-all \
  --out "${OUT_PREFIX}"

bgzip -@ "${SLURM_CPUS_PER_TASK:-16}" -f "${OUT_PREFIX}.recode.vcf"
mv -f "${OUT_PREFIX}.recode.vcf.gz" "${OUT_VCF}"
bcftools index -t "${OUT_VCF}"

echo -e "metric\tvalue" > "${OUT_PREFIX}.stats.tsv"
echo -e "n_keep_samples\t$(wc -l < ${KEEP_LIST})" >> "${OUT_PREFIX}.stats.tsv"
echo -e "n_samples_in_vcf\t$(bcftools query -l ${OUT_VCF} | wc -l)" >> "${OUT_PREFIX}.stats.tsv"
echo -e "n_variants_before_maxmiss\t$(bcftools view -H ${SUBSET_VCF} | wc -l)" >> "${OUT_PREFIX}.stats.tsv"
echo -e "n_variants_after_maxmiss\t$(bcftools view -H ${OUT_VCF} | wc -l)" >> "${OUT_PREFIX}.stats.tsv"

echo "Done. Outputs:"
ls -lh "${SUBSET_VCF}" "${SUBSET_VCF}.tbi" "${OUT_VCF}" "${OUT_VCF}.tbi" "${OUT_PREFIX}.stats.tsv"
```


### Part E — Outputs, quality controls and methodological omissions

#### Final data product

```text
/parallel_scratch/lw01165/05_analysis/05b_146_data/Ascaris.autosomes.SNPs.146.maxmiss0.90.vcf.gz
```

Final verification previously reported:

| Measure | Result |
|---|---:|
| Retained samples | **146** |
| Biallelic SNP records | **13,484,735** |
| Multiallelic records | 0 |
| Non-SNP records | 0 |

The 146-sample subset before the 90% site-call filter contained approximately **13,870,294** SNPs  



#### Base Quality Score Recalibration (BQSR) — NOT performed

BQSR is **not** present in the original processing history. Recalibration ordinarily uses a trusted set of known polymorphic sites to avoid mistaking real variation for sequencing error. A well-curated *Ascaris* known-sites resource was not available for this study, making an off-the-shelf BQSR workflow difficult to validate.

## Outputs and QC checks

- `SNPs.QCset.recode.vcf(.gz)`: biallelic, uniquely mappable SNP QC set prior to autosome split.
- `SNPs.QCset.autosomes.vcf.gz`: autosomal SNP QC set and source of final 146-sample subset.
- `SNPs.autosomes.{imiss,idepth,het}`: sample-level autosomal QC values.
- `keep.samples.txt` / `Ascaris.autosomes.SNPs.clean.vcf.gz`: **exploratory 135-sample route** only.
- `qcset.autosomes_146.txt`: original keep list used to create the final analysis cohort (should be archived).
- `Ascaris.autosomes.SNPs.146.maxmiss0.90.vcf.gz`: **canonical final 146-sample SNP VCF**.
