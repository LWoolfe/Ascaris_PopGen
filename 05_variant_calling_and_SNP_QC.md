# 05 — Variant calling and SNP QC

## Purpose

Document the original *Ascaris* **nuclear variant-calling, joint-genotyping, hard-filtering, sample-QC and final SNP-dataset construction** pathway used for the PhD thesis/paper. The principal output contains **146 samples** and **13,484,735 biallelic autosomal SNPs**. This is the historical SNP pipeline; it must not be silently replaced by the later all-sites/PIXY correction.

## Inputs

| Input | Source / provenance |
|---|---|
| `SAMPLE.markdup.bam` | Duplicate-marked, indexed Ascaris alignments from `03_mapping_and_mapping_QC.md` |
| `Ascaris_primary_scaffolds_plus_X.fasta` | Chromosome-level reference from `02_reference_preparation_and_mappability.md` |
| `mappability_unique.merged.bed` | GenMap uniquely mappable regions from reference preparation |
| `SNPs.QCset.autosomes.vcf.gz` | Autosomal, biallelic QC set produced from the filtered SNP branch |
| `qcset.autosomes_146.txt` | Final list retaining 146 samples (7 variant-QC exclusions; flagged samples retained) |

## Workflow

```text
MarkDuplicates BAMs → GATK HaplotypeCaller (GVCF per sample)
→ CombineGVCFs → GenotypeGVCFs
→ split SNP / INDEL / MIXED → hard-filter → retain PASS
→ restrict SNPs to uniquely mappable intervals
→ retain biallelic SNPs and 19 autosomes
→ individual missingness/depth/F QC
→ remove 7 variant-QC hard fails; retain 11 flagged samples
→ final 146-sample subset
→ sites called in ≥90% of samples
→ 13,484,735 biallelic SNPs
```

## Annotated scripts

The main body below preserves the previously compiled, thoroughly annotated walk-through of **the actual historical workflow**, including source commands and QC results. Some short commands are explanatory extracts from longer HPC submission scripts; comments distinguish direct code evidence from reconstructed/provenance notes. Earlier preprocessing/reference/mapping sections are cross-referenced for completeness, not additional processing steps.

### Historical pipeline overview

```text
156 paired-end WGS datasets
        |
        v
FASTQ acquisition / validation
        |
        v
fastp trimming and read filtering
        |
        v
Deacon host depletion
(human + pig reference)
        |
        v
BWA-MEM mapping
to chromosome-level Ascaris reference
(19 autosomes + 5 X chromosomes)
        |
        v
coordinate-sorted BAM + index
        |
        v
GATK MarkDuplicates
        |
        v
mapping-level QC
        |
        |-- 3 mapping-QC failures removed
        v
per-sample GATK HaplotypeCaller
GVCF mode
        |
        v
CombineGVCFs
        |
        v
GenotypeGVCFs
        |
        v
joint-called cohort VCF
        |
        v
split SNP / INDEL / MIXED
        |
        v
hard-filter variants
        |
        v
retain PASS variants
        |
        v
restrict to uniquely mappable regions
        |
        v
biallelic SNP QC dataset
        |
        v
restrict QC to 19 autosomes
        |
        v
sample-level SNP QC
missingness + mean depth + F
        |
        |-- 7 variant-QC failures removed
        v
146 retained samples
        |
        v
site missingness filter
--max-missing 0.90
        |
        v
FINAL ANALYSIS VCF

Ascaris.autosomes.SNPs.146.maxmiss0.90.vcf.gz

13,484,735 biallelic SNPs
146 samples
```

---

### 1. Raw sequencing datasets

**Status: CONFIRMED**

The study began with:

```text
28 newly generated paired-end WGS datasets
128 publicly available SRA WGS datasets
-----------------------------------------
156 samples total
```

The datasets remained separate throughout preprocessing, mapping and variant calling.

The SRA datasets were retrieved using the SRA Toolkit using:

```bash
prefetch
fasterq-dump --split-files
pigz
```

A stronger recovery script was subsequently used for incomplete/missing SRA accessions. This additionally checked:

- existence of both paired FASTQ files;
- gzip integrity;
- equal R1 and R2 read counts.

##### Main logic

```bash
prefetch "$ACC"

fasterq-dump \
    --split-files \
    "$SRA_PATH"

pigz "${ACC}_1.fastq"
pigz "${ACC}_2.fastq"
```

##### Annotation

This stage generated the paired FASTQ files used for preprocessing. It did not combine different samples.

---

### 2. Read trimming and quality filtering with fastp

**Status: CONFIRMED**

Paired reads were processed using `fastp`.

Recovered parameters were:

```text
paired-end adapter detection    enabled
trim from 5' end of R1          20 bp
trim from 5' end of R2          20 bp
minimum retained read length    50 bp
```

##### Core command

```bash
fastp \
    -i "$r1" \
    -I "$r2" \
    -o "${sample}_trimmed_R1.fastq" \
    -O "${sample}_trimmed_R2.fastq" \
    --thread "$THREADS" \
    --html "${sample}_fastp.html" \
    --detect_adapter_for_pe \
    --trim_front1 20 \
    --trim_front2 20 \
    --length_required 50
```

Outputs were compressed and retained as paired files.

##### Main outputs

```text
SAMPLE_trimmed_R1.fastq.gz
SAMPLE_trimmed_R2.fastq.gz
SAMPLE_fastp.html
```

##### QC

FastQC and MultiQC were used to assess processed read quality.

---

### 3. Host-read depletion

**Status: CONFIRMED**

Potential host-derived reads were removed using Deacon.

A combined human + pig reference was constructed from:

```text
Homo sapiens    GCF_000001405.40
Sus scrofa      GCF_000003025.6
```

The two genomic references were concatenated:

```bash
cat "$HUMAN_FASTA" "$PIG_FASTA" \
    > human_pig_genome.fna
```

A Deacon index was then generated:

```bash
deacon index build \
    human_pig_genome.fna \
    > human_pig_genome.idx
```

The fastp-processed reads were filtered against this combined host reference.

##### Main output concept

```text
host-depleted paired-end reads
```

These host-depleted reads were then used for nuclear mapping.

##### Documentation note

The exact local source filenames for the original human and pig FASTA files were not retained, although their NCBI accessions and the analytical operation are documented.

This is a **minor provenance gap**, not a gap in the biological identity of the reference genomes.

---

### 4. Ascaris reference-genome preparation

**Status: CONFIRMED**

Reference assembly:

```text
Ascaris suum ASM1343314v1
GenBank/NCBI accession GCA_013433145.1
```

The complete assembly was restricted to the 24 chromosome-level molecules:

```text
19 autosomes: CM024168.1–CM024186.1
5 X chromosomes: CM024187.1–CM024191.1
```

Unplaced/unlocalised scaffolds were excluded.

Final reference:

```text
Ascaris_primary_scaffolds_plus_X.fasta
```

The reference was prepared for mapping and GATK:

```bash
bwa index Ascaris_primary_scaffolds_plus_X.fasta

samtools faidx Ascaris_primary_scaffolds_plus_X.fasta

gatk CreateSequenceDictionary \
    --REFERENCE Ascaris_primary_scaffolds_plus_X.fasta
```

##### Main index files

```text
.fasta
.fasta.amb
.fasta.ann
.fasta.bwt
.fasta.fai
.fasta.pac
.fasta.sa
.dict
```

---

### 5. Generation of uniquely mappable regions

**Status: CONFIRMED**

GenMap was used to identify regions of the chromosome-level reference that were uniquely mappable.

Parameters:

```text
k-mer length        35 bp
permitted mismatch  1
retained score      1
```

##### Core commands

```bash
genmap index \
    -F "$REF" \
    -I "$INDEX_DIR"

genmap map \
    -K 35 \
    -E 1 \
    -I "$INDEX_DIR" \
    -O "$PREFIX" \
    -bg
```

Intervals with score `1` were extracted and ultimately converted to the merged BED mask:

```text
/parallel_scratch/lw01165/01_Reference/mappability_unique.merged.bed
```

##### Purpose

This mask was later used to restrict the filtered variant datasets to genomic regions in which reads could be mapped uniquely.

---

### 6. Read mapping with BWA-MEM

**Status: CONFIRMED**

Host-depleted paired reads were mapped against:

```text
Ascaris_primary_scaffolds_plus_X.fasta
```

using BWA-MEM.

Read-group information was added during mapping.

The alignments were then coordinate sorted and indexed with samtools.

##### Core workflow

```bash
bwa mem \
    -t "$THREADS" \
    -R "$READ_GROUP" \
    "$REF" \
    "$R1" \
    "$R2" |
samtools sort \
    -@ "$THREADS" \
    -o "${SAMPLE}.bam"

samtools index \
    -@ "$THREADS" \
    "${SAMPLE}.bam"

samtools flagstat \
    "${SAMPLE}.bam"
```

##### Main output

```text
/parallel_scratch/lw01165/03_Mapping/SAMPLE.bam
```

##### Documentation gap

The recovered full BWA submission script is specifically documented for the 128 public SRA datasets.

The newly generated samples were documented as having been mapped using the same BWA-MEM/read-group/sorting/indexing commands in earlier batches, but the exact original submission script for every one of those earlier batches was not recovered.

**Interpretation:** the analytical procedure is documented, but batch-level provenance is incomplete for the earliest mapping jobs.

---

### 7. Duplicate marking

**Status: CONFIRMED**

GATK `MarkDuplicates` was applied to each mapped BAM.

##### Core command

```bash
gatk MarkDuplicates \
    --INPUT "$IN_BAM" \
    --OUTPUT "$OUT_BAM" \
    --METRICS_FILE "$METRICS" \
    --VALIDATION_STRINGENCY LENIENT

samtools index \
    -@ "$THREADS" \
    "$OUT_BAM"
```

##### Main outputs

```text
SAMPLE.markdup.bam
SAMPLE.markdup.bam.bai
SAMPLE.metrics.txt
```

##### Important annotation

Duplicates were **marked**, not physically deleted from the BAM.

The duplicate flags could then be recognised by downstream GATK processing.

---

### 8. Mapping-level QC

**Status: CONFIRMED**

Mapping QC was calculated from the duplicate-marked BAMs.

Metrics included:

```text
median mapping quality
breadth of coverage
median depth
```

Additional 25-kb window depth and breadth summaries were also generated for genome-position/chromosome-level assessment.

##### Formal failure thresholds

```text
median depth < 5×
AND
breadth < 0.85
```

A sample was formally removed at mapping QC only when **both** criteria failed.

##### Mapping-QC exclusions

Three samples are documented as mapping-QC failures:

```text
AL12SLB_2
SRR31757526
SRR31757525
```

Reason:

```text
low breadth (<0.85)
AND
low median depth (<5×)
```

Thus:

```text
156 initial samples
- 3 mapping-QC failures
= 153 samples entering/remaining for subsequent cohort-level assessment
```

##### Note

The code also created a `Watch` category for low-ranking retained samples. `Watch` status itself was not intended to determine exclusion.

---

### 9. Per-sample GVCF calling with GATK HaplotypeCaller

**Status: CONFIRMED**

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

### 10. Combine per-sample GVCFs

**Status: CONFIRMED**

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

### 11. Joint genotyping

**Status: CONFIRMED**

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

### 12. Split joint callset by variant class

**Status: CONFIRMED**

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

### 13. Cohort variant-quality inspection

**Status: CONFIRMED**

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

### 14. SNP hard filtering

**Status: CONFIRMED**

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

##### Important annotation

`MQRankSum` and `ReadPosRankSum` were inspected and reported in QC tables, but no recovered command applies them as SNP hard-filter thresholds.

Do **not** state that they were used as hard filters unless another original script demonstrates this.

---

### 15. INDEL and MIXED hard filtering

**Status: CONFIRMED**

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

### 16. Post-filter QC

**Status: CONFIRMED**

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

### 17. Restriction to uniquely mappable regions

**Status: CONFIRMED**

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

### 18. Construction of biallelic SNP QC dataset

**Status: CONFIRMED**

A biallelic SNP dataset was generated for sample-level QC:

```bash
vcftools \
    --gzvcf merged_all_samples.SNPs.filtered.uniquemap.vcf.gz \
    --min-alleles 2 \
    --max-alleles 2 \
    --recode \
    --recode-INFO-all \
    --out SNPs.QCset
```

The source file was already a SNP-only dataset, so the allele-count restriction removed multiallelic SNP records and retained biallelic SNP sites.

---

### 19. Initial sample-level SNP QC

**Status: CONFIRMED**

Three per-individual statistics were generated:

```bash
vcftools \
    --vcf SNPs.QCset.recode.vcf \
    --missing-indv \
    --out SNPs.QCset
```

```bash
vcftools \
    --vcf SNPs.QCset.recode.vcf \
    --depth \
    --out SNPs.QCset
```

```bash
vcftools \
    --vcf SNPs.QCset.recode.vcf \
    --het \
    --out SNPs.QCset
```

These measured:

```text
per-individual genotype missingness
mean depth
heterozygosity / inbreeding coefficient F
```

---

### 20. Restriction of SNP QC to the 19 autosomes

**Status: CONFIRMED**

The biallelic SNP QC VCF was divided into autosomal and X-chromosome datasets.

Autosomes retained:

```text
CM024168.1–CM024186.1
```

X chromosomes:

```text
CM024187.1–CM024191.1
```

##### Autosome command

```bash
bcftools view \
    -Oz \
    -r CM024168.1,CM024169.1,CM024170.1,CM024171.1,CM024172.1,CM024173.1,CM024174.1,CM024175.1,CM024176.1,CM024177.1,CM024178.1,CM024179.1,CM024180.1,CM024181.1,CM024182.1,CM024183.1,CM024184.1,CM024185.1,CM024186.1 \
    SNPs.QCset.recode.vcf.gz \
    > SNPs.QCset.autosomes.vcf.gz
```

The autosomal QC metrics were then recalculated:

```bash
vcftools \
    --gzvcf SNPs.QCset.autosomes.vcf.gz \
    --missing-indv \
    --out SNPs.autosomes

vcftools \
    --gzvcf SNPs.QCset.autosomes.vcf.gz \
    --depth \
    --out SNPs.autosomes

vcftools \
    --gzvcf SNPs.QCset.autosomes.vcf.gz \
    --het \
    --out SNPs.autosomes
```

##### Annotation

The final population-genomic SNP dataset therefore focused on the 19 autosomes rather than the five X chromosomes.

---

### 21. Variant-calling sample QC thresholds

**Status: CONFIRMED**

Hard-fail lists were generated using:

##### Excess heterozygosity

```bash
awk 'NR>1 && $6 < -0.4 {print $1}' \
    SNPs.autosomes.het \
    > remove.Flt-0.4.txt
```

Criterion:

```text
F < -0.4
```

##### Low mean depth

```bash
awk 'NR>1 && $3 < 5 {print $1}' \
    SNPs.autosomes.idepth \
    > remove.depthlt5.txt
```

Criterion:

```text
mean depth < 5×
```

##### High missingness

```bash
awk 'NR>1 && $5 > 0.05 {print $1}' \
    SNPs.autosomes.imiss \
    > remove.missgt0.05.txt
```

Criterion:

```text
>5% missing genotypes
```

Hard-fail IDs were combined:

```bash
cat \
    remove.Flt-0.4.txt \
    remove.depthlt5.txt \
    remove.missgt0.05.txt |
sort -u \
    > remove.samples.txt
```

---

### 22. Flagged / borderline samples

**Status: RESOLVED FROM EARLIER PIPELINE CHAT**

Additional samples were flagged for review.

##### Moderate excess heterozygosity

```bash
awk 'NR>1 && $6 >= -0.4 && $6 < -0.2 {print $1}' \
    SNPs.autosomes.het \
    > flag.F-0.4to-0.2.txt
```

##### Lower retained depth

```bash
awk 'NR>1 && $3 >= 5 && $3 < 8 {print $1}' \
    SNPs.autosomes.idepth \
    > flag.depth5to8.txt
```

Flagged IDs were combined:

```bash
cat \
    flag.F-0.4to-0.2.txt \
    flag.depth5to8.txt |
sort -u |
comm -23 - remove.samples.txt \
    > flag.samples.txt
```

An earlier exploratory command produced:

```bash
comm -23 all.samples.txt \
    <(cat remove.samples.txt flag.samples.txt | sort -u) \
    > keep.samples.txt
```

This generated a **135-sample strict-keep list** by excluding both the seven hard-fail samples and the eleven flagged samples.

This was **not the final thesis sample set**.

The final 146-sample rebuild deliberately retained the eleven flagged samples and removed only the seven hard-fail samples from the 153-sample autosomal QC set.

The final keep list was:

```text
/parallel_scratch/lw01165/04_Variants/04d_sample_QC/qcset.autosomes_146.txt
```

The intended construction was:

```text
153 autosomal QC samples
- 7 hard-fail variant-QC samples
= 146 retained samples
```

The eleven flagged samples remained in the 146-sample dataset.

---

### 23. Documented QC exclusions

**Status: CONFIRMED**

Three samples had previously failed mapping-level QC:

```text
AL12SLB_2
SRR31757526
SRR31757525
```

Seven additional samples failed variant/sample-level QC.

##### Low depth + high missingness

```text
DE_H1a
SRR31757540
SRR31757546
```

##### Excess heterozygosity

```text
DE_WS11A
DE_WS13A
SRR31692059
UK16F13605_b
```

Therefore the final accounting was:

```text
156 initial samples
- 3 mapping-QC failures
= 153 autosomal QC samples

153 autosomal QC samples
- 7 variant/sample-QC failures
= 146 final samples
```

The eleven intermediate flagged samples were retained.

---

### 24. Final 146-sample autosomal SNP rebuild

**Status: RESOLVED FROM EARLIER PIPELINE CHAT**

The final thesis VCF was rebuilt directly from the autosomal biallelic SNP QC dataset rather than from the earlier 135-sample `clean.vcf.gz` branch.

##### Input

```text
/parallel_scratch/lw01165/04_Variants/04d_sample_QC/SNPs.QCset.autosomes.vcf.gz
```

##### Final 146-sample keep list

```text
/parallel_scratch/lw01165/04_Variants/04d_sample_QC/qcset.autosomes_146.txt
```

##### Output directory

```text
/parallel_scratch/lw01165/05_analysis/05b_146_data
```

##### Recovered final rebuild variables

```bash
IN_VCF=/parallel_scratch/lw01165/04_Variants/04d_sample_QC/SNPs.QCset.autosomes.vcf.gz

KEEP_LIST=/parallel_scratch/lw01165/04_Variants/04d_sample_QC/qcset.autosomes_146.txt

OUTDIR=/parallel_scratch/lw01165/05_analysis/05b_146_data

OUT_PREFIX=${OUTDIR}/Ascaris.autosomes.SNPs.146.maxmiss0.90

OUT_VCF=${OUT_PREFIX}.vcf.gz

SUBSET_VCF=${OUT_PREFIX}.subset.biallelicSNPs.vcf.gz
```

##### Subset to the final 146 samples and confirm biallelic SNPs

```bash
bcftools view \
    -S "${KEEP_LIST}" \
    -m2 \
    -M2 \
    -v snps \
    -Oz \
    -o "${SUBSET_VCF}" \
    "${IN_VCF}"
```

This step retained the final 146 sample identities while explicitly restricting the file to biallelic SNPs.

---

### 25. Final site-level missingness filter

**Status: CONFIRMED**

The 146-sample subset was filtered so that each SNP had genotype calls in at least 90% of samples.

```bash
vcftools \
    --gzvcf "${SUBSET_VCF}" \
    --max-missing 0.90 \
    --recode \
    --recode-INFO-all \
    --out "${OUT_PREFIX}"
```

The output VCF was then compressed.

```bash
bgzip -f "${OUT_PREFIX}.recode.vcf"
```

The final compressed file was moved to the canonical thesis filename:

```bash
mv -f \
    "${OUT_PREFIX}.recode.vcf.gz" \
    "${OUT_VCF}"
```

and indexed:

```bash
bcftools index -t "${OUT_VCF}"
```

Thus the filename transition to the final analysis VCF is now documented.

---

### 26. Earlier `Ascaris.autosomes.SNPs.clean.vcf.gz` branch

**Status: CLARIFIED**

Earlier exploratory/sample-QC code generated:

```text
Ascaris.autosomes.SNPs.clean.vcf.gz
```

from a stricter `keep.samples.txt` list.

That list excluded:

```text
7 hard-fail samples
+
11 flagged samples
```

and therefore represented a **135-sample strict-keep branch**.

It was useful during QC/missingness exploration but was **not the final thesis dataset**.

The final 146-sample analyses used the rebuilt VCF described in Sections 24–25.

---

### 27. Final analysis SNP VCF

**Status: CONFIRMED**

Final file:

```text
/parallel_scratch/lw01165/05_analysis/05b_146_data/Ascaris.autosomes.SNPs.146.maxmiss0.90.vcf.gz
```

Independent audit showed:

```text
samples:                 146
total records:           13,484,735
biallelic SNP records:   13,484,735
multiallelic records:    0
non-SNP records:         0
```

The pre-90%-call-rate 146-sample subset contained approximately:

```text
13,870,294 SNPs
```

and final site filtering retained:

```text
13,484,735 SNPs
```

---

### 28. Transition into downstream analyses

The following analyses consume the final 146-sample VCF or datasets derived from it:

```text
balanced 65-sample sensitivity dataset
PLINK conversion and LD pruning
PCA
ADMIXTURE
IBS distances
PERMANOVA and dispersion testing
NJ trees
ML SNP tree
AMOVA
Wahlund-effect checks
PIXY diversity/differentiation analyses
chromosome-level analyses
candidate high-differentiation analyses
```

The balanced 65-sample dataset was generated downstream of the primary 146-sample SNP VCF.

---

### 29. Remaining provenance/documentation issues

The major final-VCF gaps identified in the first version of this document are now resolved.

#### RESOLVED — final sample-selection logic

The final dataset retained the eleven flagged samples and removed only the seven variant-QC hard failures from the 153-sample autosomal QC set.

#### RESOLVED — final VCF construction command

The 146-sample dataset was rebuilt from:

```text
SNPs.QCset.autosomes.vcf.gz
```

using:

```text
qcset.autosomes_146.txt
```

and an explicit biallelic-SNP restriction.

#### RESOLVED — final filename transition

The 90%-call-rate VCF was compressed, moved to:

```text
Ascaris.autosomes.SNPs.146.maxmiss0.90.vcf.gz
```

and indexed.

#### REMAINING MINOR DOCUMENTATION ISSUE — software versions

Not every historical script prints the software version used. Versions should therefore be taken from retained environments/logs or the thesis methods where already recorded rather than inferred.

#### REMAINING MINOR DOCUMENTATION ISSUE — early mapping batch provenance

The complete BWA-MEM script is retained for the public SRA batch. The newly generated samples were processed using the same analytical commands in earlier batches, but not every original batch submission file has been recovered.

This is a provenance/documentation limitation rather than evidence of a different analytical method.

---

### 30. Methodological steps not present in the recovered pipeline

These are **not automatically errors**, but should not be claimed in the thesis unless supported elsewhere.

#### No recovered per-genotype DP/GQ masking step

No recovered post-genotyping script was found that explicitly masks individual genotypes based on thresholds such as:

```text
DP
GQ
```

The recovered pipeline instead applies:

- cohort/site hard filters;
- mappability restriction;
- sample-level depth/missingness/heterozygosity QC;
- final site-level call rate ≥90%.

Therefore, do not state that genotype calls were individually masked using DP/GQ thresholds unless another script demonstrates this.

---

#### No recovered BQSR step

No Base Quality Score Recalibration step is present in the recovered pipeline.

This should simply be recorded as absent from the workflow unless there is separate evidence that it was performed.

It should not be retroactively added to the pipeline description.

---

### 31. Final canonical SNP pipeline for the thesis/code archive

The most defensible concise representation is:

```text
Raw paired reads
  ↓
fastp trimming / filtering
  ↓
human + pig host depletion with Deacon
  ↓
BWA-MEM mapping to chromosome-level Ascaris reference
  ↓
coordinate sort + BAM index
  ↓
GATK MarkDuplicates
  ↓
mapping QC
  ↓
GATK HaplotypeCaller in GVCF mode
  ↓
CombineGVCFs
  ↓
GenotypeGVCFs
  ↓
split SNP / INDEL / MIXED
  ↓
hard-filter variants
  ↓
retain PASS sites
  ↓
restrict to uniquely mappable regions
  ↓
retain biallelic SNPs
  ↓
restrict QC / downstream population dataset to 19 autosomes
  ↓
sample QC:
    depth
    missingness
    heterozygosity / F
  ↓
remove 3 mapping-QC + 7 variant-QC failures
  ↓
146 samples
  ↓
retain SNP sites called in ≥90% samples
  ↓
13,484,735 autosomal biallelic SNPs
```

---

### 32. Important separation from the all-sites / PIXY amendment

The SNP branch above **already contained the uniquely mappable-region restriction**.

The later issue identified during final audit concerned the separately prepared all-sites/invariant-site input used for PIXY.

That amendment should remain documented separately and should not be rewritten into this historical SNP variant-calling pipeline.

See:

```text
99_PIXY_allsites_correction.md
```

for the separate corrective record.


---

## Outputs and QC checks — repository checklist

- Final VCF: `/parallel_scratch/lw01165/05_analysis/05b_146_data/Ascaris.autosomes.SNPs.146.maxmiss0.90.vcf.gz`
- Samples: **146**.
- Total records: **13,484,735**; all verified to be biallelic SNPs.
- No multiallelic or non-SNP records in the final audited file.
- Mapping QC removed **3** samples; variant/sample QC removed **7**; 11 intermediate flagged samples were retained.
- Do not confuse the exploratory **135-sample** `Ascaris.autosomes.SNPs.clean.vcf.gz` branch with the published/thesis **146-sample** dataset.

## Provenance and limitations

- The full reconstructed code is intentionally preserved with its original analysis parameters rather than retroactively imposing different filters.
- Sample metadata (`sample_QC_list.csv`) were **manually curated**.
- The historical source code includes some duplicate/rerun submission variants. These do not imply that every alternative command was used for every sample.
- The complete historical environment/version logs should be preserved when available.
- The **all-sites/invariant-site PIXY amendment is explicitly separate**: [`99_PIXY_allsites_correction.md`](99_PIXY_allsites_correction.md). This corrected VCF was **not** used to generate the historical thesis PIXY outputs.
