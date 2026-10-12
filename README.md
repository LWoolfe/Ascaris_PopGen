# Population genomics of *Ascaris* spp. — analysis code

## Purpose

This repository preserves the original documented read-processing, mapping, SNP calling, population-genomic and supplementary analyses used for the *Ascaris lumbricoides–A. suum* complex PhD thesis/paper. It is a **historical analytical record**, not a single push-button portable pipeline. Commands were developed on the University of Surrey HPC system and in RStudio; SLURM resource settings, HPC directories and separate R sessions reflect that environment.

The aim is **transparent reproducibility without changing what was actually done**. Where the complete original script has not been recovered, the relevant file identifies the gap instead of fabricating commands.

## Datasets

| Data stage | Sample/record summary |
|---|---|
| Original WGS samples | 156 (28 newly generated, 128 public SRA) |
| Mapping QC | 3 excluded, 153 samples remaining for variant-QC assessment |
| Variant/sample QC | 7 further excluded, **146 final samples**; 11 intermediate flagged samples retained |
| Final autosomal SNP callset | **13,484,735 biallelic SNPs, n=146**, ≥90% site call rate |
| Balanced sensitivity subset | **n=65**; retains fixed IDs and selected over-represented-group subsamples |
| Ploidy/CN branch | starts from 156 marked-duplicate BAMs; ≥5× primary, ≥10× sensitivity |
| Corrected PIXY input | 156,583,802 invariant + 13,484,410 SNP = **170,068,212** records, 146 samples |

**Do not confuse these sample universes.** In particular, `Ascaris.autosomes.SNPs.clean.vcf.gz` in older notes refers to an exploratory **135-sample** strict-keep branch; it is *not* the final n=146 analysis VCF.

## Repository organisation and running order

| File | What it documents |
|---|---|
| [`00_pipeline_overview_and_audit.md`](00_pipeline_overview_and_audit.md) | Provenance audit and outstanding code/documentation items |
| [`01_sample_download_preprocessing_and_QC.md`](01_sample_download_preprocessing_and_QC.md) | SRA retrieval, fastp, Deacon, FastQC/MultiQC |
| [`02_reference_preparation_and_mappability.md`](02_reference_preparation_and_mappability.md) | Primary-chromosome reference, GenMap, mask and QC |
| [`03_mapping_and_mapping_QC.md`](03_mapping_and_mapping_QC.md) | BWA-MEM, BAM sort/index, duplicate marking, depth/breadth |
| [`04_ploidy_CNV_and_genotype_QC.md`](04_ploidy_CNV_and_genotype_QC.md) | Separate 25-kb relative-CN and genotype-composition branches |
| [`05_variant_calling_and_SNP_QC.md`](05_variant_calling_and_SNP_QC.md) | HaplotypeCaller, joint genotyping, PASS/mappability, final 146 SNP VCF |
| [`06_population_structure.md`](06_population_structure.md) | Balanced n=65 construction, PLINK, PCA, IBS, PERMANOVA, NJ/ML, ADMIXTURE |
| [`07_population_diversity_and_differentiation.md`](07_population_diversity_and_differentiation.md) | **Historical** PIXY all-sites FST/dXY/pi/Tajima's D and bootstrap summaries |
| [`08_wahlund_effect_and_AMOVA.md`](08_wahlund_effect_and_AMOVA.md) | Wahlund check, full AMOVA input prep, balanced AMOVA models |
| [`09_chromosomal_patterns.md`](09_chromosomal_patterns.md) | Chromosome-level pi/Tajima's D and chromosome-size correlations |
| [`10_chromosome_tips_arms_and_centres.md`](10_chromosome_tips_arms_and_centres.md) | Tip/arm/centre and technical QC |
| [`11_high_differentiation_candidate_regions.md`](11_high_differentiation_candidate_regions.md) | Empirical differentiation outliers, merged intervals, annotation |
| [`99_PIXY_allsites_correction.md`](99_PIXY_allsites_correction.md) | **Separate final audit/amendment; not used for original thesis statistics** |

All files use the same **Purpose → Inputs → Workflow → Annotated scripts → Outputs and QC checks → Provenance and limitations** sections. Some retain further nested headings from the historical originals. Standalone scripts in one Markdown file must be executed **individually**, not concatenated.

## Core workflow

```text
Raw WGS → fastp → Deacon host depletion → BWA-MEM → MarkDuplicates
                                                    ├→ 25-kb mosdepth / relative CN
                                                    └→ mapping QC → HaplotypeCaller GVCFs
                                                       → CombineGVCFs → GenotypeGVCFs
                                                       → SNP hard filters → unique-mappability mask
                                                       → 19 autosomes + biallelic SNPs
                                                       → sample QC → final n=146 / 13.48M SNPs
                                                         ├→ balanced n=65 subset
                                                         ├→ PLINK/PCA/IBS/PERMANOVA/NJ/ML/ADMIXTURE
                                                         ├→ Wahlund/AMOVA
                                                         ├→ genotype-composition QC
                                                         └→ historical all-sites/PIXY windows
                                                            → chromosome patterns/tip QC
                                                            → FST/dXY outliers and candidate genes
```

**Manual metadata:** `sample_QC_list.csv` and relevant population/sample ID files were manually compiled/curated. They should be preserved **as inputs**, rather than described as automatically generated. Large FASTQ/BAM/VCF files are not distributed through this Markdown repository. A genuinely rerunnable release would also need the exact metadata, population files, exclusion lists, reference hashes and source input accession manifest.

## Confirmed important analysis settings

- Primary reference: *A. suum* ASM1343314v1 (`GCA_013433145.1`), 19 autosomes `CM024168.1`–`CM024186.1`, plus five X chromosomes for mapping.
- GenMap unique mappability: `K=35`, `E=1`; the main SNP callset was masked to these uniquely mappable intervals.
- SNP hard filtering: `QD < 2.0`, `FS > 30`, `SOR > 3`, `MQ < 40` marked; PASS retained.
- Population-structure PLINK LD pruning: `--indep-pairwise 50 10 0.15`.
- **Balanced AMOVA-specific** preparation instead used `--indep-pairwise 50 10 0.2`.
- Full n=146 AMOVA preparation used an **existing r²=0.15-pruned PLINK set**, followed by `--thin-count 100000`; no fixed random seed is visible in the recovered script.
- The 146-sample final SNP VCF was filtered for **≥90% called genotypes per SNP**, not a newly imposed individual DP/GQ mask.
- Original PIXY statistics use 10-kb windows; final genome-wide R bootstrap summaries specify 10,000 iterations (some candidate-region exploratory summaries use 1,000).

## The all-sites/PIXY audit and correction

**This must not be silently incorporated into the historical analysis.** Final audit found that the original all-sites/invariant-site PIXY branch was **not restricted to uniquely mappable regions in the same way as the curated SNP VCF**. A corrected input was subsequently built and QC-checked, containing **156,583,802 invariant sites + 13,484,410 SNPs across 146 samples**, with no invariant/SNP positional overlaps. Compared with the original curated final SNP VCF, the corrected PIXY SNP branch differed by **325 SNPs (0.0024%)**.

The historical all-sites VCF was produced by a separate GATK `GenotypeGVCFs --include-non-variant-sites` branch from the merged GVCF (command excerpt in `07`). The historical thesis PIXY statistics, downstream chromosomal patterns and candidate-region analyses were **not rerun** with the corrected VCF. The 325-SNP difference does **not** establish that pi/dXY or other statistics dependent on the original all-sites input were unaffected. Complete PIXY-derived reanalysis is reserved for publication. See [`99_PIXY_allsites_correction.md`](99_PIXY_allsites_correction.md).

