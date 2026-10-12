# 00 — Pipeline overview and provenance audit

## Purpose

An index of recovered *Ascaris* whole-genome analysis code and a conservative audit of its completeness. **Not a claim that all historical jobs were independently rerun or formally validated.**

## Inputs

Source consists of the user's original annotated Markdown files, later supplied SLURM scripts, and final pipeline QC logs. Manually curated `sample_QC_list.csv`, fixed balanced sampling lists and PIXY population files are referenced by the historical scripts but are not attached to this documentation set.

## Workflow

The numbered workflow proceeds from raw WGS through preprocessing, mapping, joint SNP calling, QC, structure and differentiation. A separate 25-kb relative-CN branch starts from marked-duplicate BAMs. The all-sites amendment is isolated at `99`.

## Annotated scripts

The numbered files contain the analytical commands. This audit itself intentionally does **not** duplicate all executables. See [README.md](README.md) for the file-by-file sequence.

## Outputs and QC checks

| Area | Recovered evidence | Status |
|---|---|---|
| FASTQ, fastp, Deacon, read QC | `01` | Covered |
| Chromosome-level reference and GenMap | `02` | Covered |
| 25-kb ≥90% mappability BED generation | Identified earlier source name `qc_25kb_mappability_070826.sub`, not recovered here | **Source gap** |
| BWA-MEM / MarkDuplicates / mapping QC | `03` | Covered; original in-house batch launcher not fully recovered |
| Per-sample mosdepth + chromosome CN + figures | `04` | Covered |
| Genotype counts + independent QC + genotype plots | `04` | Covered, full source supplied |
| GATK GVCF + joint calling + SNP filtering | `05` | Covered through audited n=146 final VCF |
| Balanced n=65 VCF | `06` | Final subsetting script and historical selection logic recovered; preserve actual fixed-ID lists |
| PCA, NJ, ML, IBS, PERMANOVA, ADMIXTURE | `06` | Covered as historical code transcript; converted from escaped source |
| All-sites GATK joint-genotyping and PIXY FST/dXY/pi/Tajima's D | `07` | Historical command excerpt and PIXY analyses; complete original all-sites SLURM wrapper not in current archive; **mappability inconsistency** |
| Wahlund | `08` | Covered |
| Full n=146 AMOVA input | `08` | Covered, no seed recorded for `--thin-count` |
| Full n=146 AMOVA R model code | Paths/models described in prior history but exact original R files not collected | **Source gap** |
| Balanced n=65 AMOVA R models | `08` | Covered |
| Chromosome summaries and tips/arms/centres | `09`, `10` | Covered, downstream of historical PIXY |
| FST/dXY candidate region and genes | `11` | Covered; some standalone R blocks depend on previously set variables |
| Corrected PIXY input construction | `99` | Covered and QC-passed; **no downstream PIXY reanalysis** |

## Provenance and limitations

### Confirmed numeric checks

- 156 WGS samples originally processed; mapping-QC and variant-QC exclusions yielded final **n=146**.
- All **13,484,735** final nuclear SNP VCF records are biallelic SNPs.
- Corrected PIXY input includes **156,583,802** invariant + **13,484,410** SNP records = **170,068,212** total, n=146, and zero invariant/SNP positional overlaps.
- The final PIXY SNP branch retained 325 fewer SNPs than the original curated SNP file, equivalent to 0.0024% of original curated SNP records.
- Mappability masking reduced original all-sites records from 210,854,819 to 189,429,567 (21,425,252 removed, all record types). Do **not** report that entire reduction as invariant sites.

### Findings that must not be misreported

1. Historical corrected PIXY VCF preparation does not demonstrate that original π, dXY, FST, Tajima's D, chromosome patterns or candidate genes were unchanged. No complete corrected PIXY rerun is included.
2. The final n=146 SNP dataset deliberately retained eleven flagged samples; earlier **135-sample** `clean.vcf.gz` was exploratory.
3. Main PLINK LD pruning r²=**0.15**; balanced AMOVA preparation r²=**0.20**. Full AMOVA 100k preparation uses the existing n=146 pruned dataset (r²=0.15), not the same r² as the balanced-specific preparation.
4. `plink2 --recode vcf bgz` is visible in the historic full-AMOVA script, while PLINK2 documentation specifies `--export vcf bgz`. Preserve original code and label before reuse.
5. Original `11` candidate R block references `outdir` before assignment; the source may depend on objects existing in the R session. Do not insert a guessed directory and claim it was used originally.
6. Genotype-count `bcftools stats` `PSC` field positions agree with the current documented schema, but software version and final count-log checks remain advisable.
7. Relative CN estimated by depth is a **relative copy-number measure**, not formal validation of absolute ploidy or aneuploidy.

### What remains to be archived before calling the repository fully reproducible

- Full original n=146 AMOVA R files and `qc_25kb_mappability_070826.sub`.
- The actual **manual** metadata CSV plus fixed keep/exclude and final population lists; account for sample-name normalization.
- Original SLURM logs and program/environment versions, and hashes/manifests for inputs and key intermediate VCFs.
- Intermediate commands retained in earlier ad hoc sessions where required. This set is an auditable analysis narrative, not an execution-tested workflow manager.
