# 08 — Wahlund-effect checks and AMOVA

## Purpose

Wahlund-effect checks and AMOVA. This file archives the commands and analysis notes from the thesis/paper workflow.

## Inputs

LD-pruned SNPs, curated metadata, n=65 and n=146 AMOVA VCFs. Source filenames and absolute HPC paths are retained in the original script blocks below.

## Workflow

See the original annotated analysis sequence below. Distinct script blocks are independent jobs/R sessions, **not one executable script**.

## Annotated scripts

### 08A. Full n=146 AMOVA input construction (recovered original commands)

**Computational constraint (confirmed by analyst):** Full-dataset AMOVA was conducted on **100,000 SNPs selected from the LD-pruned n=146 dataset**, not all ~13.48 million SNPs. Attempting AMOVA on the complete SNP set was computationally impractical and did not run successfully. The subset was therefore a deliberate feasibility step. This is distinct from the **balanced n=65 AMOVA**, which also has its own 100,000-SNP preparation workflow below. Do not describe either AMOVA as using all SNPs or imply that the 100,000-SNP subset contained every retained SNP.


Uses the **existing full-population PLINK dataset pruned at r²=0.15** in `05d_plink_146`, followed by PLINK2 `--thin-count 100000` to produce an AMOVA-specific 100k subset. **No seed is specified in this recovered script**, so exact re-selection of 100k loci from scratch is not guaranteed.

> **Compatibility warning, not a retrospective change:** the recovered job uses `plink2 --recode vcf bgz`, whereas the documented PLINK2 export flag is `--export vcf bgz`. The `.vcf.gz` was reported to exist historically, so the script may reflect a mixed executable version or a command that required manual correction. Check the historical `.out`/`.err` before claiming this exact command executed successfully.

```bash
#!/bin/bash
#SBATCH --job-name=amova_146_100k
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=08:00:00
#SBATCH --output=amova_146_100k.%j.out
#SBATCH --error=amova_146_100k.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate NJ_env

PLINKDIR=/parallel_scratch/lw01165/05_analysis/05d_plink_146
OUTDIR=/parallel_scratch/lw01165/05_analysis/05i_pixy/amova_100k_dataset

mkdir -p ${OUTDIR}
cd ${OUTDIR}

echo "Preparing 100k SNP AMOVA dataset (146 samples)"
ls ${PLINKDIR}/Ascaris.146.pruned.bed \
   ${PLINKDIR}/Ascaris.146.pruned.bim \
   ${PLINKDIR}/Ascaris.146.pruned.fam

if [ ! -f Ascaris.146.pruned.100kSNPs.bed ]; then
  cp ${PLINKDIR}/Ascaris.146.pruned.bed .
  cp ${PLINKDIR}/Ascaris.146.pruned.bim .
  cp ${PLINKDIR}/Ascaris.146.pruned.fam .

  plink2 \
    --bfile Ascaris.146.pruned \
    --thin-count 100000 \
    --make-bed \
    --out Ascaris.146.pruned.100kSNPs
fi

if [ ! -f Ascaris.146.pruned.100kSNPs.vcf.gz ]; then
  plink2 \
    --bfile Ascaris.146.pruned.100kSNPs \
    --recode vcf bgz \
    --out Ascaris.146.pruned.100kSNPs
fi

if [ ! -f Ascaris.146.pruned.100kSNPs.vcf.gz.tbi ]; then
  tabix -p vcf Ascaris.146.pruned.100kSNPs.vcf.gz
fi

wc -l Ascaris.146.pruned.100kSNPs.bim
bcftools index -n Ascaris.146.pruned.100kSNPs.vcf.gz || true

```

### 08B. Full n=146 AMOVA R models (recovered original scripts)

These original source files close the previously recorded gap in the **full n=146 AMOVA R-model code**. They load the **100,000-SNP VCF derived from the full 146-sample dataset**, run `poppr.amova()` and assess variance components using `randtest(..., nrepet = 999)`. The separately documented balanced n=65 AMOVA is a sensitivity analysis, not the primary full-dataset analysis.

| Original script | Model | Inclusion rule |
|---|---|---|
| [`amova_region.R`](scripts/amova_full_146/amova_region.R) | `~region` | Geographic region must be recorded |
| [`amova_host.R`](scripts/amova_full_146/amova_host.R) | `~host` | Host must be recorded; does not require region |
| [`amova_region_country.R`](scripts/amova_full_146/amova_region_country.R) | `~region/country` | Both region and country must be recorded |
| [`amova_region_host.R`](scripts/amova_full_146/amova_region_host.R) | `~region/host` | Both region and host must be recorded |
| [`amova_africa.R`](scripts/amova_full_146/amova_africa.R) | `~host (Africa only)` | Must be in Africa, have a host assignment, and retain ≥2 host groups |

**Important source-provenance checks:**

- **Input VCF path inconsistency:** the original `amova_region.R` points to `05i_pixy/amova_146_100k_dataset/`, whereas the four other model scripts and the recovered full-dataset preparation job point to `05i_pixy/amova_100k_dataset/`. Both paths are preserved in the original scripts; the correct historical path should be confirmed against HPC outputs/logs before reproduction.
- **Metadata and actual model sample sizes:** the `n=146` designation describes the source genotype dataset. Each R model then excludes individuals without its required grouping metadata. In particular, manually curated sample `ERR068355` lacks country and region but has host information, so it may be included in the host-only model and excluded from region-based models. The scripts write the metadata subset used for their respective models.
- **Metadata column-name reproducibility:** the original uploaded metadata CSV has a BOM and trailing whitespace in the `sample_id ` header. Original R scripts call `rename(IID = sample_id)` without first trimming column names; the normalised repository copy has an exact `sample_id` column, but the historical HPC input's exact header is not independently verified. Do not silently amend the original R scripts.
- **Execution status:** original scripts were recovered and inspected; their existence does not prove execution success or demonstrate that the old output values were reproduced. Preserve historical results, logs and source VCFs when available.

The following code blocks preserve the **original source contents without modifying parameters or paths**. Standalone, unchanged `.R` files are also archived under [`scripts/amova_full_146/`](scripts/amova_full_146/README.md).

#### Region only — `amova_region.R`

```r
# =============================
# AMOVA: region only
# =============================

.libPaths("~/R/library")

suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(stringr)
  library(adegenet)
  library(poppr)
})

# -----------------------------
# File paths
# -----------------------------
vcf_file <- "/parallel_scratch/lw01165/05_analysis/05i_pixy/amova_146_100k_dataset/Ascaris.146.pruned.100kSNPs.vcf.gz"
meta_file <- "/parallel_scratch/lw01165/05_analysis/sample_QC_list.csv"
outdir    <- "/parallel_scratch/lw01165/05_analysis/05i_pixy/amova_results_region_only"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("======================================\n")
cat("AMOVA started: region only\n")
cat("Time: ", as.character(Sys.time()), "\n", sep = "")
cat("VCF: ", vcf_file, "\n", sep = "")
cat("Metadata: ", meta_file, "\n", sep = "")
cat("Output directory: ", outdir, "\n", sep = "")
cat("======================================\n\n")

# -----------------------------
# Read data
# -----------------------------
cat("Reading VCF...\n")
vcf <- read.vcfR(vcf_file)

cat("Reading metadata...\n")
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# -----------------------------
# Clean metadata
# -----------------------------
meta <- meta %>%
  rename(IID = sample_id) %>%
  mutate(
    IID     = str_trim(IID),
    region  = str_to_lower(str_trim(region)),
    country = str_to_title(str_to_lower(str_trim(country))),
    host    = str_to_lower(str_trim(host))
  ) %>%
  mutate(
    region  = na_if(region, ""),
    country = na_if(country, ""),
    host    = na_if(host, "")
  )

# -----------------------------
# Match metadata to VCF
# -----------------------------
vcf_samples <- colnames(vcf@gt)[-1]

cat("Checking sample mismatches...\n")
vcf_not_meta <- setdiff(vcf_samples, meta$IID)
meta_not_vcf <- setdiff(meta$IID, vcf_samples)

cat("Samples in VCF but not metadata: ", length(vcf_not_meta), "\n", sep = "")
if (length(vcf_not_meta) > 0) print(vcf_not_meta)

cat("Samples in metadata but not VCF: ", length(meta_not_vcf), "\n", sep = "")
if (length(meta_not_vcf) > 0) print(meta_not_vcf)

meta <- meta %>%
  filter(IID %in% vcf_samples)

meta <- meta[match(vcf_samples, meta$IID), ]

if (!all(meta$IID == vcf_samples)) {
  stop("ERROR: metadata and VCF sample order do not match after reordering.")
}

write.csv(
  meta,
  file = file.path(outdir, "clean_metadata_matched_to_vcf.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("\nMetadata aligned to VCF successfully.\n")
cat("Total matched samples: ", nrow(meta), "\n\n", sep = "")

cat("Sample counts by region:\n")
print(table(meta$region, useNA = "ifany"))

# -----------------------------
# Keep samples for this model
# -----------------------------
keep_region <- !is.na(meta$region)

meta_region <- meta[keep_region, ]

cat("\nSamples retained for region-only AMOVA: ", nrow(meta_region), "\n", sep = "")
cat("Regions included:\n")
print(table(meta_region$region))

write.csv(
  meta_region,
  file = file.path(outdir, "metadata_region_only_used.csv"),
  row.names = FALSE,
  quote = FALSE
)

# -----------------------------
# Convert VCF to genind/genclone
# -----------------------------
cat("\nConverting VCF to genind...\n")
genind_obj <- vcfR2genind(vcf)

rm(vcf)
gc()

if (!all(indNames(genind_obj) == meta$IID)) {
  stop("ERROR: genind individual names do not match metadata IID order.")
}

cat("Subsetting genind for retained samples...\n")
genind_region <- genind_obj[keep_region]

rm(genind_obj)
gc()

if (!all(indNames(genind_region) == meta_region$IID)) {
  stop("ERROR: subset genind names do not match metadata IID order.")
}

cat("Converting genind to genclone...\n")
genclone_region <- as.genclone(genind_region)

rm(genind_region)
gc()

strata(genclone_region) <- data.frame(
  region = factor(meta_region$region)
)

# -----------------------------
# Run AMOVA
# -----------------------------
cat("\n======================================\n")
cat("Running AMOVA: region only\n")
cat("======================================\n")

amova_region <- poppr.amova(genclone_region, ~region)
print(amova_region)

cat("\nRunning permutation test (999 replicates)...\n")
region_test <- randtest(amova_region, nrepet = 999)
print(region_test)

# -----------------------------
# Save outputs
# -----------------------------
write.table(
  capture.output(amova_region),
  file = file.path(outdir, "amova_region_only.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  capture.output(region_test),
  file = file.path(outdir, "amova_region_only_randtest.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

saveRDS(amova_region, file = file.path(outdir, "amova_region_only.rds"))
saveRDS(region_test,  file = file.path(outdir, "amova_region_only_randtest.rds"))

cat("\n======================================\n")
cat("AMOVA complete: region only\n")
cat("Time: ", as.character(Sys.time()), "\n", sep = "")
cat("Results written to: ", outdir, "\n", sep = "")
cat("======================================\n")
```

#### Host only — `amova_host.R`

```r
# =============================
# AMOVA: host
# =============================

.libPaths("~/R/library")

suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(stringr)
  library(adegenet)
  library(poppr)
})

# -----------------------------
# File paths
# -----------------------------
vcf_file  <- "/parallel_scratch/lw01165/05_analysis/05i_pixy/amova_100k_dataset/Ascaris.146.pruned.100kSNPs.vcf.gz"
meta_file <- "/parallel_scratch/lw01165/05_analysis/sample_QC_list.csv"
outdir    <- "/parallel_scratch/lw01165/05_analysis/05i_pixy/amova_results_host"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("======================================\n")
cat("AMOVA started: host\n")
cat("Time: ", as.character(Sys.time()), "\n", sep = "")
cat("VCF: ", vcf_file, "\n", sep = "")
cat("Metadata: ", meta_file, "\n", sep = "")
cat("Output directory: ", outdir, "\n", sep = "")
cat("======================================\n\n")

# -----------------------------
# Read data
# -----------------------------
cat("Reading VCF...\n")
vcf <- read.vcfR(vcf_file)

cat("Reading metadata...\n")
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# -----------------------------
# Clean metadata
# -----------------------------
meta <- meta %>%
  rename(IID = sample_id) %>%
  mutate(
    IID     = str_trim(IID),
    region  = str_to_lower(str_trim(region)),
    country = str_to_title(str_to_lower(str_trim(country))),
    host    = str_to_lower(str_trim(host))
  ) %>%
  mutate(
    region  = na_if(region, ""),
    country = na_if(country, ""),
    host    = na_if(host, "")
  )

# -----------------------------
# Match metadata to VCF
# -----------------------------
vcf_samples <- colnames(vcf@gt)[-1]

cat("Checking sample mismatches...\n")
vcf_not_meta <- setdiff(vcf_samples, meta$IID)
meta_not_vcf <- setdiff(meta$IID, vcf_samples)

cat("Samples in VCF but not metadata: ", length(vcf_not_meta), "\n", sep = "")
if (length(vcf_not_meta) > 0) print(vcf_not_meta)

cat("Samples in metadata but not VCF: ", length(meta_not_vcf), "\n", sep = "")
if (length(meta_not_vcf) > 0) print(meta_not_vcf)

meta <- meta %>%
  filter(IID %in% vcf_samples)

meta <- meta[match(vcf_samples, meta$IID), ]

if (!all(meta$IID == vcf_samples)) {
  stop("ERROR: metadata and VCF sample order do not match after reordering.")
}

write.csv(
  meta,
  file = file.path(outdir, "clean_metadata_matched_to_vcf.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("\nMetadata aligned to VCF successfully.\n")
cat("Total matched samples: ", nrow(meta), "\n\n", sep = "")

cat("Sample counts by host:\n")
print(table(meta$host, useNA = "ifany"))

# -----------------------------
# Keep samples for this model
# -----------------------------
keep_host <- !is.na(meta$host)

meta_host <- meta[keep_host, ]

cat("\nSamples retained for host AMOVA: ", nrow(meta_host), "\n", sep = "")
cat("Hosts included:\n")
print(table(meta_host$host))

write.csv(
  meta_host,
  file = file.path(outdir, "metadata_host_used.csv"),
  row.names = FALSE,
  quote = FALSE
)

# -----------------------------
# Convert VCF to genind/genclone
# -----------------------------
cat("\nConverting VCF to genind...\n")
genind_obj <- vcfR2genind(vcf)

rm(vcf)
gc()

if (!all(indNames(genind_obj) == meta$IID)) {
  stop("ERROR: genind individual names do not match metadata IID order.")
}

cat("Subsetting genind for retained samples...\n")
genind_host <- genind_obj[keep_host]

rm(genind_obj)
gc()

if (!all(indNames(genind_host) == meta_host$IID)) {
  stop("ERROR: subset genind names do not match metadata IID order.")
}

cat("Converting genind to genclone...\n")
genclone_host <- as.genclone(genind_host)

rm(genind_host)
gc()

strata(genclone_host) <- data.frame(
  host = factor(meta_host$host)
)

# -----------------------------
# Run AMOVA
# -----------------------------
cat("\n======================================\n")
cat("Running AMOVA: host\n")
cat("======================================\n")

amova_host <- poppr.amova(genclone_host, ~host)
print(amova_host)

cat("\nRunning permutation test (999 replicates)...\n")
host_test <- randtest(amova_host, nrepet = 999)
print(host_test)

# -----------------------------
# Save outputs
# -----------------------------
write.table(
  capture.output(amova_host),
  file = file.path(outdir, "amova_host.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  capture.output(host_test),
  file = file.path(outdir, "amova_host_randtest.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

saveRDS(amova_host, file = file.path(outdir, "amova_host.rds"))
saveRDS(host_test,  file = file.path(outdir, "amova_host_randtest.rds"))

cat("\n======================================\n")
cat("AMOVA complete: host\n")
cat("Time: ", as.character(Sys.time()), "\n", sep = "")
cat("Results written to: ", outdir, "\n", sep = "")
cat("======================================\n")
```

#### Region / country — `amova_region_country.R`

```r
# =============================
# AMOVA: region / country
# =============================

.libPaths("~/R/library")

suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(stringr)
  library(adegenet)
  library(poppr)
})

# -----------------------------
# File paths
# -----------------------------
vcf_file <- "/parallel_scratch/lw01165/05_analysis/05i_pixy/amova_100k_dataset/Ascaris.146.pruned.100kSNPs.vcf.gz"
meta_file <- "/parallel_scratch/lw01165/05_analysis/sample_QC_list.csv"
outdir    <- "/parallel_scratch/lw01165/05_analysis/05i_pixy/amova_results_region_country"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("======================================\n")
cat("AMOVA started: region / country\n")
cat("Time: ", as.character(Sys.time()), "\n", sep = "")
cat("VCF: ", vcf_file, "\n", sep = "")
cat("Metadata: ", meta_file, "\n", sep = "")
cat("Output directory: ", outdir, "\n", sep = "")
cat("======================================\n\n")

# -----------------------------
# Read data
# -----------------------------
cat("Reading VCF...\n")
vcf <- read.vcfR(vcf_file)

cat("Reading metadata...\n")
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# -----------------------------
# Clean metadata
# -----------------------------
meta <- meta %>%
  rename(IID = sample_id) %>%
  mutate(
    IID     = str_trim(IID),
    region  = str_to_lower(str_trim(region)),
    country = str_to_title(str_to_lower(str_trim(country))),
    host    = str_to_lower(str_trim(host))
  ) %>%
  mutate(
    region  = na_if(region, ""),
    country = na_if(country, ""),
    host    = na_if(host, "")
  )

# -----------------------------
# Match metadata to VCF
# -----------------------------
vcf_samples <- colnames(vcf@gt)[-1]

cat("Checking sample mismatches...\n")
vcf_not_meta <- setdiff(vcf_samples, meta$IID)
meta_not_vcf <- setdiff(meta$IID, vcf_samples)

cat("Samples in VCF but not metadata: ", length(vcf_not_meta), "\n", sep = "")
if (length(vcf_not_meta) > 0) print(vcf_not_meta)

cat("Samples in metadata but not VCF: ", length(meta_not_vcf), "\n", sep = "")
if (length(meta_not_vcf) > 0) print(meta_not_vcf)

meta <- meta %>%
  filter(IID %in% vcf_samples)

meta <- meta[match(vcf_samples, meta$IID), ]

if (!all(meta$IID == vcf_samples)) {
  stop("ERROR: metadata and VCF sample order do not match after reordering.")
}

write.csv(
  meta,
  file = file.path(outdir, "clean_metadata_matched_to_vcf.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("\nMetadata aligned to VCF successfully.\n")
cat("Total matched samples: ", nrow(meta), "\n\n", sep = "")

cat("Sample counts by region:\n")
print(table(meta$region, useNA = "ifany"))

cat("\nSample counts by country:\n")
print(table(meta$country, useNA = "ifany"))

# -----------------------------
# Keep samples for this model
# -----------------------------
keep_geo <- !is.na(meta$region) & !is.na(meta$country)

meta_geo <- meta[keep_geo, ]

cat("\nSamples retained for region/country AMOVA: ", nrow(meta_geo), "\n", sep = "")
cat("Regions included:\n")
print(table(meta_geo$region))
cat("Countries included:\n")
print(table(meta_geo$country))

write.csv(
  meta_geo,
  file = file.path(outdir, "metadata_region_country_used.csv"),
  row.names = FALSE,
  quote = FALSE
)

# -----------------------------
# Convert VCF to genind/genclone
# -----------------------------
cat("\nConverting VCF to genind...\n")
genind_obj <- vcfR2genind(vcf)

rm(vcf)
gc()

if (!all(indNames(genind_obj) == meta$IID)) {
  stop("ERROR: genind individual names do not match metadata IID order.")
}

cat("Subsetting genind for retained samples...\n")
genind_geo <- genind_obj[keep_geo]

rm(genind_obj)
gc()

if (!all(indNames(genind_geo) == meta_geo$IID)) {
  stop("ERROR: subset genind names do not match metadata IID order.")
}

cat("Converting genind to genclone...\n")
genclone_geo <- as.genclone(genind_geo)

rm(genind_geo)
gc()

strata(genclone_geo) <- data.frame(
  region  = factor(meta_geo$region),
  country = factor(meta_geo$country)
)

# -----------------------------
# Run AMOVA
# -----------------------------
cat("\n======================================\n")
cat("Running AMOVA: region / country\n")
cat("======================================\n")

amova_geo <- poppr.amova(genclone_geo, ~region/country)
print(amova_geo)

cat("\nRunning permutation test (999 replicates)...\n")
geo_test <- randtest(amova_geo, nrepet = 999)
print(geo_test)

# -----------------------------
# Save outputs
# -----------------------------
write.table(
  capture.output(amova_geo),
  file = file.path(outdir, "amova_region_country.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  capture.output(geo_test),
  file = file.path(outdir, "amova_region_country_randtest.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

saveRDS(amova_geo, file = file.path(outdir, "amova_region_country.rds"))
saveRDS(geo_test,  file = file.path(outdir, "amova_region_country_randtest.rds"))

cat("\n======================================\n")
cat("AMOVA complete: region / country\n")
cat("Time: ", as.character(Sys.time()), "\n", sep = "")
cat("Results written to: ", outdir, "\n", sep = "")
cat("======================================\n")
```

#### Region / host — `amova_region_host.R`

```r
# =============================
# AMOVA: region / host
# =============================

.libPaths("~/R/library")

suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(stringr)
  library(adegenet)
  library(poppr)
})

# -----------------------------
# File paths
# -----------------------------
vcf_file <- "/parallel_scratch/lw01165/05_analysis/05i_pixy/amova_100k_dataset/Ascaris.146.pruned.100kSNPs.vcf.gz"
meta_file <- "/parallel_scratch/lw01165/05_analysis/sample_QC_list.csv"
outdir    <- "/parallel_scratch/lw01165/05_analysis/05i_pixy/amova_results_region_host"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("======================================\n")
cat("AMOVA started: region / host\n")
cat("Time: ", as.character(Sys.time()), "\n", sep = "")
cat("VCF: ", vcf_file, "\n", sep = "")
cat("Metadata: ", meta_file, "\n", sep = "")
cat("Output directory: ", outdir, "\n", sep = "")
cat("======================================\n\n")

# -----------------------------
# Read data
# -----------------------------
cat("Reading VCF...\n")
vcf <- read.vcfR(vcf_file)

cat("Reading metadata...\n")
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# -----------------------------
# Clean metadata
# -----------------------------
meta <- meta %>%
  rename(IID = sample_id) %>%
  mutate(
    IID     = str_trim(IID),
    region  = str_to_lower(str_trim(region)),
    country = str_to_title(str_to_lower(str_trim(country))),
    host    = str_to_lower(str_trim(host))
  ) %>%
  mutate(
    region  = na_if(region, ""),
    country = na_if(country, ""),
    host    = na_if(host, "")
  )

# -----------------------------
# Match metadata to VCF
# -----------------------------
vcf_samples <- colnames(vcf@gt)[-1]

cat("Checking sample mismatches...\n")
vcf_not_meta <- setdiff(vcf_samples, meta$IID)
meta_not_vcf <- setdiff(meta$IID, vcf_samples)

cat("Samples in VCF but not metadata: ", length(vcf_not_meta), "\n", sep = "")
if (length(vcf_not_meta) > 0) print(vcf_not_meta)

cat("Samples in metadata but not VCF: ", length(meta_not_vcf), "\n", sep = "")
if (length(meta_not_vcf) > 0) print(meta_not_vcf)

meta <- meta %>%
  filter(IID %in% vcf_samples)

meta <- meta[match(vcf_samples, meta$IID), ]

if (!all(meta$IID == vcf_samples)) {
  stop("ERROR: metadata and VCF sample order do not match after reordering.")
}

write.csv(
  meta,
  file = file.path(outdir, "clean_metadata_matched_to_vcf.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("\nMetadata aligned to VCF successfully.\n")
cat("Total matched samples: ", nrow(meta), "\n\n", sep = "")

cat("Sample counts by region:\n")
print(table(meta$region, useNA = "ifany"))

cat("\nSample counts by host:\n")
print(table(meta$host, useNA = "ifany"))

# -----------------------------
# Keep samples for this model
# -----------------------------
keep_host_geo <- !is.na(meta$region) & !is.na(meta$host)

meta_host_geo <- meta[keep_host_geo, ]

cat("\nSamples retained for region/host AMOVA: ", nrow(meta_host_geo), "\n", sep = "")
cat("Regions included:\n")
print(table(meta_host_geo$region))
cat("Hosts included:\n")
print(table(meta_host_geo$host))

write.csv(
  meta_host_geo,
  file = file.path(outdir, "metadata_region_host_used.csv"),
  row.names = FALSE,
  quote = FALSE
)

# -----------------------------
# Convert VCF to genind/genclone
# -----------------------------
cat("\nConverting VCF to genind...\n")
genind_obj <- vcfR2genind(vcf)

rm(vcf)
gc()

if (!all(indNames(genind_obj) == meta$IID)) {
  stop("ERROR: genind individual names do not match metadata IID order.")
}

cat("Subsetting genind for retained samples...\n")
genind_host_geo <- genind_obj[keep_host_geo]

rm(genind_obj)
gc()

if (!all(indNames(genind_host_geo) == meta_host_geo$IID)) {
  stop("ERROR: subset genind names do not match metadata IID order.")
}

cat("Converting genind to genclone...\n")
genclone_host_geo <- as.genclone(genind_host_geo)

rm(genind_host_geo)
gc()

strata(genclone_host_geo) <- data.frame(
  region = factor(meta_host_geo$region),
  host   = factor(meta_host_geo$host)
)

# -----------------------------
# Run AMOVA
# -----------------------------
cat("\n======================================\n")
cat("Running AMOVA: region / host\n")
cat("======================================\n")

amova_host_geo <- poppr.amova(genclone_host_geo, ~region/host)
print(amova_host_geo)

cat("\nRunning permutation test (999 replicates)...\n")
host_geo_test <- randtest(amova_host_geo, nrepet = 999)
print(host_geo_test)

# -----------------------------
# Save outputs
# -----------------------------
write.table(
  capture.output(amova_host_geo),
  file = file.path(outdir, "amova_region_host.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  capture.output(host_geo_test),
  file = file.path(outdir, "amova_region_host_randtest.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

saveRDS(amova_host_geo, file = file.path(outdir, "amova_region_host.rds"))
saveRDS(host_geo_test,  file = file.path(outdir, "amova_region_host_randtest.rds"))

cat("\n======================================\n")
cat("AMOVA complete: region / host\n")
cat("Time: ", as.character(Sys.time()), "\n", sep = "")
cat("Results written to: ", outdir, "\n", sep = "")
cat("======================================\n")
```

#### Host within Africa — `amova_africa.R`

```r
# =============================
# AMOVA: within Africa, host only
# =============================

.libPaths("~/R/library")

suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(stringr)
  library(adegenet)
  library(poppr)
})

# -----------------------------
# File paths
# -----------------------------
vcf_file <- "/parallel_scratch/lw01165/05_analysis/05i_pixy/amova_100k_dataset/Ascaris.146.pruned.100kSNPs.vcf.gz"
meta_file <- "/parallel_scratch/lw01165/05_analysis/sample_QC_list.csv"
outdir    <- "/parallel_scratch/lw01165/05_analysis/05i_pixy/amova_results_africa_host"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("======================================\n")
cat("AMOVA started: within Africa / host\n")
cat("Time: ", as.character(Sys.time()), "\n", sep = "")
cat("VCF: ", vcf_file, "\n", sep = "")
cat("Metadata: ", meta_file, "\n", sep = "")
cat("Output directory: ", outdir, "\n", sep = "")
cat("======================================\n\n")

# -----------------------------
# Read data
# -----------------------------
cat("Reading VCF...\n")
vcf <- read.vcfR(vcf_file)

cat("Reading metadata...\n")
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# -----------------------------
# Clean metadata
# -----------------------------
meta <- meta %>%
  rename(IID = sample_id) %>%
  mutate(
    IID     = str_trim(IID),
    region  = str_to_lower(str_trim(region)),
    country = str_to_title(str_to_lower(str_trim(country))),
    host    = str_to_lower(str_trim(host))
  ) %>%
  mutate(
    region  = na_if(region, ""),
    country = na_if(country, ""),
    host    = na_if(host, "")
  )

# -----------------------------
# Match metadata to VCF
# -----------------------------
vcf_samples <- colnames(vcf@gt)[-1]

cat("Checking sample mismatches...\n")
vcf_not_meta <- setdiff(vcf_samples, meta$IID)
meta_not_vcf <- setdiff(meta$IID, vcf_samples)

cat("Samples in VCF but not metadata: ", length(vcf_not_meta), "\n", sep = "")
if (length(vcf_not_meta) > 0) print(vcf_not_meta)

cat("Samples in metadata but not VCF: ", length(meta_not_vcf), "\n", sep = "")
if (length(meta_not_vcf) > 0) print(meta_not_vcf)

meta <- meta %>%
  filter(IID %in% vcf_samples)

meta <- meta[match(vcf_samples, meta$IID), ]

if (!all(meta$IID == vcf_samples)) {
  stop("ERROR: metadata and VCF sample order do not match after reordering.")
}

write.csv(
  meta,
  file = file.path(outdir, "clean_metadata_matched_to_vcf.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("\nMetadata aligned to VCF successfully.\n")
cat("Total matched samples: ", nrow(meta), "\n\n", sep = "")

cat("Sample counts by region:\n")
print(table(meta$region, useNA = "ifany"))

cat("\nSample counts by host:\n")
print(table(meta$host, useNA = "ifany"))

# -----------------------------
# Keep Africa samples with host metadata
# -----------------------------
keep_africa_host <- !is.na(meta$region) & !is.na(meta$host) & meta$region == "africa"

meta_africa_host <- meta[keep_africa_host, ]

cat("\nSamples retained for within-Africa host AMOVA: ", nrow(meta_africa_host), "\n", sep = "")
cat("Countries included:\n")
print(table(meta_africa_host$country, useNA = "ifany"))
cat("Hosts included:\n")
print(table(meta_africa_host$host, useNA = "ifany"))

write.csv(
  meta_africa_host,
  file = file.path(outdir, "metadata_africa_host_used.csv"),
  row.names = FALSE,
  quote = FALSE
)

# Optional safety check
if (length(unique(meta_africa_host$host)) < 2) {
  stop("ERROR: Fewer than 2 host groups remain within Africa. Cannot run host AMOVA.")
}

# -----------------------------
# Convert VCF to genind/genclone
# -----------------------------
cat("\nConverting VCF to genind...\n")
genind_obj <- vcfR2genind(vcf)

rm(vcf)
gc()

if (!all(indNames(genind_obj) == meta$IID)) {
  stop("ERROR: genind individual names do not match metadata IID order.")
}

cat("Subsetting genind for retained Africa samples...\n")
genind_africa_host <- genind_obj[keep_africa_host]

rm(genind_obj)
gc()

if (!all(indNames(genind_africa_host) == meta_africa_host$IID)) {
  stop("ERROR: subset genind names do not match metadata IID order.")
}

cat("Converting genind to genclone...\n")
genclone_africa_host <- as.genclone(genind_africa_host)

rm(genind_africa_host)
gc()

strata(genclone_africa_host) <- data.frame(
  host = factor(meta_africa_host$host)
)

# -----------------------------
# Run AMOVA
# -----------------------------
cat("\n======================================\n")
cat("Running AMOVA: within Africa / host\n")
cat("======================================\n")

amova_africa_host <- poppr.amova(genclone_africa_host, ~host)
print(amova_africa_host)

cat("\nRunning permutation test (999 replicates)...\n")
africa_host_test <- randtest(amova_africa_host, nrepet = 999)
print(africa_host_test)

# -----------------------------
# Save outputs
# -----------------------------
write.table(
  capture.output(amova_africa_host),
  file = file.path(outdir, "amova_africa_host.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  capture.output(africa_host_test),
  file = file.path(outdir, "amova_africa_host_randtest.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

saveRDS(amova_africa_host, file = file.path(outdir, "amova_africa_host.rds"))
saveRDS(africa_host_test,  file = file.path(outdir, "amova_africa_host_randtest.rds"))

cat("\n======================================\n")
cat("AMOVA complete: within Africa / host\n")
cat("Time: ", as.character(Sys.time()), "\n", sep = "")
cat("Results written to: ", outdir, "\n", sep = "")
cat("======================================\n")
```

---

### 08C. Historical Wahlund and balanced n=65 AMOVA scripts

# Supplementary heterozygosity-pooling check and AMOVA

This annotated section covers two separate but related checks:

1. **Wahlund-effect / heterozygosity pooling check**: tests whether pooling countries into broad regions inflates expected heterozygosity.
2. **AMOVA**: partitions genetic variance among hierarchical groups such as region/country and region/host.

These are supplementary because both address the same broader concern: host and geography are confounded in the sampling design, and broad groups may contain substructure.


## Wahlund-effect check: make country-level keep files

Purpose: create PLINK keep files for countries within regions. The main analysis uses countries with n ≥ 5; the n ≥ 2 version is descriptive only.

```r
# ============================================================
# Make country-level PLINK keep files for Wahlund-effect check
# ============================================================

setwd("S:/Research/Lauren_PopGen/Wahlund_effect")

base_dir <- "S:/Research/Lauren_PopGen/Wahlund_effect"

metadata_file <- file.path(base_dir, "metadata_146.csv")
fam_file      <- file.path(base_dir, "Ascaris.146.pruned.fam")

out_dir <- file.path(base_dir, "wahlund_country_keep")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Read input files
# -----------------------------

meta <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

fam <- read.table(
  fam_file,
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(fam) <- c("FID", "IID", "father", "mother", "sex", "phenotype")

# -----------------------------
# Clean metadata
# -----------------------------

# Clean column names
colnames(meta) <- trimws(colnames(meta))

# Check required columns exist
required_cols <- c("sample_id", "country", "region", "host")
missing_cols <- setdiff(required_cols, colnames(meta))

if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing required columns in metadata:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# Trim whitespace from all character columns
for (col in colnames(meta)) {
  if (is.character(meta[[col]])) {
    meta[[col]] <- trimws(meta[[col]])
  }
}

# Remove unknown/empty country or region rows
meta <- meta[
  !(meta$country %in% c("/", "", NA)) &
    !(meta$region %in% c("/", "", NA)),
]

# -----------------------------
# Match metadata to PLINK .fam
# -----------------------------

meta2 <- merge(
  meta,
  fam[, c("FID", "IID")],
  by.x = "sample_id",
  by.y = "IID",
  all.x = TRUE
)

missing_samples <- meta2[is.na(meta2$FID), ]

if (nrow(missing_samples) > 0) {
  print(missing_samples[, c("sample_id", "country", "region", "host")])
  stop("Some metadata sample IDs did not match the PLINK .fam IID column.")
}

# -----------------------------
# Count samples by region/country
# -----------------------------

country_counts <- as.data.frame(
  table(meta2$region, meta2$country),
  stringsAsFactors = FALSE
)

colnames(country_counts) <- c("region", "country", "N")
country_counts <- country_counts[country_counts$N > 0, ]

country_counts <- country_counts[order(country_counts$region,
                                       country_counts$country), ]

write.table(
  country_counts,
  file.path(base_dir, "wahlund_country_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print("Country counts:")
print(country_counts)

# -----------------------------
# Select eligible countries
# -----------------------------

# Main analysis: only countries with n >= 5
min_n <- 5

eligible <- country_counts[country_counts$N >= min_n, ]

# Keep only regions with at least two eligible countries
region_counts <- as.data.frame(
  table(eligible$region),
  stringsAsFactors = FALSE
)

colnames(region_counts) <- c("region", "n_eligible_countries")

eligible_regions <- region_counts$region[
  region_counts$n_eligible_countries >= 2
]

eligible <- eligible[eligible$region %in% eligible_regions, ]

print("Eligible countries for main Wahlund analysis:")
print(eligible)

# -----------------------------
# Make PLINK keep files
# -----------------------------

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_]", "_", x)
}

group_table <- data.frame(
  group_id = character(),
  region = character(),
  country = character(),
  n_samples = integer(),
  keep_file = character(),
  stringsAsFactors = FALSE
)

if (nrow(eligible) == 0) {
  warning("No eligible regions/countries found using min_n = 5.")
} else {
  
  for (i in seq_len(nrow(eligible))) {
    
    reg <- eligible$region[i]
    cty <- eligible$country[i]
    
    group_id <- paste(
      safe_name(reg),
      safe_name(cty),
      sep = "__"
    )
    
    samples <- meta2[
      meta2$region == reg & meta2$country == cty,
      c("FID", "sample_id")
    ]
    
    colnames(samples) <- c("FID", "IID")
    
    keep_file <- file.path(out_dir, paste0(group_id, ".keep"))
    
    write.table(
      samples,
      keep_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    
    group_table <- rbind(
      group_table,
      data.frame(
        group_id = group_id,
        region = reg,
        country = cty,
        n_samples = nrow(samples),
        keep_file = keep_file,
        stringsAsFactors = FALSE
      )
    )
  }
}

write.table(
  group_table,
  file.path(base_dir, "wahlund_country_groups_n5.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print("Final group table:")
print(group_table)

print("Done. Keep files written to:")
print(out_dir)

# ============================================================
# Additional descriptive-only keep files: countries with n >= 2
# ============================================================

out_dir_n2 <- file.path(base_dir, "wahlund_country_keep_n2_descriptive")
dir.create(out_dir_n2, showWarnings = FALSE, recursive = TRUE)

min_n_desc <- 2

eligible_desc <- country_counts[country_counts$N >= min_n_desc, ]

# Keep only regions with at least two eligible countries
region_counts_desc <- as.data.frame(
  table(eligible_desc$region),
  stringsAsFactors = FALSE
)

colnames(region_counts_desc) <- c("region", "n_eligible_countries")

eligible_regions_desc <- region_counts_desc$region[
  region_counts_desc$n_eligible_countries >= 2
]

eligible_desc <- eligible_desc[
  eligible_desc$region %in% eligible_regions_desc,
]

print("Eligible countries for descriptive-only Wahlund analysis, n >= 2:")
print(eligible_desc)

group_table_desc <- data.frame(
  group_id = character(),
  region = character(),
  country = character(),
  n_samples = integer(),
  keep_file = character(),
  analysis_type = character(),
  stringsAsFactors = FALSE
)

if (nrow(eligible_desc) == 0) {
  warning("No eligible regions/countries found using min_n = 2.")
} else {
  
  for (i in seq_len(nrow(eligible_desc))) {
    
    reg <- eligible_desc$region[i]
    cty <- eligible_desc$country[i]
    
    group_id <- paste(
      safe_name(reg),
      safe_name(cty),
      sep = "__"
    )
    
    samples <- meta2[
      meta2$region == reg & meta2$country == cty,
      c("FID", "sample_id")
    ]
    
    colnames(samples) <- c("FID", "IID")
    
    keep_file <- file.path(out_dir_n2, paste0(group_id, ".keep"))
    
    write.table(
      samples,
      keep_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    
    group_table_desc <- rbind(
      group_table_desc,
      data.frame(
        group_id = group_id,
        region = reg,
        country = cty,
        n_samples = nrow(samples),
        keep_file = keep_file,
        analysis_type = "descriptive_n2",
        stringsAsFactors = FALSE
      )
    )
  }
}

write.table(
  group_table_desc,
  file.path(base_dir, "wahlund_country_groups_n2_descriptive.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print("Final descriptive-only group table, n >= 2:")
print(group_table_desc)

print("Descriptive n >= 2 keep files written to:")
print(out_dir_n2)
```


## Cluster code: PLINK allele counts for Wahlund check

# Purpose: use the LD-pruned full dataset to generate allele counts for each eligible country group.

```bash
#!/bin/bash
#SBATCH --job-name=wahlund_acount
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --time=06:00:00
#SBATCH --output=wahlund_acount.%j.out
#SBATCH --error=wahlund_acount.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate NJ_env

BASE=/parallel_scratch/lw01165/05_analysis/05d_plink_146
BFILE=${BASE}/Ascaris.146.pruned

echo "Starting Wahlund allele-count analysis"
date

echo "Using PLINK files:"
ls -lh ${BFILE}.bed ${BFILE}.bim ${BFILE}.fam

run_acount () {

  GROUP_FILE=$1
  KEEP_DIR=$2
  OUTDIR=$3
  LABEL=$4

  echo "----------------------------------------"
  echo "Running: ${LABEL}"
  echo "Group file: ${GROUP_FILE}"
  echo "Keep dir: ${KEEP_DIR}"
  echo "Output dir: ${OUTDIR}"
  echo "----------------------------------------"

  if [[ ! -f "${GROUP_FILE}" ]]; then
    echo "ERROR: group file not found: ${GROUP_FILE}"
    exit 1
  fi

  if [[ ! -d "${KEEP_DIR}" ]]; then
    echo "ERROR: keep directory not found: ${KEEP_DIR}"
    exit 1
  fi

  mkdir -p "${OUTDIR}"

  tail -n +2 "${GROUP_FILE}" | while IFS=$'\t' read -r group_id region country n_samples keep_file extra; do

    if [[ -z "${group_id}" ]]; then
      continue
    fi

    KEEP=${KEEP_DIR}/${group_id}.keep

    if [[ ! -f "${KEEP}" ]]; then
      echo "ERROR: keep file not found for ${group_id}: ${KEEP}"
      exit 1
    fi

    echo "Running PLINK allele counts for ${group_id} | ${region} | ${country} | n=${n_samples}"

    plink2 \
      --bfile "${BFILE}" \
      --allow-extra-chr \
      --keep "${KEEP}" \
      --freq counts \
      --out "${OUTDIR}/${group_id}"

  done

  echo "Finished ${LABEL}"
  date
}

# Main analysis: countries with n >= 5
run_acount \
  "${BASE}/wahlund_country_groups_n5.tsv" \
  "${BASE}/wahlund_country_keep" \
  "${BASE}/wahlund_country_freq_n5" \
  "main n >= 5 Wahlund check"

# Descriptive-only analysis: countries with n >= 2
run_acount \
  "${BASE}/wahlund_country_groups_n2_descriptive.tsv" \
  "${BASE}/wahlund_country_keep_n2_descriptive" \
  "${BASE}/wahlund_country_freq_n2_descriptive" \
  "descriptive n >= 2 Wahlund check"

echo "All PLINK allele counts completed"
date
```


## RStudio code: calculate HT, HS and proportional pooling inflation

Purpose: calculate pooled expected heterozygosity (HT), weighted within-country expected heterozygosity (HS), the difference HT − HS, and proportional inflation.

```r

# ============================================================
# After PLINK run in cluster. Wahlund-effect HT/HS analysis
# ============================================================

setwd("S:/Research/Lauren_PopGen/Wahlund_effect")

base_dir <- "S:/Research/Lauren_PopGen/Wahlund_effect"

min_call_rate <- 0.80

# ------------------------------------------------------------
# Helper function: read one PLINK2 .acount file
# ------------------------------------------------------------

read_acount_file <- function(file, group_id, region, country, n_samples) {
  
  if (!file.exists(file)) {
    stop(paste("Missing .acount file:", file))
  }
  
  x <- read.table(
    file,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""   # IMPORTANT: prevents R skipping #CHROM header
  )
  
  # PLINK2 usually names first column "#CHROM"
  if ("#CHROM" %in% colnames(x)) {
    colnames(x)[colnames(x) == "#CHROM"] <- "CHROM"
  }
  
  print(paste("Reading:", file))
  print(colnames(x))
  
  required_cols <- c("CHROM", "ID", "ALT_CTS", "OBS_CT")
  missing_cols <- setdiff(required_cols, colnames(x))
  
  if (length(missing_cols) > 0) {
    stop(
      paste(
        "Missing columns in", file, ":",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  x <- x[, required_cols]
  
  x$ALT_CTS <- as.numeric(x$ALT_CTS)
  x$OBS_CT  <- as.numeric(x$OBS_CT)
  x$n_samples <- as.numeric(n_samples)
  
  # Remove sites with no observations
  x <- x[x$OBS_CT > 0, ]
  
  # Require minimum call rate within this country group
  min_obs_required <- 2 * as.numeric(n_samples) * min_call_rate
  x <- x[x$OBS_CT >= min_obs_required, ]
  
  x$p_alt <- x$ALT_CTS / x$OBS_CT
  x$HE_country <- 2 * x$p_alt * (1 - x$p_alt)
  
  suffix <- paste(region, country, sep = "__")
  suffix <- gsub("[^A-Za-z0-9_]", "_", suffix)
  
  colnames(x)[colnames(x) == "ALT_CTS"]    <- paste0("ALT_CTS__", suffix)
  colnames(x)[colnames(x) == "OBS_CT"]     <- paste0("OBS_CT__", suffix)
  colnames(x)[colnames(x) == "HE_country"] <- paste0("HE_country__", suffix)
  
  x <- x[, c(
    "CHROM",
    "ID",
    paste0("ALT_CTS__", suffix),
    paste0("OBS_CT__", suffix),
    paste0("HE_country__", suffix)
  )]
  
  return(x)
}

# ------------------------------------------------------------
# Main function: calculate HT, HS and pooling inflation
# ------------------------------------------------------------

calculate_wahlund <- function(groups_file,
                              freq_dir,
                              out_prefix,
                              analysis_label,
                              min_call_rate = 0.80,
                              write_per_snp = TRUE) {
  
  cat("\n========================================\n")
  cat("Running:", analysis_label, "\n")
  cat("Groups file:", groups_file, "\n")
  cat("Frequency directory:", freq_dir, "\n")
  cat("========================================\n")
  
  if (!file.exists(groups_file)) {
    stop(paste("Groups file not found:", groups_file))
  }
  
  if (!dir.exists(freq_dir)) {
    stop(paste("Frequency directory not found:", freq_dir))
  }
  
  groups <- read.table(
    groups_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Clean whitespace
  for (col in colnames(groups)) {
    if (is.character(groups[[col]])) {
      groups[[col]] <- trimws(groups[[col]])
    }
  }
  
  required_group_cols <- c("group_id", "region", "country", "n_samples")
  missing_group_cols <- setdiff(required_group_cols, colnames(groups))
  
  if (length(missing_group_cols) > 0) {
    stop(
      paste(
        "Missing columns in group file:",
        paste(missing_group_cols, collapse = ", ")
      )
    )
  }
  
  groups$n_samples <- as.numeric(groups$n_samples)
  
  regions <- unique(groups$region)
  
  all_region_summaries <- list()
  all_per_snp <- list()
  
  for (reg in regions) {
    
    cat("\nProcessing region:", reg, "\n")
    
    reg_groups <- groups[groups$region == reg, ]
    
    if (nrow(reg_groups) < 2) {
      warning(paste("Skipping", reg, "- fewer than two countries."))
      next
    }
    
    region_data_list <- list()
    
    for (i in seq_len(nrow(reg_groups))) {
      
      group_id <- reg_groups$group_id[i]
      country <- reg_groups$country[i]
      n_samples <- reg_groups$n_samples[i]
      
      acount_file <- file.path(freq_dir, paste0(group_id, ".acount"))
      
      cat("  Reading:", group_id, "|", country, "| n =", n_samples, "\n")
      
      region_data_list[[group_id]] <- read_acount_file(
        file = acount_file,
        group_id = group_id,
        region = reg,
        country = country,
        n_samples = n_samples
      )
    }
    
    # Merge all countries in this region by CHROM and ID
    merged <- Reduce(
      function(x, y) merge(x, y, by = c("CHROM", "ID"), all = FALSE),
      region_data_list
    )
    
    cat("  SNPs retained after requiring all countries:", nrow(merged), "\n")
    
    if (nrow(merged) == 0) {
      warning(paste("No SNPs retained for region:", reg))
      next
    }
    
    alt_cols <- grep("^ALT_CTS__", colnames(merged), value = TRUE)
    obs_cols <- grep("^OBS_CT__", colnames(merged), value = TRUE)
    he_cols  <- grep("^HE_country__", colnames(merged), value = TRUE)
    
    total_alt <- rowSums(merged[, alt_cols, drop = FALSE], na.rm = TRUE)
    total_obs <- rowSums(merged[, obs_cols, drop = FALSE], na.rm = TRUE)
    
    p_region <- total_alt / total_obs
    HT <- 2 * p_region * (1 - p_region)
    
    # HS = weighted mean country HE, weighted by observed allele counts
    obs_matrix <- as.matrix(merged[, obs_cols, drop = FALSE])
    he_matrix  <- as.matrix(merged[, he_cols, drop = FALSE])
    
    weights <- obs_matrix / rowSums(obs_matrix, na.rm = TRUE)
    HS <- rowSums(weights * he_matrix, na.rm = TRUE)
    
    delta_HT_HS <- HT - HS
    pooling_inflation <- ifelse(HT > 0, delta_HT_HS / HT, NA)
    
    per_snp <- data.frame(
      analysis = analysis_label,
      region = reg,
      CHROM = merged$CHROM,
      ID = merged$ID,
      HT = HT,
      HS = HS,
      delta_HT_HS = delta_HT_HS,
      pooling_inflation = pooling_inflation,
      stringsAsFactors = FALSE
    )
    
    # Remove pooled monomorphic sites
    per_snp <- per_snp[per_snp$HT > 0, ]
    
    countries_included <- paste(
      paste0(reg_groups$country, " (n=", reg_groups$n_samples, ")"),
      collapse = "; "
    )
    
    summary <- data.frame(
      analysis = analysis_label,
      region = reg,
      countries_included = countries_included,
      total_samples = sum(reg_groups$n_samples),
      n_snps = nrow(per_snp),
      mean_HT = mean(per_snp$HT, na.rm = TRUE),
      mean_HS = mean(per_snp$HS, na.rm = TRUE),
      mean_delta_HT_HS = mean(per_snp$delta_HT_HS, na.rm = TRUE),
      proportional_inflation =
        (mean(per_snp$HT, na.rm = TRUE) - mean(per_snp$HS, na.rm = TRUE)) /
        mean(per_snp$HT, na.rm = TRUE),
      mean_locus_inflation = mean(per_snp$pooling_inflation, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    all_region_summaries[[reg]] <- summary
    all_per_snp[[reg]] <- per_snp
    
    print(summary)
  }
  
  final_summary <- do.call(rbind, all_region_summaries)
  
  summary_csv <- file.path(base_dir, paste0(out_prefix, "_summary.csv"))
  summary_tsv <- file.path(base_dir, paste0(out_prefix, "_summary.tsv"))
  
  write.csv(final_summary, summary_csv, row.names = FALSE)
  write.table(
    final_summary,
    summary_tsv,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  cat("\nWrote summary files:\n")
  cat(summary_csv, "\n")
  cat(summary_tsv, "\n")
  
  if (write_per_snp) {
    
    final_per_snp <- do.call(rbind, all_per_snp)
    
    per_snp_csv <- file.path(base_dir, paste0(out_prefix, "_per_snp.csv"))
    per_snp_tsv <- file.path(base_dir, paste0(out_prefix, "_per_snp.tsv"))
    
    write.csv(final_per_snp, per_snp_csv, row.names = FALSE)
    write.table(
      final_per_snp,
      per_snp_tsv,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    cat("\nWrote per-SNP files:\n")
    cat(per_snp_csv, "\n")
    cat(per_snp_tsv, "\n")
  }
  
  return(final_summary)
}

# ------------------------------------------------------------
# Run main n >= 5 analysis
# ------------------------------------------------------------

summary_n5 <- calculate_wahlund(
  groups_file = file.path(base_dir, "wahlund_country_groups_n5.tsv"),
  freq_dir = file.path(base_dir, "wahlund_country_freq_n5"),
  out_prefix = "wahlund_n5_main",
  analysis_label = "Main country-level check: n >= 5",
  min_call_rate = min_call_rate,
  write_per_snp = TRUE
)

# ------------------------------------------------------------
# Run descriptive n >= 2 analysis
# ------------------------------------------------------------

summary_n2 <- calculate_wahlund(
  groups_file = file.path(base_dir, "wahlund_country_groups_n2_descriptive.tsv"),
  freq_dir = file.path(base_dir, "wahlund_country_freq_n2_descriptive"),
  out_prefix = "wahlund_n2_descriptive",
  analysis_label = "Descriptive sensitivity check: n >= 2",
  min_call_rate = min_call_rate,
  write_per_snp = TRUE
)

# ------------------------------------------------------------
# Combine summaries
# ------------------------------------------------------------

combined_summary <- rbind(summary_n5, summary_n2)

combined_csv <- file.path(base_dir, "wahlund_combined_summary.csv")
combined_tsv <- file.path(base_dir, "wahlund_combined_summary.tsv")

write.csv(combined_summary, combined_csv, row.names = FALSE)
write.table(
  combined_summary,
  combined_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nCombined summary written to:\n")
cat(combined_csv, "\n")
cat(combined_tsv, "\n")

print(combined_summary)

cat("\nDone.\n")
```

# ------------------------------------------------------------
## AMOVA input preparation and 100k SNP selection

Purpose: prepare a 100,000-SNP VCF for AMOVA from SNP data.
# ------------------------------------------------------------

```bash
#!/bin/bash
#SBATCH --job-name=prep_amova65_100k
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --output=prep_amova65_100k.%j.out
#SBATCH --error=prep_amova65_100k.%j.err

set -eo pipefail

module purge
source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate NJ_env

INDIR=/parallel_scratch/lw01165/05_analysis/05a_Balanced_dataset
OUTDIR=/parallel_scratch/lw01165/05_analysis/05i_pixy/amova_65_100k_dataset
VCF=${INDIR}/Ascaris.autosomes.SNPs.146.maxmiss0.90.balanced.n65.vcf.gz

mkdir -p ${OUTDIR}
cd ${OUTDIR}

echo "Starting balanced 65-sample 100k SNP preparation"
date
echo "Input VCF: ${VCF}"

# Convert VCF to PLINK
plink --vcf ${VCF} \
      --double-id \
      --allow-extra-chr \
      --set-missing-var-ids @:# \
      --make-bed \
      --out Ascaris.65

# LD pruning
plink --bfile Ascaris.65 \
      --allow-extra-chr \
      --indep-pairwise 50 10 0.2 \
      --out Ascaris.65.pruned

# Count available pruned SNPs
echo "Number of pruned SNPs available:"
wc -l Ascaris.65.pruned.prune.in

# Extract all pruned SNPs into a bed set
plink --bfile Ascaris.65 \
      --allow-extra-chr \
      --extract Ascaris.65.pruned.prune.in \
      --make-bed \
      --out Ascaris.65.pruned

# Randomly select up to 100,000 SNPs
shuf Ascaris.65.pruned.prune.in | head -100000 > Ascaris.65.pruned.100kSNPs.list

echo "Number of SNPs selected for AMOVA:"
wc -l Ascaris.65.pruned.100kSNPs.list

# Make final 100k SNP PLINK set
plink --bfile Ascaris.65 \
      --allow-extra-chr \
      --extract Ascaris.65.pruned.100kSNPs.list \
      --make-bed \
      --out Ascaris.65.pruned.100kSNPs

# Export bgzipped VCF
plink --bfile Ascaris.65.pruned.100kSNPs \
      --allow-extra-chr \
      --recode vcf bgz \
      --out Ascaris.65.pruned.100kSNPs

tabix -p vcf Ascaris.65.pruned.100kSNPs.vcf.gz

echo "Finished balanced 65-sample 100k SNP preparation"
date
```


## Cluster wrapper: run AMOVA R script

Purpose: run the AMOVA models on the cluster using R.

```bash
#!/bin/bash
#SBATCH --job-name=amova65_balanced
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=24:00:00
#SBATCH --output=amova65_balanced.%j.out
#SBATCH --error=amova65_balanced.%j.err

set -eo pipefail

module load GCC/13.3.0
module load R/4.4.2

Rscript /parallel_scratch/lw01165/05_analysis/05i_pixy/amova_65_balanced_all_models.R
```


## RStudio/cluster R code: AMOVA models

Purpose: run the balanced n=65 AMOVA sensitivity models: region only, region/country, host only, and region/host. The full n=146 models are preserved separately in Section 08B. 


```r
# =============================
# AMOVA on balanced 65-sample dataset
# Models:
#   1. region only
#   2. region / country
#   3. host only
#   4. region / host
#
# Fixes included:
#   - Cleans duplicated VCF sample names such as:
#       AL17SLB_merge_AL17SLB_merge -> AL17SLB_merge
#   - Removes name attributes from sample vectors before identical() checks
#   - Writes UPDATED_ outputs
# =============================

.libPaths("~/R/library")

suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(stringr)
  library(adegenet)
  library(poppr)
})

# -----------------------------
# Helper: collapse duplicated sample IDs
# -----------------------------
collapse_duplicated_sample_id <- function(x) {
  out <- vapply(x, function(id) {
    id <- str_trim(id)

    parts <- strsplit(id, "_", fixed = TRUE)[[1]]

    if (length(parts) >= 2) {
      for (i in seq_len(length(parts) - 1)) {
        left  <- paste(parts[1:i], collapse = "_")
        right <- paste(parts[(i + 1):length(parts)], collapse = "_")

        if (identical(left, right)) {
          return(left)
        }
      }
    }

    return(id)
  }, character(1), USE.NAMES = FALSE)

  return(unname(out))
}

# -----------------------------
# Helper: safe identical check
# -----------------------------
same_ids <- function(x, y) {
  identical(
    unname(as.character(x)),
    unname(as.character(y))
  )
}

# -----------------------------
# File paths
# -----------------------------
vcf_file <- "/parallel_scratch/lw01165/05_analysis/05i_pixy/amova_65_100k_dataset/Ascaris.65.pruned.100kSNPs.vcf.gz"
meta_file <- "/parallel_scratch/lw01165/05_analysis/sample_QC_list.csv"
outdir   <- "/parallel_scratch/lw01165/05_analysis/05i_pixy/amova_results_65_balanced"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("======================================\n")
cat("UPDATED balanced 65-sample AMOVA started\n")
cat("Time: ", as.character(Sys.time()), "\n", sep = "")
cat("VCF: ", vcf_file, "\n", sep = "")
cat("Metadata: ", meta_file, "\n", sep = "")
cat("Output directory: ", outdir, "\n", sep = "")
cat("All output files will begin with: UPDATED_\n")
cat("======================================\n\n")

# -----------------------------
# Read data
# -----------------------------
cat("Reading VCF...\n")
vcf <- read.vcfR(vcf_file)

cat("Reading metadata...\n")
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# -----------------------------
# Clean metadata
# -----------------------------
required_cols <- c("sample_id", "region", "country", "host")
missing_cols <- setdiff(required_cols, names(meta))

if (length(missing_cols) > 0) {
  stop(
    "ERROR: metadata file is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

meta <- meta %>%
  rename(IID = sample_id) %>%
  mutate(
    IID     = str_trim(IID),
    region  = str_to_lower(str_trim(region)),
    country = str_to_title(str_to_lower(str_trim(country))),
    host    = str_to_lower(str_trim(host))
  ) %>%
  mutate(
    region  = na_if(region, ""),
    country = na_if(country, ""),
    host    = na_if(host, "")
  )

# -----------------------------
# Clean VCF sample names
# -----------------------------
cat("Cleaning VCF sample names...\n")

vcf_samples_raw <- colnames(vcf@gt)[-1]
vcf_samples_raw <- unname(str_trim(vcf_samples_raw))

vcf_samples <- collapse_duplicated_sample_id(vcf_samples_raw)
vcf_samples <- unname(as.character(vcf_samples))

sample_name_map <- data.frame(
  raw_vcf_sample     = vcf_samples_raw,
  cleaned_vcf_sample = vcf_samples,
  changed            = vcf_samples_raw != vcf_samples,
  stringsAsFactors   = FALSE
)

write.csv(
  sample_name_map,
  file = file.path(outdir, "UPDATED_vcf_sample_name_cleanup_map.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("VCF sample name cleanup complete.\n")
cat("Number of VCF samples renamed: ", sum(sample_name_map$changed), "\n", sep = "")
cat("First cleaned sample names:\n")
print(head(sample_name_map, 10))

# Update VCF genotype column names to cleaned IDs
gt_names <- colnames(vcf@gt)
gt_names[-1] <- vcf_samples
colnames(vcf@gt) <- gt_names

# -----------------------------
# Check sample IDs
# -----------------------------
cat("\nChecking sample mismatches after VCF name cleanup...\n")

if (any(is.na(vcf_samples) | vcf_samples == "")) {
  stop("ERROR: VCF contains missing or blank sample names after cleaning.")
}

if (any(duplicated(vcf_samples))) {
  cat("ERROR: duplicated sample IDs in VCF after cleaning:\n")
  print(vcf_samples[duplicated(vcf_samples)])
  stop("VCF contains duplicated sample IDs after cleaning.")
}

bad_meta_iid <- is.na(meta$IID) | meta$IID == ""

if (any(bad_meta_iid)) {
  cat("ERROR: metadata contains missing or blank sample IDs:\n")
  print(meta[bad_meta_iid, ])
  stop("Metadata contains missing or blank sample IDs.")
}

meta_overlap <- meta %>% filter(IID %in% vcf_samples)

if (any(duplicated(meta_overlap$IID))) {
  cat("ERROR: duplicated metadata sample IDs among VCF samples:\n")
  print(meta_overlap$IID[duplicated(meta_overlap$IID)])
  stop("Metadata contains duplicated sample IDs for samples present in the VCF.")
}

vcf_not_meta <- setdiff(vcf_samples, meta$IID)
meta_not_vcf <- setdiff(meta$IID, vcf_samples)

cat("Samples in VCF but not metadata: ", length(vcf_not_meta), "\n", sep = "")
if (length(vcf_not_meta) > 0) print(vcf_not_meta)

cat("Samples in metadata but not VCF: ", length(meta_not_vcf), "\n", sep = "")
if (length(meta_not_vcf) > 0) print(meta_not_vcf)

writeLines(
  c(
    "UPDATED sample matching report",
    paste0("Time: ", as.character(Sys.time())),
    paste0("VCF samples after cleanup: ", length(vcf_samples)),
    paste0("Metadata rows: ", nrow(meta)),
    paste0("Samples in VCF but not metadata: ", length(vcf_not_meta)),
    paste0("Samples in metadata but not VCF: ", length(meta_not_vcf)),
    paste0("VCF sample names renamed: ", sum(sample_name_map$changed))
  ),
  con = file.path(outdir, "UPDATED_sample_matching_report.txt")
)

write.csv(
  data.frame(sample_id = vcf_not_meta),
  file = file.path(outdir, "UPDATED_samples_in_vcf_not_metadata.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(sample_id = meta_not_vcf),
  file = file.path(outdir, "UPDATED_samples_in_metadata_not_vcf.csv"),
  row.names = FALSE,
  quote = FALSE
)

if (length(vcf_not_meta) > 0) {
  stop("ERROR: some cleaned VCF samples are still missing from metadata. See UPDATED_samples_in_vcf_not_metadata.csv")
}

# -----------------------------
# Reorder metadata to match VCF order
# -----------------------------
match_index <- match(vcf_samples, meta$IID)

if (any(is.na(match_index))) {
  cat("ERROR: unmatched VCF samples after match():\n")
  print(vcf_samples[is.na(match_index)])
  stop("Metadata matching failed.")
}

meta <- meta[match_index, , drop = FALSE]
rownames(meta) <- NULL

if (!same_ids(meta$IID, vcf_samples)) {
  cat("First metadata IDs:\n")
  print(head(meta$IID, 20))

  cat("First VCF sample IDs:\n")
  print(head(vcf_samples, 20))

  mismatch_rows <- which(as.character(meta$IID) != as.character(vcf_samples))

  cat("\nRows where metadata and VCF sample IDs differ:\n")
  print(head(data.frame(
    row = mismatch_rows,
    metadata_IID = meta$IID[mismatch_rows],
    vcf_sample = vcf_samples[mismatch_rows]
  ), 30))

  stop("ERROR: metadata and VCF sample order do not match after reordering.")
}

write.csv(
  meta,
  file = file.path(outdir, "UPDATED_clean_metadata_matched_to_vcf.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("\nMetadata aligned to cleaned VCF successfully.\n")
cat("Total matched VCF samples: ", nrow(meta), "\n\n", sep = "")

cat("Sample counts by region:\n")
print(table(meta$region, useNA = "ifany"))

cat("\nSample counts by country:\n")
print(table(meta$country, useNA = "ifany"))

cat("\nSample counts by host:\n")
print(table(meta$host, useNA = "ifany"))

cat("\nCross-tab: region x host\n")
print(table(meta$region, meta$host, useNA = "ifany"))

# -----------------------------
# Convert VCF to genind once
# -----------------------------
cat("\nConverting VCF to genind...\n")
genind_obj <- vcfR2genind(vcf)

rm(vcf)
gc()

indNames(genind_obj) <- unname(str_trim(indNames(genind_obj)))

if (!same_ids(indNames(genind_obj), meta$IID)) {
  cat("First genind IDs:\n")
  print(head(indNames(genind_obj), 20))

  cat("First metadata IDs:\n")
  print(head(meta$IID, 20))

  mismatch_rows <- which(as.character(indNames(genind_obj)) != as.character(meta$IID))

  cat("\nRows where genind and metadata IDs differ:\n")
  print(head(data.frame(
    row = mismatch_rows,
    genind_ID = indNames(genind_obj)[mismatch_rows],
    metadata_IID = meta$IID[mismatch_rows]
  ), 30))

  stop("ERROR: genind individual names do not match metadata IID order.")
}

# ============================================================
# Helper function to run and save AMOVA models
# ============================================================

run_amova_model <- function(model_label,
                            keep_samples,
                            meta_used,
                            strata_data,
                            formula_text,
                            metadata_filename,
                            output_stub) {

  cat("\n======================================\n")
  cat("Running UPDATED AMOVA: ", model_label, "\n", sep = "")
  cat("======================================\n")

  keep_samples <- as.logical(keep_samples)

  if (length(keep_samples) != nrow(meta)) {
    stop(paste0("ERROR: keep_samples length does not match metadata for ", model_label))
  }

  if (nrow(meta_used) < 2) {
    stop(paste0("ERROR: fewer than two samples retained for ", model_label))
  }

  if (nrow(strata_data) != nrow(meta_used)) {
    stop(paste0("ERROR: strata_data row count does not match metadata for ", model_label))
  }

  rownames(strata_data) <- meta_used$IID

  cat("Samples retained: ", nrow(meta_used), "\n", sep = "")

  write.csv(
    meta_used,
    file = file.path(outdir, paste0("UPDATED_", metadata_filename)),
    row.names = FALSE,
    quote = FALSE
  )

  genind_subset <- genind_obj[keep_samples]

  if (!same_ids(indNames(genind_subset), meta_used$IID)) {
    cat("First subset genind IDs:\n")
    print(head(indNames(genind_subset), 20))

    cat("First metadata IDs:\n")
    print(head(meta_used$IID, 20))

    mismatch_rows <- which(as.character(indNames(genind_subset)) != as.character(meta_used$IID))

    cat("\nRows where subset genind and metadata IDs differ:\n")
    print(head(data.frame(
      row = mismatch_rows,
      genind_ID = indNames(genind_subset)[mismatch_rows],
      metadata_IID = meta_used$IID[mismatch_rows]
    ), 30))

    stop(paste0("ERROR: ", model_label, " genind names do not match metadata."))
  }

  genclone_subset <- as.genclone(genind_subset)

  strata(genclone_subset) <- strata_data

  amova_result <- poppr.amova(
    genclone_subset,
    as.formula(formula_text)
  )

  print(amova_result)

  cat("\nRunning permutation test (999 replicates)...\n")
  amova_test <- randtest(amova_result, nrepet = 999)

  print(amova_test)

  write.table(
    capture.output(amova_result),
    file = file.path(outdir, paste0("UPDATED_", output_stub, ".txt")),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  write.table(
    capture.output(amova_test),
    file = file.path(outdir, paste0("UPDATED_", output_stub, "_randtest.txt")),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  saveRDS(
    amova_result,
    file = file.path(outdir, paste0("UPDATED_", output_stub, ".rds"))
  )

  saveRDS(
    amova_test,
    file = file.path(outdir, paste0("UPDATED_", output_stub, "_randtest.rds"))
  )

  cat("\nUPDATED output files written for: ", model_label, "\n", sep = "")
  cat("Text file: UPDATED_", output_stub, ".txt\n", sep = "")
  cat("Randtest file: UPDATED_", output_stub, "_randtest.txt\n", sep = "")

  rm(genind_subset, genclone_subset, amova_result, amova_test)
  gc()
}

# ============================================================
# 1. AMOVA: region only
# ============================================================

keep_region <- !is.na(meta$region)
meta_region <- meta[keep_region, , drop = FALSE]

cat("\nRegion-only model sample counts:\n")
print(table(meta_region$region, useNA = "ifany"))

run_amova_model(
  model_label       = "region only",
  keep_samples      = keep_region,
  meta_used         = meta_region,
  strata_data       = data.frame(region = factor(meta_region$region)),
  formula_text      = "~region",
  metadata_filename = "metadata_region_only_used.csv",
  output_stub       = "amova_region_only"
)

# ============================================================
# 2. AMOVA: region / country
# ============================================================

keep_region_country <- !is.na(meta$region) & !is.na(meta$country)
meta_region_country <- meta[keep_region_country, , drop = FALSE]

cat("\nRegion / country model sample counts:\n")
print(table(meta_region_country$region, meta_region_country$country, useNA = "ifany"))

run_amova_model(
  model_label       = "region / country",
  keep_samples      = keep_region_country,
  meta_used         = meta_region_country,
  strata_data       = data.frame(
    region  = factor(meta_region_country$region),
    country = factor(meta_region_country$country)
  ),
  formula_text      = "~region/country",
  metadata_filename = "metadata_region_country_used.csv",
  output_stub       = "amova_region_country"
)

# ============================================================
# 3. AMOVA: host only
# ============================================================

keep_host <- !is.na(meta$host)
meta_host <- meta[keep_host, , drop = FALSE]

cat("\nHost-only model sample counts:\n")
print(table(meta_host$host, useNA = "ifany"))

run_amova_model(
  model_label       = "host only",
  keep_samples      = keep_host,
  meta_used         = meta_host,
  strata_data       = data.frame(host = factor(meta_host$host)),
  formula_text      = "~host",
  metadata_filename = "metadata_host_used.csv",
  output_stub       = "amova_host"
)

# ============================================================
# 4. AMOVA: region / host
# ============================================================

keep_region_host <- !is.na(meta$region) & !is.na(meta$host)
meta_region_host <- meta[keep_region_host, , drop = FALSE]

cat("\nRegion / host model sample counts:\n")
print(table(meta_region_host$region, meta_region_host$host, useNA = "ifany"))

run_amova_model(
  model_label       = "region / host",
  keep_samples      = keep_region_host,
  meta_used         = meta_region_host,
  strata_data       = data.frame(
    region = factor(meta_region_host$region),
    host   = factor(meta_region_host$host)
  ),
  formula_text      = "~region/host",
  metadata_filename = "metadata_region_host_used.csv",
  output_stub       = "amova_region_host"
)

# -----------------------------
# Finish
# -----------------------------
rm(genind_obj)
gc()

cat("\n======================================\n")
cat("UPDATED balanced 65-sample AMOVA complete\n")
cat("Time: ", as.character(Sys.time()), "\n", sep = "")
cat("Results written to: ", outdir, "\n", sep = "")
cat("UPDATED files created/overwritten:\n")
print(list.files(outdir, pattern = "^UPDATED_", full.names = FALSE))
cat("======================================\n")
```

---

## Outputs and QC checks

Wahlund H_T / H_S tests, AMOVA models with permutation tests. Consult recorded logs and intermediate tables to distinguish expected settings from actual execution success.

## Provenance and limitations

- Original full n=146 R scripts are now archived verbatim for five models (region only, host only, region/country, region/host and within-Africa host), alongside the full 100,000-SNP preparation job. The source-path discrepancy in `amova_region.R` and metadata-header compatibility require verification against historical inputs.
- Sample identifiers and host/country metadata were manually curated in `sample_QC_list.csv`; it is not generated by an algorithm in these files.
- Commands, versions, seeds and paths are retained where present; annotations are not evidence that all jobs have been freshly rerun.
- For the discovery and correction of the PIXY all-sites input issue, see [`99_PIXY_allsites_correction.md`](99_PIXY_allsites_correction.md).
