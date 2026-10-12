# 06 — Population structure and balanced-sample sensitivity analysis

## Purpose

Population structure and balanced-sample sensitivity analysis. This file archives the commands and analysis notes from the thesis/paper workflow.

## Inputs

Final n=146 SNP VCF, curated metadata, balanced-sample selection lists. Source filenames and absolute HPC paths are retained in the original script blocks below.

## Workflow

See the original annotated analysis sequence below. Distinct script blocks are independent jobs/R sessions, **not one executable script**.

## Annotated scripts

### 06A. Balanced sample-list construction — recovered historical selection logic

The fixed SRR3 and SRR8 lists plus deterministic `shuf` draws were used to create `balanced_samples.txt`. This is distinct from the later VCF intersection. **Required fixed lists and the final 65 IDs are inputs that should be archived as files, not regenerated from assumptions.**

```bash
#!/usr/bin/env bash
set -euo pipefail

# Prepare fixed SRR3 and SRR8 sample lists as used in the historical analysis.
# all_samples.txt contains sample IDs for the final QC set.
# The required fixed lists must be retained with this repository to recreate the exact draw.
# bcftools query -l Ascaris.autosomes.SNPs.146.maxmiss0.90.vcf.gz > all_samples.txt

grep "^SRR3" all_samples.txt > all_SRR3.txt
grep "^SRR8" all_samples.txt > all_SRR8.txt
grep -v -E "^SRR3|^SRR8" all_samples.txt > non_SRR38_samples.txt
grep -v -F -f keep_SRR3_fixed.txt all_SRR3.txt > SRR3_pool.txt
grep -v -F -f keep_SRR8_fixed.txt all_SRR8.txt > SRR8_pool.txt
shuf --random-source=<(yes 42) -n 10 SRR3_pool.txt > SRR3_random10.txt
shuf --random-source=<(yes 84) -n 15 SRR8_pool.txt > SRR8_random15.txt
cat keep_SRR3_fixed.txt SRR3_random10.txt > final_SRR3_samples.txt
cat keep_SRR8_fixed.txt SRR8_random15.txt > final_SRR8_samples.txt
cat non_SRR38_samples.txt final_SRR3_samples.txt final_SRR8_samples.txt > balanced_samples.txt
wc -l balanced_samples.txt
sort balanced_samples.txt | uniq -d || true

```

### 06B. Construct balanced 65-sample VCF — user-supplied source

The complete historical job was supplied as a message. For readability this copy omits some `echo` status-only statements but preserves all substantive commands and parameters; retain the original `.sub` separately if exact byte-level replication is important.

```bash
#!/bin/bash
#SBATCH --job-name=bal146
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --output=bal146.%j.out
#SBATCH --error=bal146.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env

VCF=/parallel_scratch/lw01165/05_analysis/05b_146_data/Ascaris.autosomes.SNPs.146.maxmiss0.90.vcf.gz
BAL=/parallel_scratch/lw01165/05_analysis/05a_Balanced_dataset/balanced_samples.txt
OUTDIR=/parallel_scratch/lw01165/05_analysis/05a_Balanced_dataset

OUTVCF=${OUTDIR}/Ascaris.autosomes.SNPs.146.maxmiss0.90.balanced.n65.vcf.gz
KEEP=${OUTDIR}/balanced_samples.in146.txt
STATS=${OUTDIR}/Ascaris.autosomes.SNPs.146.maxmiss0.90.balanced.n65.stats.tsv

TMPVCF=${OUTVCF%.vcf.gz}.tmp.vcf.gz
TMPSTATS=${STATS%.tsv}.tmp.tsv

mkdir -p "$OUTDIR"

if [ ! -f "${VCF}.tbi" ]; then
  bcftools index -t "${VCF}"
fi

comm -12 <(bcftools query -l "$VCF" | sort) <(sed '/^$/d' "$BAL" | sort -u) > "$KEEP"

n_keep=$(wc -l < "$KEEP")
echo "n_keep_intersection=${n_keep}"
EXPECTED=65
if [ "$n_keep" -ne "$EXPECTED" ]; then
  echo "ERROR: Expected $EXPECTED balanced samples, but intersection has $n_keep."
  comm -23 <(sed '/^$/d' "$BAL" | sort -u) <(bcftools query -l "$VCF" | sort) | head -n 50
  exit 1
fi

bcftools view -S "$KEEP" -Oz -o "$TMPVCF" "$VCF"
bcftools index -t "$TMPVCF"

echo -e "metric\tvalue" > "$TMPSTATS"
echo -e "n_samples_balanced\t$(bcftools query -l $TMPVCF | wc -l)" >> "$TMPSTATS"
echo -e "n_variants\t$(bcftools view -H $TMPVCF | wc -l)" >> "$TMPSTATS"

mv -f "$TMPVCF" "$OUTVCF"
mv -f "$TMPVCF.tbi" "$OUTVCF.tbi"
mv -f "$TMPSTATS" "$STATS"

echo "Done. Final outputs:"
ls -lh "$OUTVCF" "${OUTVCF}.tbi" "$KEEP" "$STATS"
echo -n "Final n_samples in balanced VCF: "
bcftools query -l "$OUTVCF" | wc -l

```

### 06C. Population-structure analysis transcript

**Source-format issue:** The historic source was exported with escaped Markdown characters (`\#`, `\_`, `&#x20;`) and contains multiple independent Bash and R blocks. The escapes have been decoded, and redundant blank lines reduced, for readable display. Each section remains an **archival source transcript**, not a single executable Bash/R script. Extract individual jobs by their shebangs (`#!/bin/bash`) or R package starts before execution.

#### PLINK LD pruning and PCA

```text
# ==================================================================================================

# Annotated population-structure code notes

# --------------------------------------------------------------------------------------------------

#

# Main workflow:

#   1) Convert filtered VCFs to PLINK binary format and LD-prune SNPs.

#   2) Run PCA on the LD-pruned full and balanced datasets, including subset-specific PCAs.

#   3) Calculate pairwise 1-IBS genetic distance matrices from LD-pruned SNP data.

#   4) Use the 1-IBS genetic distance matrices for PERMANOVA and dispersion tests.

#   5) Use the same 1-IBS genetic distance matrices for NJ trees and IBS similarity networks.

#   6) Build a SNP-based maximum-likelihood tree from the balanced dataset.

#   7) Run ADMIXTURE on a 100k LD-pruned SNP subset and visualise K = 2–4.

#

# ==================================================================================================

5.2.3 Population structure Analysis

# Section purpose: analyse broad population structure using complementary methods: PCA, PERMANOVA,

# tree-based clustering, IBS similarity, maximum-likelihood tree building, and ADMIXTURE.

I) This script takes the full n = 146 dataset and the balanced n = 65 dataset  and converts them into PLINK format, performs LD pruning, and then creates new LD-pruned PLINK datasets.

# Purpose: create genotype datasets suitable for PCA, NJ trees, IBS distances, and ADMIXTURE.

# Input: autosomal SNP VCFs for the full dataset and the balanced downsampled dataset.

# Output: PLINK .bed/.bim/.fam files before and after LD pruning.

# Rationale: LD pruning reduces the influence of clusters of highly correlated nearby SNPs, so

# downstream structure analyses are not dominated by dense linked regions.

#!/bin/bash

#SBATCH --job-name=plink_prune_146_65

#SBATCH --partition=high_mem

#SBATCH --cpus-per-task=16

#SBATCH --mem=120G

#SBATCH --time=12:00:00

#SBATCH --output=plink_prune_146_65.%j.out

#SBATCH --error=plink_prune_146_65.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh

conda activate NJ_env

# -----------------------------

# Input files

# -----------------------------

VCF146=/parallel_scratch/lw01165/05_analysis/05b_146_data/Ascaris.autosomes.SNPs.146.maxmiss0.90.vcf.gz

VCF65=/parallel_scratch/lw01165/05_analysis/05a_Balanced_dataset/Ascaris.autosomes.SNPs.146.maxmiss0.90.balanced.n65.vcf.gz

CHROM_FILE=/parallel_scratch/lw01165/chrom_order.txt

# chrom_order.txt defines the autosomes to retain and preserves the intended chromosome order.

# -----------------------------

# Output directories

# -----------------------------

OUT146=/parallel_scratch/lw01165/05_analysis/05d_plink_146

OUT65=/parallel_scratch/lw01165/05_analysis/05e_plink_65

mkdir -p ${OUT146}

mkdir -p ${OUT65}

# Convert chromosome list to comma-separated string for PLINK.

# PLINK expects --chr to receive a comma-separated list, so paste -sd, collapses one chromosome per

# line into one comma-delimited string.

CHR_LIST=$(paste -sd, ${CHROM_FILE})

echo "Chromosomes used:"

echo "${CHR_LIST}"

echo "Starting PLINK conversion and LD pruning..."

date

# -----------------------------

# 146-sample dataset

# -----------------------------

# Convert the full VCF to PLINK binary format. Variant IDs are reset to CHROM:POS so that variants

# have consistent IDs for later extraction and pruning.

# NOTE: if duplicate variants exist at the same coordinate, CHROM:POS IDs can become duplicated.

echo "Converting 146-sample VCF to PLINK..."

plink2 \
  --vcf ${VCF146} \
  --chr ${CHR_LIST} \
  --make-bed \
  --allow-extra-chr \
  --set-all-var-ids @:# \
  --out ${OUT146}/Ascaris.146

# LD pruning for the full dataset. Parameters are window = 50 variants, step = 10 variants, r^2 = 0.15.

# SNP pairs with r^2 above 0.15 within the sliding window are treated as highly correlated; one SNP is

# removed from correlated pairs/sets.

echo "LD pruning 146-sample dataset..."

plink2 \
  --bfile ${OUT146}/Ascaris.146 \
  --allow-extra-chr \
  --indep-pairwise 50 10 0.15 \
  --out ${OUT146}/Ascaris.146

# Create the final full LD-pruned PLINK dataset by extracting the SNPs retained in prune.in.

echo "Making pruned 146-sample PLINK dataset..."

plink2 \
  --bfile ${OUT146}/Ascaris.146 \
  --allow-extra-chr \
  --extract ${OUT146}/Ascaris.146.prune.in \
  --make-bed \
  --out ${OUT146}/Ascaris.146.pruned

echo "Calculating allele frequencies for pruned 146-sample dataset..."

plink2 \
  --bfile ${OUT146}/Ascaris.146.pruned \
  --allow-extra-chr \
  --freq \
  --out ${OUT146}/Ascaris.146.pruned

# -----------------------------

# 65-sample dataset

# -----------------------------

# Convert the balanced VCF to PLINK binary format using the same chromosome list and variant ID scheme.

echo "Converting 65-sample VCF to PLINK..."

plink2 \
  --vcf ${VCF65} \
  --chr ${CHR_LIST} \
  --make-bed \
  --allow-extra-chr \
  --set-all-var-ids @:# \
  --out ${OUT65}/Ascaris.65

# LD pruning for the balanced dataset using the same pruning settings as the full dataset.

echo "LD pruning 65-sample dataset..."

plink2 \
  --bfile ${OUT65}/Ascaris.65 \
  --allow-extra-chr \
  --indep-pairwise 50 10 0.15 \
  --out ${OUT65}/Ascaris.65

# Create the final balanced LD-pruned PLINK dataset by extracting the SNPs retained in prune.in.

echo "Making pruned 65-sample PLINK dataset..."

plink2 \
  --bfile ${OUT65}/Ascaris.65 \
  --allow-extra-chr \
  --extract ${OUT65}/Ascaris.65.prune.in \
  --make-bed \
  --out ${OUT65}/Ascaris.65.pruned

echo "Finished PLINK conversion and pruning."

date

IIa) Principal Component Analysis in Cluster

# Purpose: run PLINK PCA on the LD-pruned full and balanced datasets.

# Output: .eigenvec files containing sample PC coordinates and .eigenval files containing eigenvalues.

#!/bin/bash

#SBATCH --job-name=pca_146_65

#SBATCH --partition=high_mem

#SBATCH --cpus-per-task=16

#SBATCH --mem=80G

#SBATCH --time=04:00:00

#SBATCH --output=pca_146_65.%j.out

#SBATCH --error=pca_146_65.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh

conda activate NJ_env

IN146=/parallel_scratch/lw01165/05_analysis/05d_plink_146/Ascaris.146.pruned

IN65=/parallel_scratch/lw01165/05_analysis/05e_plink_65/Ascaris.65.pruned

OUT146=/parallel_scratch/lw01165/05_analysis/05f_pca_146

OUT65=/parallel_scratch/lw01165/05_analysis/05g_pca_65

mkdir -p ${OUT146}

mkdir -p ${OUT65}

# Run PCA on the LD-pruned full dataset. PCA summarises major axes of genotype variation.

echo "Running PCA on 146-sample pruned dataset..."

plink2 \
  --bfile ${IN146} \
  --allow-extra-chr \
  --pca \
  --out ${OUT146}/Ascaris.146.pca

# Run PCA on the LD-pruned balanced dataset, useful for checking whether the full-dataset structure is

# driven by sample-size imbalance.

echo "Running PCA on 65-sample pruned dataset..."

plink2 \
  --bfile ${IN65} \
  --allow-extra-chr \
  --pca \
  --out ${OUT65}/Ascaris.65.pca

echo "PCA complete."

date

IIb) Principal Component Analysis using RStudio. The following code was also repeated for the full dataset. 

# Purpose: merge PCA coordinates with metadata, clean plotting factors, calculate variance explained,

# and save publication-style PCA plots.

# Note: PCA axes can flip sign between runs/datasets; sign flips do not change the biological meaning.

setwd("S:/Research/Lauren_PopGen/PCA/65")

# ----------------------------

# Load libraries

# ----------------------------

library(dplyr)

library(stringr)

library(ggplot2)

library(patchwork)

# ----------------------------

# Read metadata and PCA

# ----------------------------

samples <- read.delim("sample_146_list.txt", header = TRUE, stringsAsFactors = FALSE)

# Metadata table. For the balanced PCA, inner_join below keeps only samples present in the 65-sample

# PCA output.

pca <- read.table("Ascaris.65.pca.eigenvec", header = FALSE)

eig <- read.table("Ascaris.65.pca.eigenval", header = FALSE)

colnames(pca) <- c("FID", "IID", paste0("PC", 1:(ncol(pca) - 2)))

# ----------------------------

# Merge metadata + PCA coordinates

# Keep only samples present in the balanced PCA

# ----------------------------

pca_meta <- inner_join(samples, pca, by = c("sample_id" = "IID"))

# inner_join is appropriate here because the metadata file may contain all 146 samples, while the PCA

# file contains only the samples included in the current PLINK run.

# ----------------------------

# Align balanced PCA orientation to full dataset

# ----------------------------

pca_meta <- pca_meta %>%

  mutate(

    PC2 = -PC2

  )

# PC2 is flipped so the balanced PCA is visually aligned with the full dataset. PCA sign is arbitrary,

# so this does not alter distances or clustering.

# ----------------------------

# Clean metadata

# ----------------------------

pca_meta <- pca_meta %>%

  mutate(

    region  = str_trim(region),

    country = str_trim(country),

    host    = str_trim(host),

    region  = ifelse(is.na(region)  | region == "",  "Unknown", str_to_title(region)),

    country = ifelse(is.na(country) | country == "", "Unknown", str_to_title(country)),

    host    = ifelse(is.na(host)    | host == "",    "Unknown", str_to_title(host))

  )

pca_meta$region <- factor(

  pca_meta$region,

  levels = c("Africa", "Asia", "Europe", "Unknown")

)

pca_meta$country <- factor(

  pca_meta$country,

  levels = c("Ethiopia", "Germany", "Hungary", "India", "Ireland",

             "Kenya", "Philippines", "Thailand", "Uganda", "Unknown")

)

pca_meta$host <- factor(

  pca_meta$host,

  levels = c("Human", "Pig", "Unknown")

)

# ----------------------------

# Variance explained

# ----------------------------

var_exp <- eig$V1 / sum(eig$V1) * 100

# Convert eigenvalues into the percentage of total PCA variance explained by each PC.

pc1_var <- round(var_exp[1], 1)

pc2_var <- round(var_exp[2], 1)

pc3_var <- round(var_exp[3], 1)

pc4_var <- round(var_exp[4], 1)

# ----------------------------

# Shared theme

# ----------------------------

theme_pca <- theme_classic(base_size = 15) +

  theme(

    legend.title = element_text(face = "bold", size = 13),

    legend.text  = element_text(size = 11),

    legend.position = "right",

    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),

    axis.title = element_text(size = 15, face = "bold"),

    axis.text  = element_text(size = 12)

  )

# ----------------------------

# Helper for axes

# ----------------------------

# Helper functions keep axis limits consistent while adding a small amount of white space around the

# plotted points.

axis_pad <- function(x, pad = 0.025) {

  r <- range(x, na.rm = TRUE)

  w <- diff(r)

  c(r[1] - w * pad, r[2] + w * pad)

}

coord_pca <- function(data, xvar, yvar, pad_x = 0.025, pad_y = 0.035) {

  coord_cartesian(

    xlim = axis_pad(data[[xvar]], pad_x),

    ylim = axis_pad(data[[yvar]], pad_y)

  )

}

# ----------------------------

# Shared colours and shapes

# ----------------------------

# Manual colour/shape definitions ensure that the same metadata groups are displayed consistently

# across figures.

region_cols <- c(

  Africa  = "#0047AB",

  Asia    = "#F8766D",

  Europe  = "#B79F00",

  Unknown = "grey50"

)

country_cols <- c(

  Ethiopia    = "#00BFC4",

  Germany     = "#E69F00",

  Hungary     = "gold",

  India       = "#4CAF50",

  Ireland     = "green",

  Kenya       = "#0047AB",

  Philippines = "#F8766D",

  Thailand    = "#8B80F9",

  Uganda      = "#F564E3",

  Unknown     = "grey50"

)

host_cols <- c(

  Human   = "#8702b0",

  Pig     = "darkgreen",

  Unknown = "grey50"

)

host_shapes <- c(

  Human   = 21,

  Pig     = 24,

  Unknown = 22

)

# ========================================================

# MAIN BALANCED PCA FIGURE

# A. Country

# B. Host

# C. Region + host

# D. PC1 vs PC3 by region + host

# ========================================================

p_country <- ggplot(pca_meta, aes(PC1, PC2, fill = country)) +

  geom_point(

    shape = 21,

    size = 2.8,

    alpha = 0.85,

    colour = "grey10",

    stroke = 0.25

  ) +

  scale_fill_manual(

    values = country_cols,

    drop = FALSE,

    guide = guide_legend(override.aes = list(shape = 21))

  ) +

  coord_pca(pca_meta, "PC1", "PC2") +

  labs(

    x = paste0("PC1 (", pc1_var, "%)"),

    y = paste0("PC2 (", pc2_var, "%)"),

    fill = "Country",

    title = "A. By country"

  ) +

  theme_pca

p_host <- ggplot(pca_meta, aes(PC1, PC2, fill = host, shape = host)) +

  geom_point(

    size = 2.8,

    alpha = 0.85,

    colour = "grey10",

    stroke = 0.25

  ) +

  scale_fill_manual(values = host_cols, drop = FALSE) +

  scale_shape_manual(values = host_shapes, drop = FALSE) +

  coord_pca(pca_meta, "PC1", "PC2") +

  labs(

    x = paste0("PC1 (", pc1_var, "%)"),

    y = paste0("PC2 (", pc2_var, "%)"),

    fill = "Host",

    shape = "Host",

    title = "B. By host"

  ) +

  theme_pca

p_region_host <- ggplot(pca_meta, aes(PC1, PC2, fill = region, shape = host)) +

  geom_point(

    size = 2.8,

    alpha = 0.85,

    colour = "grey10",

    stroke = 0.25

  ) +

  scale_fill_manual(

    values = region_cols,

    drop = FALSE,

    guide = guide_legend(override.aes = list(shape = 21))

  ) +

  scale_shape_manual(values = host_shapes, drop = FALSE) +

  coord_pca(pca_meta, "PC1", "PC2") +

  labs(

    x = paste0("PC1 (", pc1_var, "%)"),

    y = paste0("PC2 (", pc2_var, "%)"),

    fill = "Region",

    shape = "Host",

    title = "C. By region and host"

  ) +

  theme_pca

p_pc1_pc3 <- ggplot(pca_meta, aes(PC1, PC3, fill = region, shape = host)) +

  geom_point(

    size = 2.8,

    alpha = 0.85,

    colour = "grey10",

    stroke = 0.25

  ) +

  scale_fill_manual(

    values = region_cols,

    drop = FALSE,

    guide = guide_legend(override.aes = list(shape = 21))

  ) +

  scale_shape_manual(values = host_shapes, drop = FALSE) +

  coord_pca(pca_meta, "PC1", "PC3") +

  labs(

    x = paste0("PC1 (", pc1_var, "%)"),

    y = paste0("PC3 (", pc3_var, "%)"),

    fill = "Region",

    shape = "Host",

    title = "D. PC1 vs PC3 by region and host"

  ) +

  theme_pca

final_pca_fig <- (p_country | p_host) /

  (p_region_host | p_pc1_pc3)

print(final_pca_fig)

ggsave(

  filename = "Ascaris_PCA_balanced_4panel.png",

  plot = final_pca_fig,

  width = 14,

  height = 11,

  dpi = 600

)

# ========================================================

# SCREE PLOT

# ========================================================

pca_var <- data.frame(

  PC = 1:length(var_exp),

  Variance = var_exp

)

p_scree <- ggplot(pca_var[1:10, ], aes(x = PC, y = Variance)) +

  geom_point(size = 3) +

  geom_line() +

  scale_x_continuous(breaks = 1:10) +

  labs(

    x = "Principal component",

    y = "Variance explained (%)",

    title = "Balanced dataset scree plot"

  ) +

  theme_pca +

  theme(legend.position = "none")

print(p_scree)

ggsave(

  filename = "Ascaris_PCA_balanced_scree_plot.png",

  plot = p_scree,

  width = 7,

  height = 5,

  dpi = 600

)

# ========================================================

# PHILIPPINES-EXCLUDED VISUALISATION

# Purpose: replot the same PCA coordinates after removing Philippines samples, so that structure among

# the remaining samples is easier to see. The PCA itself is not recalculated here.

# PCA axes are still from the balanced PCA

# ========================================================

pca_meta_no_philippines <- pca_meta %>%

  filter(country != "Philippines") %>%

  droplevels()

n_no_philippines <- nrow(pca_meta_no_philippines)

# ========================================================

# Country zoom excluding Philippines

# ========================================================

p_country_zoom_no_philippines <- ggplot(

  pca_meta_no_philippines,

  aes(PC1, PC2, fill = country)

) +

  geom_point(

    shape = 21,

    size = 3,

    alpha = 0.9,

    colour = "grey10",

    stroke = 0.25

  ) +

  scale_fill_manual(

    values = country_cols,

    drop = TRUE,

    guide = guide_legend(override.aes = list(shape = 21))

  ) +

  coord_pca(

    pca_meta_no_philippines,

    "PC1",

    "PC2",

    pad_x = 0.04,

    pad_y = 0.04

  ) +

  labs(

    x = paste0("PC1 (", pc1_var, "%)"),

    y = paste0("PC2 (", pc2_var, "%)"),

    fill = "Country",

    title = paste0("Country PCA zoom excluding Philippines (n=", n_no_philippines, ")")

  ) +

  theme_pca

print(p_country_zoom_no_philippines)

ggsave(

  filename = "Ascaris_PCA_balanced_country_zoom_no_Philippines.png",

  plot = p_country_zoom_no_philippines,

  width = 8.5,

  height = 6.5,

  dpi = 600

)

# ========================================================

# PC3 vs PC4 excluding Philippines

# ========================================================

p_pc3_pc4_country_no_philippines <- ggplot(

  pca_meta_no_philippines,

  aes(PC3, PC4, fill = country, shape = host)

) +

  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +

  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +

  geom_point(

    size = 3,

    alpha = 0.9,

    colour = "grey10",

    stroke = 0.25

  ) +

  scale_fill_manual(

    values = country_cols,

    drop = TRUE,

    guide = guide_legend(override.aes = list(shape = 21))

  ) +

  scale_shape_manual(

    values = host_shapes,

    drop = FALSE

  ) +

  coord_pca(

    pca_meta_no_philippines,

    "PC3",

    "PC4",

    pad_x = 0.08,

    pad_y = 0.08

  ) +

  labs(

    x = paste0("PC3 (", pc3_var, "%)"),

    y = paste0("PC4 (", pc4_var, "%)"),

    fill = "Country",

    shape = "Host",

    title = paste0("PC3 vs PC4 by country and host, excluding Philippines (n=", n_no_philippines, ")")

  ) +

  theme_pca

print(p_pc3_pc4_country_no_philippines)

ggsave(

  filename = "Ascaris_PCA_balanced_PC3_PC4_country_host_no_Philippines.png",

  plot = p_pc3_pc4_country_no_philippines,

  width = 8.5,

  height = 6.5,

  dpi = 600

)

# ========================================================

# PC3/PC4 summary table excluding Philippines

# ========================================================

# Summarise PC3/PC4 positions by metadata group. This provides a numerical companion to the PC3/PC4

# visualisation.

pc3_pc4_summary <- pca_meta_no_philippines %>%

  group_by(region, country, host) %>%

  summarise(

    n = n(),

    mean_PC3 = mean(PC3, na.rm = TRUE),

    mean_PC4 = mean(PC4, na.rm = TRUE),

    sd_PC3 = sd(PC3, na.rm = TRUE),

    sd_PC4 = sd(PC4, na.rm = TRUE),

    .groups = "drop"

  ) %>%

  arrange(region, country, host)

print(pc3_pc4_summary)

write.table(

  pc3_pc4_summary,

  file = "Ascaris_PCA_balanced_PC3_PC4_summary_no_Philippines.tsv",

  sep = "\\t",

  quote = FALSE,

  row.names = FALSE

)

# ========================================================

# Three-panel supplementary PCA figure

# A. Country zoom excluding Philippines

# B. PC3 vs PC4 excluding Philippines

# C. Scree plot

# ========================================================

p_country_zoom_panel <- p_country_zoom_no_philippines +

  labs(title = paste0("A. Country PCA zoom excluding Philippines (n=", n_no_philippines, ")")) +

  theme(legend.position = "none")

p_pc3_pc4_country_panel <- p_pc3_pc4_country_no_philippines +

  labs(title = paste0("B. PC3 vs PC4 by country and host, excluding Philippines (n=", n_no_philippines, ")")) +

  guides(

    shape = guide_legend(order = 1),

    fill  = guide_legend(

      order = 2,

      override.aes = list(shape = 21)

    )

  )

p_scree_panel <- p_scree +

  labs(title = "C. Balanced dataset scree plot") +

  theme(legend.position = "none")

pca_extra_three_panel <- 

  (p_country_zoom_panel | p_pc3_pc4_country_panel) /

  p_scree_panel +

  plot_layout(heights = c(1.1, 0.8))

print(pca_extra_three_panel)

ggsave(

  filename = "Ascaris_PCA_balanced_zoom_PC3_PC4_scree_three_panel.png",

  plot = pca_extra_three_panel,

  width = 14,

  height = 11,

  dpi = 600

)

IIc) Principal Component Analysis on subset data on Cluster

# Purpose: rerun PCA separately for pairwise comparison subsets. This can reveal structure that may be

# compressed or hidden in the global PCA.

# Input: the full n = 146 LD-pruned PLINK dataset and popfiles listing samples to keep for each subset.

# Output: one .eigenvec/.eigenval PCA set per comparison.

#!/bin/bash

#SBATCH --job-name=subset_pca_146

#SBATCH --partition=high_mem

#SBATCH --cpus-per-task=16

#SBATCH --mem=80G

#SBATCH --time=06:00:00

#SBATCH --output=subset_pca_146.%j.out

#SBATCH --error=subset_pca_146.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh

conda activate NJ_env

BFILE=/parallel_scratch/lw01165/05_analysis/05d_plink_146/Ascaris.146.pruned

FREQ=/parallel_scratch/lw01165/05_analysis/05d_plink_146/Ascaris.146.pruned.afreq

POPDIR=/parallel_scratch/lw01165/05_analysis/05i_pixy/popfiles

OUTBASE=/parallel_scratch/lw01165/05_analysis/05k_subset_PCA

mkdir -p "${OUTBASE}"

for comp in africa_vs_asia africa_vs_europe asia_vs_europe kenya_vs_ethiopia human_vs_pig

do

    echo "Processing ${comp}..."

    COMPDIR="${OUTBASE}/${comp}"

    mkdir -p "${COMPDIR}"

    # Build a PLINK keep file from the popfile. The first column is assumed to contain sample IDs, and

    # FID is set to 0 to match the PLINK files generated from VCF.

    awk 'NR>1 {print 0, $1}' "${POPDIR}/${comp}.txt" > "${COMPDIR}/${comp}.keep"

    # Run PCA for the subset. --read-freq uses allele frequencies from the full pruned dataset.

  

    plink2 \
      --bfile "${BFILE}" \
      --keep "${COMPDIR}/${comp}.keep" \
      --read-freq "${FREQ}" \
      --pca \
      --out "${COMPDIR}/${comp}"

done

echo "All subset PCAs complete."

IId) Principal Component Analysis of sub-set data in RStudio

# Purpose: read the subset PCA outputs, merge them with metadata, and create a multi-panel figure for

# the pairwise regional/country PCA comparisons.

# ----------------------------

# Load libraries

# ----------------------------

library(dplyr)

library(stringr)

library(ggplot2)

library(patchwork)

# ----------------------------

# Read and clean metadata

# ----------------------------

meta <- read.csv("S:/Research/Lauren_PopGen/sample_QC_list.csv", stringsAsFactors = FALSE)

names(meta) <- str_trim(names(meta))

meta <- meta %>% rename(label = sample_id)

meta <- meta %>%

  mutate(

    region  = str_to_lower(str_trim(region)),

    country = str_to_title(str_trim(country)),

    host    = str_to_lower(str_trim(host))

  ) %>%

  mutate(

    host = recode(host,

                  "human" = "Human",

                  "pig"   = "Pig",

                  .default = str_to_title(host)

    ),

    region = str_to_title(region)

  )

# ----------------------------

# Function to read PCA

# ----------------------------

# Function to avoid repeating the same import/merge/variance-explained code for each subset PCA.

read_pca <- function(prefix, meta) {

  pca <- read.table(paste0(prefix, ".eigenvec"), header = FALSE)

  colnames(pca) <- c("FID", "IID", paste0("PC", 1:(ncol(pca) - 2)))

  

  eig <- read.table(paste0(prefix, ".eigenval"), header = FALSE)

  var_exp <- eig$V1 / sum(eig$V1) * 100

# Convert eigenvalues into the percentage of total PCA variance explained by each PC.

  

  x <- paste0("PC1 (", round(var_exp[1], 1), "%)")

  y <- paste0("PC2 (", round(var_exp[2], 1), "%)")

  

  pca <- pca %>%

    left_join(meta, by = c("IID" = "label"))

  

  list(

    pca = pca,

    eig = eig,

    var_exp = var_exp,

    xlab = x,

    ylab = y

  )

}

# ----------------------------

# Read PCA subsets

# ----------------------------

# Each object stores the PCA table, eigenvalues, variance explained, and axis labels.

af_as <- read_pca("S:/Research/Lauren_PopGen/PCA/Subset_PCA/Africa_vs_Asia/africa_vs_asia", meta)

af_eu <- read_pca("S:/Research/Lauren_PopGen/PCA/Subset_PCA/Africa_vs_Europe/africa_vs_europe", meta)

ke_et <- read_pca("S:/Research/Lauren_PopGen/PCA/Subset_PCA/Kenya_vs_Ethiopia/kenya_vs_ethiopia", meta)

asia_eu <- read_pca("S:/Research/Lauren_PopGen/PCA/Subset_PCA/Asia_vs_Europe/asia_vs_europe", meta)

# ----------------------------

# Assign grouping variable

# ----------------------------

af_as$pca <- af_as$pca %>% mutate(group_colour = region)

af_eu$pca <- af_eu$pca %>% mutate(group_colour = region)

asia_eu$pca <- asia_eu$pca %>% mutate(group_colour = region)

ke_et$pca <- ke_et$pca %>% mutate(group_colour = country)

# ----------------------------

# Shared colour palette (UPDATED)

# ----------------------------

group_colours <- c(

  "Africa"   = "#0047AB",

  "Asia"     = "#D10000",

  "Europe"   = "#B79F00",

  "Kenya"    = "#0047AB",

  "Ethiopia" = "#00BFC4",

  "Unknown"  = "#00BA38"

)

group_levels <- c("Africa", "Asia", "Europe", "Kenya", "Ethiopia")

# ----------------------------

# Shapes

# ----------------------------

host_shapes <- c(

  "Human" = 21,

  "Pig"   = 24

)

# ----------------------------

# Theme

# ----------------------------

theme_pca <- theme_bw() +

  theme(

    panel.grid = element_blank(),

    axis.text = element_text(size = 11, colour = "grey10"),

    axis.title = element_text(size = 12),

    legend.title = element_text(size = 11),

    legend.text = element_text(size = 10),

    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)

  )

# ----------------------------

# Plot A

# ----------------------------

p_af_as <- ggplot(af_as$pca, aes(PC1, PC2, fill = group_colour, shape = host)) +

  geom_point(size = 2.5, alpha = 0.8, colour = "grey10", stroke = 0.25) +

  scale_fill_manual(

    values = group_colours,

    breaks = group_levels,

    drop = FALSE,

    name = "Population",

    guide = guide_legend(override.aes = list(shape = 21))

  ) +

  scale_shape_manual(values = host_shapes, name = "Host", drop = FALSE) +

  labs(title = "A. Africa vs Asia", x = af_as$xlab, y = af_as$ylab) +

  coord_cartesian(

    xlim = range(af_as$pca$PC1) + c(-0.01, 0.01),

    ylim = range(af_as$pca$PC2) + c(-0.01, 0.01)

  ) +

  theme_pca

# ----------------------------

# Plot B

# ----------------------------

p_af_eu <- ggplot(af_eu$pca, aes(PC1, PC2, fill = group_colour, shape = host)) +

  geom_point(size = 2.5, alpha = 0.6, colour = "grey10", stroke = 0.25) +

  scale_fill_manual(

    values = group_colours,

    breaks = group_levels,

    drop = FALSE,

    name = "Population",

    guide = guide_legend(override.aes = list(shape = 21))

  ) +

  scale_shape_manual(values = host_shapes, name = "Host", drop = FALSE) +

  labs(title = "B. Africa vs Europe", x = af_eu$xlab, y = af_eu$ylab) +

  coord_cartesian(

    xlim = range(af_eu$pca$PC1) + c(-0.01, 0.01),

    ylim = range(af_eu$pca$PC2) + c(-0.01, 0.01)

  ) +

  theme_pca

# ----------------------------

# Plot C

# ----------------------------

p_asia_eu <- ggplot(asia_eu$pca, aes(PC1, PC2, fill = group_colour, shape = host)) +

  geom_point(size = 2.5, alpha = 0.6, colour = "grey10", stroke = 0.25) +

  scale_fill_manual(

    values = group_colours,

    breaks = group_levels,

    drop = FALSE,

    name = "Population",

    guide = guide_legend(override.aes = list(shape = 21))

  ) +

  scale_shape_manual(values = host_shapes, name = "Host", drop = FALSE) +

  labs(title = "C. Asia vs Europe", x = asia_eu$xlab, y = asia_eu$ylab) +

  coord_cartesian(

    xlim = range(asia_eu$pca$PC1) + c(-0.01, 0.01),

    ylim = range(asia_eu$pca$PC2) + c(-0.01, 0.01)

  ) +

  theme_pca

# ----------------------------

# Plot D

# ----------------------------

p_ke_et <- ggplot(ke_et$pca, aes(PC1, PC2, fill = group_colour, shape = host)) +

  geom_point(size = 2.5, alpha = 0.6, colour = "grey10", stroke = 0.25) +

  scale_fill_manual(

    values = group_colours,

    breaks = group_levels,

    drop = FALSE,

    name = "Population",

    guide = guide_legend(override.aes = list(shape = 21))

  ) +

  scale_shape_manual(values = host_shapes, name = "Host", drop = FALSE) +

  labs(title = "D. Kenya vs Ethiopia", x = ke_et$xlab, y = ke_et$ylab) +

  coord_cartesian(

    xlim = range(ke_et$pca$PC1) + c(-0.01, 0.01),

    ylim = range(ke_et$pca$PC2) + c(-0.01, 0.01)

  ) +

  theme_pca

# ----------------------------

# Plot E

# ----------------------------

p_ke_et_pc34 <- ggplot(ke_et$pca, aes(PC3, PC4, fill = country, shape = host)) +

  geom_point(size = 2.5, alpha = 0.6, colour = "grey10", stroke = 0.25) +

  scale_fill_manual(

    values = group_colours,

    breaks = group_levels,

    drop = FALSE,

    name = "Population",

    guide = guide_legend(override.aes = list(shape = 21))

  ) +

  scale_shape_manual(values = host_shapes, name = "Host", drop = FALSE) +

  labs(

    title = "E. Kenya vs Ethiopia (PC3 vs PC4)",

    x = paste0("PC3 (", round(ke_et$var_exp[3], 1), "%)"),

    y = paste0("PC4 (", round(ke_et$var_exp[4], 1), "%)")

  ) +

  coord_cartesian(

    xlim = range(ke_et$pca$PC3) + c(-0.005, 0.005),

    ylim = range(ke_et$pca$PC4) + c(-0.005, 0.005)

  ) +

  theme_pca

# ----------------------------

# Combine panels 

# ----------------------------

# Combine subset PCA panels into a single figure and save the output.

pca_panel_5 <- (

  p_af_as | p_af_eu | p_asia_eu

) / (

  p_ke_et | p_ke_et_pc34 | plot_spacer()

) +

  plot_layout(guides = "collect") \&

  theme(legend.position = "right")

pca_panel_5

ggsave(

  filename = "Ascaris_subset_PCA_5panel.png",

  plot = pca_panel_5,

  width = 16,

  height = 10,

  dpi = 600

)
```

#### IBS distance calculations

```text
III) Pairwise 1-IBS genetic distance matrices

# Purpose: calculate pairwise genetic distances from the LD-pruned PLINK datasets using 1 - IBS.

# These distance matrices are used downstream for PERMANOVA, neighbour-joining trees, and IBS

# similarity networks.

# Note: PLINK --distance square 1-ibs outputs 1 - IBS distances, so lower values indicate greater

# genotype similarity.

#!/bin/bash

#SBATCH --job-name=plink_ibsdist_146_65

#SBATCH --partition=high_mem

#SBATCH --cpus-per-task=8

#SBATCH --mem=60G

#SBATCH --time=06:00:00

#SBATCH --output=plink_ibsdist_146_65.%j.out

#SBATCH --error=plink_ibsdist_146_65.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh

conda activate NJ_env

OUT146=/parallel_scratch/lw01165/05_analysis/05d_plink_146

OUT65=/parallel_scratch/lw01165/05_analysis/05e_plink_65

date

echo "Calculating pairwise 1-IBS distance matrices for PERMANOVA, NJ trees and IBS networks"

echo "Running 146 dataset..."

# Produces a square matrix of pairwise 1-IBS distances. Lower distance means higher genotype similarity.

plink \
  --bfile ${OUT146}/Ascaris.146.pruned \
  --distance square 1-ibs \
  --allow-extra-chr \
  --threads 8 \
  --out ${OUT146}/Ascaris.146.pruned

echo "Running 65 dataset..."

# Same IBS-distance calculation for the balanced dataset.

plink \
  --bfile ${OUT65}/Ascaris.65.pruned \
  --distance square 1-ibs \
  --allow-extra-chr \
  --threads 8 \
  --out ${OUT65}/Ascaris.65.pruned

echo "Done."

date
```

#### PERMANOVA and dispersion

```text
IV) PERMANOVA using genome-wide 1-IBS genetic distances

# ============================================================

# PERMANOVA using genome-wide 1-IBS genetic distances

# ============================================================

# This analysis uses the LD-pruned PLINK 1-IBS distance matrices

# generated in the previous section.

#

# Unlike the earlier PCA-based PERMANOVA, this analysis is based

# on pairwise genetic distances across the LD-pruned SNP dataset,

# rather than Euclidean distances calculated from PC1 and PC2 only.

#

# Region and host are tested as explanatory variables. Because

# host and region are partially confounded in the sampling design,

# both single-factor models and combined marginal models are run.

#

# Homogeneity of multivariate dispersion is tested using betadisper.

# This is important because significant PERMANOVA results can reflect

# differences in group centroids, differences in within-group spread,

# or both.

 

# ============================================================

# PERMANOVA on LD-pruned genome-wide 1-IBS genetic distances

# Full n = 146 and balanced n = 65 datasets

# ============================================================

library(dplyr)

library(stringr)

library(readr)

library(vegan)

library(tibble)

# Working directory containing PLINK .mdist and .mdist.id files

setwd("S:/Research/Lauren_PopGen/NJ_tree/ibs")

# Input metadata

meta <- read.csv(

  "S:/Research/Lauren_PopGen/sample_QC_list.csv",

  stringsAsFactors = FALSE

)

names(meta) <- str_trim(names(meta))

meta <- meta %>%

  rename(IID = sample_id) %>%

  mutate(

    IID     = str_trim(IID),

    region  = str_to_lower(str_trim(region)),

    country = str_to_title(str_trim(country)),

    host    = str_to_lower(str_trim(host))

  ) %>%

  mutate(

    region  = na_if(region, ""),

    country = na_if(country, ""),

    host    = na_if(host, "")

  ) %>%

  mutate(

    region = recode(

      region,

      "africa" = "Africa",

      "asia"   = "Asia",

      "europe" = "Europe",

      .default = str_to_title(region)

    ),

    host = recode(

      host,

      "human" = "Human",

      "pig"   = "Pig",

      .default = str_to_title(host)

    )

  ) %>%

  filter(

    !is.na(IID),

    !is.na(region),

    !is.na(host),

    region != "Unknown",

    host != "Unknown"

  )

# Function to read PLINK square 1-IBS distance matrix

read_plink_mdist_matrix <- function(dist_file, id_file) {

  

  ids <- read.table(id_file, header = FALSE, stringsAsFactors = FALSE)

  colnames(ids) <- c("FID", "IID")

  

  dist_mat <- as.matrix(read.table(dist_file, header = FALSE))

  

  rownames(dist_mat) <- ids$IID

  colnames(dist_mat) <- ids$IID

  

  return(dist_mat)

}

# Function to run PERMANOVA + dispersion tests

run_ibs_permanova <- function(dist_file, id_file, meta, dataset_label) {

  

  # Read PLINK 1-IBS distance matrix

  dist_mat <- read_plink_mdist_matrix(dist_file, id_file)

  

  # Keep only samples present in both metadata and distance matrix

  common_ids <- intersect(rownames(dist_mat), meta$IID)

  

  dist_mat <- dist_mat[common_ids, common_ids]

  

  meta_sub <- meta %>%

    filter(IID %in% common_ids) %>%

    arrange(match(IID, common_ids))

  

  # Check order is identical

  stopifnot(all(meta_sub$IID == rownames(dist_mat)))

  

  # Convert square matrix to R distance object

  genetic_dist <- as.dist(dist_mat)

  

  # Sample counts

  region_counts <- table(meta_sub$region)

  host_counts <- table(meta_sub$host)

  cross_counts <- table(meta_sub$region, meta_sub$host)

  

  print(paste("Dataset:", dataset_label))

  print(region_counts)

  print(host_counts)

  print(cross_counts)

  

  # ----------------------------------------------------------

  # PERMANOVA models

  # ----------------------------------------------------------

  

  # Single-factor models

  permanova_region <- adonis2(

    genetic_dist \~ region,

    data = meta_sub,

    permutations = 999

  )

  

  permanova_host <- adonis2(

    genetic_dist \~ host,

    data = meta_sub,

    permutations = 999

  )

  

  # Combined model: sequential terms

  # This tests region first, then host after region.

  permanova_region_host_terms <- adonis2(

    genetic_dist \~ region + host,

    data = meta_sub,

    permutations = 999,

    by = "terms"

  )

  

  # Combined model with reversed order

  # This is useful because region and host are confounded in your dataset.

  permanova_host_region_terms <- adonis2(

    genetic_dist \~ host + region,

    data = meta_sub,

    permutations = 999,

    by = "terms"

  )

  

  # Marginal model

  # This tests each term after accounting for the other term.

  permanova_margin <- adonis2(

    genetic_dist \~ region + host,

    data = meta_sub,

    permutations = 999,

    by = "margin"

  )

  

  # ----------------------------------------------------------

  # Dispersion tests

  # These test whether groups differ in spread around centroids.

  # Significant dispersion means PERMANOVA may partly reflect

  # unequal within-group spread, not only centroid separation.

  # ----------------------------------------------------------

  

  disp_region <- betadisper(genetic_dist, meta_sub$region)

  disp_region_perm <- permutest(disp_region, permutations = 999)

  

  disp_host <- betadisper(genetic_dist, meta_sub$host)

  disp_host_perm <- permutest(disp_host, permutations = 999)

  

  # Helper functions to tidy outputs

  tidy_adonis <- function(x, dataset, model_label) {

    out <- as.data.frame(x)

    out$term <- rownames(out)

    rownames(out) <- NULL

    

    out %>%

      rename(

        df = Df,

        sumsq = SumOfSqs,

        r2 = R2,

        f = F,

        p = `Pr(>F)`

      ) %>%

      mutate(

        dataset = dataset,

        model = model_label

      ) %>%

      select(dataset, model, term, df, sumsq, r2, f, p)

  }

  

  tidy_disp <- function(x, dataset, factor_label) {

    tab <- as.data.frame(x$tab)

    

    tibble(

      dataset = dataset,

      analysis = "Dispersion",

      factor = factor_label,

      df = tab$Df[1],

      f = tab$F[1],

      p = tab$`Pr(>F)`[1]

    )

  }

  

  permanova_table <- bind_rows(

    tidy_adonis(permanova_region, dataset_label, "single_factor_region"),

    tidy_adonis(permanova_host, dataset_label, "single_factor_host"),

    tidy_adonis(permanova_region_host_terms, dataset_label, "combined_region_plus_host_terms"),

    tidy_adonis(permanova_host_region_terms, dataset_label, "combined_host_plus_region_terms"),

    tidy_adonis(permanova_margin, dataset_label, "combined_region_plus_host_margin")

  )

  

  dispersion_table <- bind_rows(

    tidy_disp(disp_region_perm, dataset_label, "region"),

    tidy_disp(disp_host_perm, dataset_label, "host")

  )

  

  return(list(

    meta = meta_sub,

    dist_matrix = dist_mat,

    permanova_region = permanova_region,

    permanova_host = permanova_host,

    permanova_region_host_terms = permanova_region_host_terms,

    permanova_host_region_terms = permanova_host_region_terms,

    permanova_margin = permanova_margin,

    disp_region = disp_region_perm,

    disp_host = disp_host_perm,

    permanova_table = permanova_table,

    dispersion_table = dispersion_table,

    region_counts = region_counts,

    host_counts = host_counts,

    cross_counts = cross_counts

  ))

}

# Run analyses

res_146 <- run_ibs_permanova(

  dist_file = "Ascaris.146.pruned.mdist",

  id_file   = "Ascaris.146.pruned.mdist.id",

  meta = meta,

  dataset_label = "Full_146"

)

res_65 <- run_ibs_permanova(

  dist_file = "Ascaris.65.pruned.mdist",

  id_file   = "Ascaris.65.pruned.mdist.id",

  meta = meta,

  dataset_label = "Balanced_65"

)

# Print raw outputs

print(res_146$permanova_region)

print(res_146$permanova_host)

print(res_146$permanova_region_host_terms)

print(res_146$permanova_host_region_terms)

print(res_146$permanova_margin)

print(res_146$disp_region)

print(res_146$disp_host)

print(res_65$permanova_region)

print(res_65$permanova_host)

print(res_65$permanova_region_host_terms)

print(res_65$permanova_host_region_terms)

print(res_65$permanova_margin)

print(res_65$disp_region)

print(res_65$disp_host)

# ============================================================

# Combine and clean tables

# ============================================================

permanova_results <- bind_rows(

  res_146$permanova_table,

  res_65$permanova_table

)

dispersion_results <- bind_rows(

  res_146$dispersion_table,

  res_65$dispersion_table

)

permanova_results_clean <- permanova_results %>%

  filter(term %in% c("region", "host", "Model")) %>%

  mutate(

    r2 = round(r2, 4),

    f = round(f, 3),

    p = ifelse(

      is.na(p),

      NA,

      ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

    )

  )

dispersion_results_clean <- dispersion_results %>%

  mutate(

    f = round(f, 3),

    p = ifelse(

      is.na(p),

      NA,

      ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

    )

  )

summary_table <- bind_rows(

  permanova_results_clean %>%

    transmute(

      dataset,

      analysis = "PERMANOVA_1IBS",

      model,

      factor = term,

      df,

      r2,

      f,

      p

    ),

  dispersion_results_clean %>%

    transmute(

      dataset,

      analysis = "Dispersion",

      model = NA_character_,

      factor,

      df,

      r2 = NA_real_,

      f,

      p

    )

)

print(summary_table)

# ============================================================

# Save outputs

# ============================================================

write.csv(

  permanova_results_clean,

  "PERMANOVA_1IBS_results_clean.csv",

  row.names = FALSE

)

write.csv(

  dispersion_results_clean,

  "PERMANOVA_1IBS_dispersion_results_clean.csv",

  row.names = FALSE

)

write.csv(

  summary_table,

  "PERMANOVA_1IBS_summary_table.csv",

  row.names = FALSE

)
```

#### NJ trees and IBS networks

```text
V) Neighbour-joining trees and IBS network analyses

# Purpose: use the pairwise 1-IBS distance matrices generated above for distance-based tree building

# and exploratory IBS similarity network visualisation.

# A) Neighbour-joining tree analysis

# Purpose: convert PLINK distance matrices into NJ trees, midpoint-root them, and plot sample metadata

# on the tree tips.

library(ape)

library(ggtree)

library(ggplot2)

library(dplyr)

library(stringr)

library(readr)

library(phytools)   # midpoint.root()

setwd("S:/Research/Lauren_PopGen/NJ_tree/ibs")

# ----------------------------

# Metadata

# ----------------------------

meta <- read.csv(

  "S:/Research/Lauren_PopGen/sample_QC_list.csv",

  stringsAsFactors = FALSE

)

names(meta) <- str_trim(names(meta))

meta <- meta %>%

  rename(label = sample_id)

meta <- meta %>%

  mutate(

    label   = str_trim(label),

    region  = str_to_lower(str_trim(region)),

    country = str_trim(country),

    host    = str_to_lower(str_trim(host))

  ) %>%

  mutate(

    region  = na_if(region, ""),

    country = na_if(country, ""),

    host    = na_if(host, "")

  ) %>%

  mutate(

    region = recode(

      region,

      "africa" = "Africa",

      "asia"   = "Asia",

      "europe" = "Europe",

      .default = str_to_title(region)

    ),

    region = if_else(is.na(region), "Unknown", region),

    host   = recode(

      host,

      "human" = "Human",

      "pig"   = "Pig",

      .default = "Unknown"

    )

  ) %>%

  mutate(

    region_country = if_else(region == "Africa", country, region)

  )

af_countries <- sort(unique(meta$country[meta$region == "Africa" \& !is.na(meta$country)]))

meta$region_country <- factor(

  meta$region_country,

  levels = c("Asia", "Europe", "Unknown", af_countries)

)

meta$host <- factor(

  meta$host,

  levels = c("Human", "Pig", "Unknown")

)

# ----------------------------

# Read PLINK mdist

# ----------------------------

read_plink_mdist <- function(dist_file, id_file) {

  ids <- read.table(id_file, header = FALSE, stringsAsFactors = FALSE)

  colnames(ids) <- c("FID", "IID")

  

  dist_mat <- as.matrix(read.table(dist_file, header = FALSE))

  rownames(dist_mat) <- ids$IID

  colnames(dist_mat) <- ids$IID

  

  as.dist(dist_mat)

}

# ----------------------------

# Input files

# ----------------------------

dist_file_146 <- "Ascaris.146.pruned.mdist"

id_file_146   <- "Ascaris.146.pruned.mdist.id"

# ----------------------------

# Build tree

# ----------------------------

dist_146 <- read_plink_mdist(dist_file_146, id_file_146)

# Build a neighbour-joining tree from pairwise IBS distances. This is a distance-based clustering tree,

# not a model-based phylogeny.

tree_146 <- nj(dist_146)

# midpoint root + ladderize

tree_146 <- midpoint.root(tree_146)

tree_146 <- ladderize(tree_146, right = FALSE)

# ----------------------------

# Keep matching metadata only

# ----------------------------

meta_146 <- meta %>%

  filter(label %in% tree_146$tip.label)

p_base <- ggtree(tree_146, size = 0.25) %<+% meta_146 +

  geom_tiplab(

    align = TRUE,

    linetype = "dotted",

    linesize = 0.2,

    size = 1.2,

    offset = 0.0015

  )

# extract plotted tree data

tree_dat <- p_base$data

tip_dat <- tree_dat %>%

  filter(isTip) %>%

  mutate(x_meta = max(x) + 0.0035)

# ----------------------------

# Plot

# ----------------------------

p146 <- ggtree(tree_146, size = 0.25) %<+% meta_146 +

  

  # aligned tip labels

  geom_tiplab(

    align = TRUE,

    linetype = "dotted",

    linesize = 0.15,

    size = 2.0,

    offset = 0.0015

  ) +

  

  # metadata points moved to the right of the aligned labels

  geom_tippoint(

    aes(color = region_country, shape = host),

    size = 1.8,

    na.rm = TRUE,

    offset = 0.002

  ) +

  

  scale_colour_manual(

    values = c(

      "Asia"     = "#F8766D",

      "Europe"   = "#B79F00",

      "Unknown"  = "#00BA38",

      "Kenya"    = "#0047AB",

      "Ethiopia" = "#00BFC4",

      "Uganda"   = "#F564E3"

    )) +

  

  scale_shape_manual(values = c(

    "Human"   = 16,

    "Pig"     = 17,

    "Unknown" = 15

  )) +

  

  theme_tree2() +

  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +

  coord_cartesian(xlim = c(0, 0.09), clip = "off") +

  

  labs(

    title = "Neighbour-joining tree (IBS distance, n = 146)",

    color = "Region / country",

    shape = "Host"

  ) +

  

  theme(

    plot.title = element_text(hjust = 0.5, face = "bold"),

    plot.subtitle = element_text(hjust = 0.5),

    

    # move legend inside, closer to tree

    legend.position = c(0.80, 0.52),

    legend.justification = c(0, 0.5),

    

    # tidy legend spacing

    legend.background = element_blank(),

    legend.key = element_blank(),

    legend.box.background = element_blank(),

    

    # less wasted space on the right

    plot.margin = margin(10, 20, 10, 10)

  )

p146

ggsave(

  filename = "NJ_tree_146_aligned_labels.png",

  plot = p146,

  width = 12,

  height = 14,

  units = "in",

  dpi = 600,

  bg = "white"

)

# ----------------------------

# Input files: balanced dataset (65)

# ----------------------------

dist_file_65 <- "Ascaris.65.pruned.mdist"

id_file_65   <- "Ascaris.65.pruned.mdist.id"

# ----------------------------

# Build tree

# ----------------------------

dist_65 <- read_plink_mdist(dist_file_65, id_file_65)

# Build the same distance-based NJ tree for the balanced dataset.

tree_65 <- nj(dist_65)

# midpoint root + ladderize

tree_65 <- midpoint.root(tree_65)

# Midpoint rooting is used because no external outgroup is included. Root position should therefore be

# interpreted as a visualisation choice, not direct evidence of ancestry direction.

tree_65 <- ladderize(tree_65, right = FALSE)

# ----------------------------

# Keep matching metadata only

# ----------------------------

meta_65 <- meta %>%

  filter(label %in% tree_65$tip.label)

# ----------------------------

# Plot

# ----------------------------

p65 <- ggtree(tree_65, size = 0.3) %<+% meta_65 +

  

  geom_tiplab(

    align = TRUE,

    linetype = "dotted",

    linesize = 0.2,

    size = 2.2,

    offset = 0.0015

  ) +

  

  geom_tippoint(

    aes(color = region_country, shape = host),

    size = 2.2,

    na.rm = TRUE,

    offset = 0.002

  ) +

  

  scale_colour_manual(

    values = c(

      "Asia"     = "#F8766D",

      "Europe"   = "#B79F00",

      "Unknown"  = "#00BA38",

      "Ethiopia" = "#00BFC4",

      "Kenya"    = "#0047AB",

      "Uganda"   = "#F564E3"

    )

  ) +

  

  scale_shape_manual(values = c(

    "Human"   = 16,

    "Pig"     = 17,

    "Unknown" = 15

  )) +

  

  theme_tree2() +

  coord_cartesian(xlim = c(0, 0.095), clip = "off") +

  

  labs(

    title = "Neighbour-joining tree (IBS distance, n = 65)",

    color = "Region / country",

    shape = "Host"

  ) +

  

  theme(

    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),

    

    axis.text.y = element_blank(),

    axis.ticks.y = element_blank(),

    

    legend.position = "right",

    legend.box = "vertical",

    legend.title = element_text(size = 12),

    legend.text = element_text(size = 10),

    legend.key = element_blank(),

    legend.background = element_blank(),

    

    plot.margin = margin(15, 140, 15, 15)

  ) +

  

  guides(

    shape = guide_legend(order = 1),

    colour = guide_legend(order = 2)

  )

p65

ggsave(

  filename = "NJ_tree_65_aligned_labels.png",

  plot = p65,

  width = 12,

  height = 16,

  units = "in",

  dpi = 600,

  bg = "white"

)

library(patchwork)

combined_plot <- p146 + p65 +

  plot_layout(ncol = 2, guides = "collect") \&

  theme(

    legend.position = "right"

  )

combined_plot

ggsave(

  filename = "NJ_tree_combined_146_65.png",

  plot = combined_plot,

  width = 18,

  height = 14,

  units = "in",

  dpi = 600,

  bg = "white"

)

#B) IBS network analysis 

# Purpose: visualise the strongest pairwise genotype similarities as networks. Nodes are isolates;

# edges connect pairs above IBS percentile thresholds.

library(dplyr)

library(readr)

library(stringr)

library(tidyr)

library(igraph)

library(ggraph)

library(ggplot2)

library(patchwork)

setwd("S:/Research/Lauren_PopGen/NJ_tree/ibs")

# -----------------------------

# Metadata

# -----------------------------

meta <- read_csv("S:/Research/Lauren_PopGen/sample_QC_list.csv", show_col_types = FALSE)

names(meta) <- str_trim(names(meta))

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

  ) %>%

  mutate(

    region = recode(

      region,

      "africa" = "Africa",

      "asia"   = "Asia",

      "europe" = "Europe",

      .default = str_to_title(region)

    ),

    region = if_else(is.na(region), "Unknown", region),

    host = recode(

      host,

      "human" = "Human",

      "pig"   = "Pig",

      .default = "Unknown"

    ),

    plot_group = if_else(region == "Africa" \& !is.na(country), country, region)

  )

# -----------------------------

# Read PLINK mdist

# -----------------------------

read_plink_mdist_matrix <- function(dist_file, id_file) {

  ids <- read.table(id_file, header = FALSE, stringsAsFactors = FALSE)

  colnames(ids) <- c("FID", "IID")

  

  dist_mat <- as.matrix(read.table(dist_file, header = FALSE))

  rownames(dist_mat) <- ids$IID

  colnames(dist_mat) <- ids$IID

  

  dist_mat

}

dist_mat <- read_plink_mdist_matrix(

  "Ascaris.146.pruned.mdist",

  "Ascaris.146.pruned.mdist.id"

)

# -----------------------------

# Convert to IBS and edge list

# -----------------------------

# Convert 1-IBS distance back to IBS similarity, so larger values indicate more similar genotypes.

ibs_mat <- 1 - dist_mat

diag(ibs_mat) <- NA

edges <- as.data.frame(as.table(ibs_mat), stringsAsFactors = FALSE)

colnames(edges) <- c("from", "to", "ibs")

edges <- edges %>%

  filter(!is.na(ibs)) %>%

  rowwise() %>%

  mutate(pair = paste(sort(c(from, to)), collapse = "__")) %>%

  ungroup() %>%

  distinct(pair, .keep_all = TRUE) %>%

  select(-pair)

# -----------------------------

# Thresholds

# -----------------------------

# Percentile thresholds retain only the strongest pairwise IBS similarities. These are exploratory,

# data-relative cut-offs rather than biological relatedness thresholds.

thr80 <- quantile(edges$ibs, 0.80, na.rm = TRUE)  # top 20%

thr90 <- quantile(edges$ibs, 0.90, na.rm = TRUE)  # top 10%

thr98 <- quantile(edges$ibs, 0.98, na.rm = TRUE)  # top 2%

thr99 <- quantile(edges$ibs, 0.99, na.rm = TRUE)  # top 1%

cat("80th percentile threshold (top 20%):", round(thr80, 4), "\\n")

cat("90th percentile threshold (top 10%):", round(thr90, 4), "\\n")

cat("98th percentile threshold (top 2%):", round(thr98, 4), "\\n")

cat("99th percentile threshold (top 1%):", round(thr99, 4), "\\n")

# -----------------------------

# Edge sets

# -----------------------------

edges_20 <- edges %>% filter(ibs >= thr80)

edges_10 <- edges %>% filter(ibs >= thr90)

edges_2  <- edges %>% filter(ibs >= thr98)

edges_1  <- edges %>% filter(ibs >= thr99)

cat("Top 20% edges:", nrow(edges_20), "\\n")

cat("Top 10% edges:", nrow(edges_10), "\\n")

cat("Top 2% edges:",  nrow(edges_2),  "\\n")

cat("Top 1% edges:",  nrow(edges_1),  "\\n")

# -----------------------------

# Graph builder

# -----------------------------

make_graph <- function(edge_df, meta_df) {

  nodes <- meta_df %>%

    distinct(IID, plot_group, host) %>%

    rename(name = IID)

  

  graph_from_data_frame(

    d = edge_df,

    vertices = nodes,

    directed = FALSE

  )

}

g_20 <- make_graph(edges_20, meta)

g_10 <- make_graph(edges_10, meta)

g_2  <- make_graph(edges_2,  meta)

g_1  <- make_graph(edges_1,  meta)

# remove isolates

g_20 <- delete_vertices(g_20, degree(g_20) == 0)

g_10 <- delete_vertices(g_10, degree(g_10) == 0)

g_2  <- delete_vertices(g_2,  degree(g_2)  == 0)

g_1  <- delete_vertices(g_1,  degree(g_1)  == 0)

# -----------------------------

# Plot function

# Each panel gets its own layout

# -----------------------------

plot_ibs_network <- function(g, title, subtitle, show_legend = FALSE) {

  set.seed(123)

  

  ggraph(g, layout = "fr") +

    geom_edge_link(

      alpha = 0.16,

      colour = "grey20",

      linewidth = 0.45

    ) +

    geom_node_point(

      aes(colour = plot_group, shape = host),

      size = 2.8

    ) +

    scale_colour_manual(

      values = c(

        "Asia"     = "#F8766D",

        "Europe"   = "gold",

        "Unknown"  = "#00BA38",

        "Kenya"    = "#0047AB",

        "Ethiopia" = "#00BFC4",

        "Uganda"   = "#F564E3"

      ),

      drop = FALSE

    ) +

    scale_shape_manual(

      values = c(

        "Human"   = 16,

        "Pig"     = 17,

        "Unknown" = 15

      ),

      drop = FALSE

    ) +

    theme_void() +

    labs(

      title = title,

      subtitle = subtitle,

      colour = "Region / country",

      shape  = "Host"

    ) +

    theme(

      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),

      plot.subtitle = element_text(hjust = 0.5, size = 9),

      legend.position = if (show_legend) "right" else "none"

    )

}

# -----------------------------

# Individual panels

# -----------------------------

p_20 <- plot_ibs_network(

  g_20,

  title = "A. Top 20% of IBS values",

  subtitle = paste0("Threshold = ", round(thr80, 3), "; edges = ", nrow(edges_20)),

  show_legend = TRUE

)

p_10 <- plot_ibs_network(

  g_10,

  title = "B. Top 10% of IBS values",

  subtitle = paste0("Threshold = ", round(thr90, 3), "; edges = ", nrow(edges_10)),

  show_legend = FALSE

)

p_2 <- plot_ibs_network(

  g_2,

  title = "C. Top 2% of IBS values",

  subtitle = paste0("Threshold = ", round(thr98, 3), "; edges = ", nrow(edges_2)),

  show_legend = FALSE

)

p_1 <- plot_ibs_network(

  g_1,

  title = "D. Top 1% of IBS values",

  subtitle = paste0("Threshold = ", round(thr99, 3), "; edges = ", nrow(edges_1)),

  show_legend = FALSE

)

# -----------------------------

# Combine panel

# -----------------------------

panel_plot <- ((p_20 | p_10) / plot_spacer() / (p_2 | p_1)) +

  plot_layout(

    guides = "collect",

    heights = c(1, 0.25, 1)

  ) \&

  theme(

    legend.position = "right",

    legend.box = "vertical",

    legend.spacing.y = unit(0.35, "cm")

  )

panel_plot <- panel_plot +

  plot_annotation(

    title = "IBS similarity networks across thresholds",

    subtitle = "Nodes represent Ascaris isolates; edges represent pairwise IBS similarity",

    theme = theme(

      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),

      plot.subtitle = element_text(hjust = 0.5, size = 11)

    )

  )

panel_plot

# -----------------------------

# Save

# -----------------------------

ggsave(

  filename = "IBS_network_panel_20_10_2_1_final.png",

  plot = panel_plot,

  width = 15,

  height = 12,

  units = "in",

  dpi = 600,

  bg = "white"

)
```

#### Maximum-likelihood phylogeny

```text
VI) Maximum Likelihood Tree 

# Purpose: build a SNP-based maximum-likelihood tree for the balanced dataset.

# Workflow: clean to biallelic SNPs -> randomly sample 100,000 SNPs -> convert genotypes to a FASTA

# alignment -> convert FASTA to PHYLIP -> run IQ-TREE -> plot with metadata in R.

#a) ensure data are clean

#!/bin/bash

#SBATCH --job-name=ascaris_clean_vcf

#SBATCH --partition=high_mem

#SBATCH --cpus-per-task=6

#SBATCH --mem=24G

#SBATCH --time=12:00:00

#SBATCH --output=ascaris_clean_vcf.%j.out

#SBATCH --error=ascaris_clean_vcf.%j.err

set -eo pipefail

echo "===== JOB START: CLEAN VCF ====="

date

echo "Host: $(hostname)"

source /users/lw01165/miniconda3/etc/profile.d/conda.sh

conda activate tools_env

WORKDIR=/parallel_scratch/lw01165/05_analysis/05M_ML_tree

mkdir -p "${WORKDIR}"

cd "${WORKDIR}"

# Input VCF: balanced 65-sample dataset

VCF_IN=/parallel_scratch/lw01165/05_analysis/05a_Balanced_dataset/Ascaris.autosomes.SNPs.146.maxmiss0.90.balanced.n65.vcf.gz

# Output cleaned VCF

VCF_CLEAN=Ascaris.65.biallelic.snps.vcf.gz

echo "===== CHECK INPUTS ====="

[[ -f "${VCF_IN}" ]] || { echo "Missing input VCF: ${VCF_IN}"; exit 1; }

echo "Input VCF: ${VCF_IN}"

echo "Output cleaned VCF: ${VCF_CLEAN}"

echo "===== CHECK SOFTWARE ====="

which bcftools || { echo "bcftools not found"; exit 1; }

which tabix || { echo "tabix not found"; exit 1; }

bcftools --version | head -n 1

echo "===== STEP 1: FILTER TO BIALLELIC SNPs ONLY ====="

date

# Keep only biallelic SNPs. This removes indels and multiallelic variants before constructing the

# SNP alignment.

bcftools view \
  -m2 -M2 -v snps \
  "${VCF_IN}" \
  -Oz -o "${VCF_CLEAN}"

echo "===== STEP 2: INDEX CLEANED VCF ====="

date

tabix -f -p vcf "${VCF_CLEAN}"

echo "===== CHECK OUTPUTS ====="

date

ls -lh "${VCF_CLEAN}" "${VCF_CLEAN}.tbi"

echo "===== JOB COMPLETE: CLEAN VCF ====="

date

#B) Thin to 100k SNPS

#!/bin/bash

#SBATCH --job-name=ascaris_thin_exact100k

#SBATCH --partition=high_mem

#SBATCH --cpus-per-task=16

#SBATCH --mem=80G

#SBATCH --time=04:00:00

#SBATCH --output=ascaris_thin_exact100k.%j.out

#SBATCH --error=ascaris_thin_exact100k.%j.err

set -eo pipefail

echo "===== JOB START: THIN TO EXACT 100K SNPs ====="

date

echo "Host: $(hostname)"

source /users/lw01165/miniconda3/etc/profile.d/conda.sh

conda activate tools_env

WORKDIR=/parallel_scratch/lw01165/05_analysis/05M_ML_tree

cd "${WORKDIR}"

VCF_IN="${WORKDIR}/Ascaris.65.biallelic.snps.vcf.gz"

VCF_OUT="${WORKDIR}/Ascaris.65.exact100000.vcf.gz"

TMP_HEADER="${WORKDIR}/tmp.header.${SLURM_JOB_ID}.vcf"

TMP_BODY="${WORKDIR}/tmp.body.${SLURM_JOB_ID}.vcf"

TMP_IDX="${WORKDIR}/tmp.keepidx.${SLURM_JOB_ID}.txt"

TARGET=100000

# Number of SNPs to sample for the ML alignment. The random seed below makes this subset reproducible.

echo "===== CHECK INPUT ====="

[[ -f "${VCF_IN}" ]] || { echo "Missing input: ${VCF_IN}"; exit 1; }

which bcftools >/dev/null 2>\&1 || { echo "bcftools not found"; exit 1; }

which bgzip >/dev/null 2>\&1 || { echo "bgzip not found"; exit 1; }

which tabix >/dev/null 2>\&1 || { echo "tabix not found"; exit 1; }

TOTAL=$(bcftools index -n "${VCF_IN}")

echo "Total SNPs in input: ${TOTAL}"

if [[ "${TOTAL}" -lt "${TARGET}" ]]; then

    echo "Input has fewer than ${TARGET} SNPs; cannot sample exactly ${TARGET}"

    exit 1

fi

rm -f "${VCF_OUT}" "${VCF_OUT}.tbi" "${TMP_HEADER}" "${TMP_BODY}" "${TMP_IDX}"

echo "===== MAKE HEADER ====="

bcftools view -h "${VCF_IN}" > "${TMP_HEADER}"

echo "===== CHOOSE EXACT RECORD INDICES ====="

python3 - <<PY

import random

random.seed(42)

# Fixed seed makes the 100,000-SNP random subset reproducible.

total = ${TOTAL}

target = ${TARGET}

idx = sorted(random.sample(range(1, total + 1), target))

with open("${TMP_IDX}", "w") as f:

    for i in idx:

        f.write(f"{i}\\n")

PY

echo "===== EXTRACT EXACT SNP SET ====="

bcftools view -H "${VCF_IN}" | awk 'NR==FNR{keep[$1]=1; next} keep[FNR]' "${TMP_IDX}" - > "${TMP_BODY}"

KEPT=$(wc -l < "${TMP_BODY}")

echo "Records retained: ${KEPT}"

[[ "${KEPT}" -eq "${TARGET}" ]] || { echo "Did not retain exactly ${TARGET} SNPs"; exit 1; }

echo "===== WRITE OUTPUT ====="

cat "${TMP_HEADER}" "${TMP_BODY}" | bgzip -c > "${VCF_OUT}"

tabix -f -p vcf "${VCF_OUT}"

echo "===== FINAL CHECK ====="

bcftools index -n "${VCF_OUT}"

ls -lh "${VCF_OUT}" "${VCF_OUT}.tbi"

rm -f "${TMP_HEADER}" "${TMP_BODY}" "${TMP_IDX}"

echo "===== JOB COMPLETE ====="

date

# C) turn VCF to .fasta format 

#!/bin/bash

#SBATCH --job-name=ascaris_vcf_to_fasta_clean

#SBATCH --partition=high_mem

#SBATCH --cpus-per-task=8

#SBATCH --mem=50G

#SBATCH --time=04:00:00

#SBATCH --output=ascaris_vcf_to_fasta_clean.%j.out

#SBATCH --error=ascaris_vcf_to_fasta_clean.%j.err

set -eo pipefail

echo "===== JOB START: VCF TO CLEAN FASTA ====="

date

echo "Host: $(hostname)"

source /users/lw01165/miniconda3/etc/profile.d/conda.sh

conda activate tools_env

WORKDIR=/parallel_scratch/lw01165/05_analysis/05M_ML_tree

cd "${WORKDIR}"

VCF_IN="${WORKDIR}/Ascaris.65.exact100000.vcf.gz"

FASTA_OUT="${WORKDIR}/Ascaris.65.exact100000.clean.fasta"

echo "===== CHECK INPUT ====="

[[ -f "${VCF_IN}" ]] || { echo "Missing input: ${VCF_IN}"; exit 1; }

which bcftools >/dev/null 2>\&1 || { echo "bcftools not found"; exit 1; }

python3 - <<'PY'

import gzip

import re

vcf = "Ascaris.65.exact100000.vcf.gz"

fasta = "Ascaris.65.exact100000.clean.fasta"

# Convert diploid genotypes to a single consensus base for tree building. Homozygous reference and

# homozygous alternate calls are encoded as their nucleotide; heterozygous and missing calls become N.

def gt_to_base(gt, ref, alt):

    gt = gt.replace("|", "/")

    if gt in ("0/0",):

        return ref

    elif gt in ("1/1",):

        return alt

    elif gt in ("0/1", "1/0"):

        return "N"

    elif gt in ("./.", "."):

        return "N"

    else:

        return "N"

samples = []

seqs = {}

with gzip.open(vcf, "rt") as f:

    for line in f:

        if line.startswith("##"):

            continue

        if line.startswith("#CHROM"):

            header = line.rstrip("\\n").split("\\t")

            samples = header[9:]

            seqs = {s: [] for s in samples}

            continue

        fields = line.rstrip("\\n").split("\\t")

        ref = fields[3]

        alt = fields[4]

        # enforce simple biallelic SNPs only

        if len(ref) != 1 or len(alt) != 1 or "," in alt:

            continue

        genos = fields[9:]

        for sample, entry in zip(samples, genos):

            gt = entry.split(":")[0]

            seqs[sample].append(gt_to_base(gt, ref, alt))

with open(fasta, "w") as out:

    for s in samples:

        out.write(f">{s}\\n")

        seq = "".join(seqs[s])

        for i in range(0, len(seq), 80):

            out.write(seq[i:i+80] + "\\n")

print(f"Wrote FASTA for {len(samples)} samples")

print(f"Alignment length: {len(next(iter(seqs.values())))}")

PY

echo "===== CHECK OUTPUT ====="

ls -lh "${FASTA_OUT}"

head -4 "${FASTA_OUT}"

echo "===== JOB COMPLETE ====="

date

#D) .fasta to Phylip format

#!/bin/bash

#SBATCH --job-name=ascaris_fasta_to_phylip_clean

#SBATCH --partition=high_mem

#SBATCH --cpus-per-task=4

#SBATCH --mem=30G

#SBATCH --time=01:00:00

#SBATCH --output=ascaris_fasta_to_phylip_clean.%j.out

#SBATCH --error=ascaris_fasta_to_phylip_clean.%j.err

set -eo pipefail

echo "===== JOB START: FASTA TO PHYLIP ====="

date

echo "Host: $(hostname)"

source /users/lw01165/miniconda3/etc/profile.d/conda.sh

conda activate tools_env

WORKDIR=/parallel_scratch/lw01165/05_analysis/05M_ML_tree

cd "${WORKDIR}"

FASTA_IN="${WORKDIR}/Ascaris.65.exact100000.clean.fasta"

PHY_OUT="${WORKDIR}/Ascaris.65.exact100000.clean.phy"

python3 - <<'PY'

fasta = "Ascaris.65.exact100000.clean.fasta"

phy = "Ascaris.65.exact100000.clean.phy"

seqs = {}

name = None

with open(fasta) as f:

    for line in f:

        line = line.strip()

        if not line:

            continue

        if line.startswith(">"):

            name = line[1:]

            seqs[name] = []

        else:

            seqs[name].append(line)

seqs = {k: "".join(v) for k, v in seqs.items()}

names = list(seqs.keys())

nseq = len(names)

lengths = {len(v) for v in seqs.values()}

if len(lengths) != 1:

    raise ValueError("Sequences are not all the same length")

nsite = lengths.pop()

with open(phy, "w") as out:

    out.write(f"{nseq} {nsite}\\n")

    for name in names:

        out.write(f"{name} {seqs[name]}\\n")

print(f"Wrote PHYLIP: {nseq} taxa, {nsite} sites")

PY

echo "===== CHECK OUTPUT ====="

ls -lh "${PHY_OUT}"

head -3 "${PHY_OUT}"

echo "===== JOB COMPLETE ====="

date

#E) Run IQ Tree

#!/bin/bash

#SBATCH --job-name=ascaris_iqtree_ASC

#SBATCH --partition=high_mem

#SBATCH --cpus-per-task=20

#SBATCH --mem=120G

#SBATCH --time=3-00:00:00

#SBATCH --output=ascaris_iqtree_ASC.%j.out

#SBATCH --error=ascaris_iqtree_ASC.%j.err

set -eo pipefail

echo "===== JOB START: IQ-TREE ASC-CORRECTED SNP ALIGNMENT ====="

date

echo "Host: $(hostname)"

source /users/lw01165/miniconda3/etc/profile.d/conda.sh

conda activate iqtree

WORKDIR=/parallel_scratch/lw01165/05_analysis/05M_ML_tree

cd "${WORKDIR}"

PHY_IN="${WORKDIR}/Ascaris.65.exact100000.clean.phy"

PREFIX="Ascaris.65.exact100000.clean.iqtree.ASC"

[[ -f "${PHY_IN}" ]] || { echo "Missing PHYLIP: ${PHY_IN}"; exit 1; }

if command -v iqtree2 >/dev/null 2>\&1; then

    IQT=$(command -v iqtree2)

elif command -v iqtree >/dev/null 2>\&1; then

    IQT=$(command -v iqtree)

else

    echo "No IQ-TREE binary found"

    exit 1

fi

echo "Using IQ-TREE: ${IQT}"

echo "Input PHYLIP: ${PHY_IN}"

echo "Output prefix: ${PREFIX}"

rm -f "${PREFIX}".*

# Run IQ-TREE with ModelFinder, ultrafast bootstrap, and SH-aLRT branch support.

# Because this alignment contains SNPs only and excludes invariant sites, ModelFinder is run

# with +ASC so that ascertainment-bias-corrected models are considered. This avoids treating

# the SNP-only alignment as though it represented a full genome alignment containing constant sites.

"${IQT}" \
  -s "${PHY_IN}" \
  -st DNA \
  -m MFP+ASC \
  -nt ${SLURM_CPUS_PER_TASK} \
  -bb 1000 \
  -alrt 1000 \
  -pre "${PREFIX}"

echo "===== CHECK OUTPUT ====="

ls -lh "${PREFIX}".*

echo "===== JOB COMPLETE ====="

date

#F) Analysis in RStudio. 

# Purpose: midpoint-root the IQ-TREE output, join tip labels to metadata, extract branch support values,

# and save a labelled ML tree figure.

library(ape)

library(ggtree)

library(ggplot2)

library(dplyr)

library(stringr)

library(phytools)

library(tibble)

setwd("S:/Research/Lauren_PopGen/ML_tree")

meta <- read.csv(

  "S:/Research/Lauren_PopGen/sample_QC_list.csv",

  stringsAsFactors = FALSE,

  check.names = FALSE

)

names(meta) <- str_trim(names(meta))

meta <- meta %>%

  mutate(

    label   = str_trim(sample_id),

    region  = str_to_lower(str_trim(region)),

    country = str_trim(country),

    host    = str_to_lower(str_trim(host))

  ) %>%

  mutate(

    region  = na_if(region, ""),

    country = na_if(country, ""),

    host    = na_if(host, "")

  ) %>%

  mutate(

    region = recode(

      region,

      "africa" = "Africa",

      "asia"   = "Asia",

      "europe" = "Europe",

      .default = str_to_title(region)

    ),

    region = if_else(is.na(region), "Unknown", region),

    host   = recode(

      host,

      "human" = "Human",

      "pig"   = "Pig",

      .default = "Unknown"

    ),

    region_country = if_else(region == "Africa", country, region)

  )

af_countries <- sort(unique(meta$country[meta$region == "Africa" \& !is.na(meta$country)]))

meta$region_country <- factor(

  meta$region_country,

  levels = c("Asia", "Europe", "Unknown", af_countries)

)

meta$host <- factor(

  meta$host,

  levels = c("Human", "Pig", "Unknown")

)

tree_65 <- read.tree("Ascaris.65.exact100000.clean.iqtree.ASC.treefile")

tree_65 <- midpoint.root(tree_65)

# Midpoint rooting is used because no external outgroup is included. Root position should therefore be

# interpreted as a visualisation choice, not direct evidence of ancestry direction.

meta_65 <- meta %>%

  filter(label %in% tree_65$tip.label) %>%

  select(label, region_country, host)

# IQ-TREE labels produced by -alrt and -bb are commonly formatted as SH-aLRT/UFBoot; the code below

# extracts the UFBoot value after the slash.

bs_df <- tibble(

  node = seq(length(tree_65$tip.label) + 1,

             length(tree_65$tip.label) + tree_65$Nnode),

  node_label_raw = tree_65$node.label

) %>%

  mutate(

    ufboot = suppressWarnings(as.numeric(str_extract(node_label_raw, "(?<=/)[^/]+$")))

  )

p0 <- ggtree(tree_65, size = 0.35)

tree_dat <- p0$data %>%

  left_join(meta_65, by = "label") %>%

  left_join(bs_df, by = "node")

bs_plot <- tree_dat %>%

  filter(!isTip, !is.na(ufboot), ufboot >= 80, !is.na(parent)) %>%

  select(node, parent, x, y, ufboot) %>%

  left_join(

    tree_dat %>%

      select(node, x_parent = x),

    by = c("parent" = "node")

  ) %>%

  mutate(

    x_mid = x_parent + 0.5 * (x - x_parent),

    y_mid = y - 0.25

  )

pML <- ggtree(tree_65, size = 0.35) %<+% tree_dat +

  geom_tiplab(

    align = TRUE,

    linetype = "dotted",

    linesize = 0.12,

    colour = "black",

    size = 3.4,

    offset = 0.002

  ) +

  geom_tippoint(

    aes(color = region_country, shape = host),

    size = 3.2,

    na.rm = TRUE

  ) +

  geom_text(

    data = bs_plot,

    mapping = aes(x = x_mid, y = y_mid, label = round(ufboot, 0)),

    inherit.aes = FALSE,

    size = 2.5,

    colour = "black"

  ) +

  scale_colour_manual(

    values = c(

      "Asia"     = "#F8766D",

      "Europe"   = "#B79F00",

      "Unknown"  = "#00BA38",

      "Ethiopia" = "#00BFC4",

      "Kenya"    = "#0047AB",

      "Uganda"   = "#F564E3"

    ),

    na.value = "grey50"

  ) +

  scale_shape_manual(

    values = c(

      "Human"   = 16,

      "Pig"     = 17,

      "Unknown" = 15

    )

  ) +

  theme_tree2() +

  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +

  coord_cartesian(xlim = c(0, 0.17), clip = "off") +

  labs(

    title = "Maximum likelihood tree",

    colour = "Region",

    shape = "Host"

  ) +

  theme(

    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),

    legend.position = c(1.02, 0.55),

    legend.justification = c(0, 0.5),

    legend.background = element_blank(),

    legend.key = element_blank(),

    plot.margin = margin(10, 180, 10, 10)

  )

pML

ggsave(

  "Ascaris_65_ML_tree_ASC.png",

  plot = pML,

  width = 16,

  height = 12,

  dpi = 600

)
```

#### ADMIXTURE

```text
VII) ADMIXTURE

# Purpose: infer ancestry components from the balanced LD-pruned dataset using ADMIXTURE across K = 1–10

# with replicated runs, then plot selected K values.

# Input: balanced n = 65 LD-pruned PLINK files.

# Output: .Q ancestry-proportion files, .P allele-frequency files, and ADMIXTURE log files.

#!/bin/bash

#SBATCH --job-name=admix_100k_reps

#SBATCH --partition=high_mem

#SBATCH --cpus-per-task=20

#SBATCH --mem=100G

#SBATCH --time=72:00:00

#SBATCH --output=admix_100k_reps.%j.out

#SBATCH --error=admix_100k_reps.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh

conda activate NJ_env

PLINKDIR=/parallel_scratch/lw01165/05_analysis/05e_plink_65

OUTDIR=/parallel_scratch/lw01165/05_analysis/05h_admixture_65/100K_SNPs

CHROM_FILE=/parallel_scratch/lw01165/chrom_order.txt

# chrom_order.txt defines the autosomes to retain and preserves the intended chromosome order.

mkdir -p ${OUTDIR}

cd ${OUTDIR}

echo "Checking for existing 100k SNP dataset..."

date

if [ ! -f Ascaris.65.pruned.100kSNPs.bed ]; then

    echo "100k dataset not found — creating it now"

    cp ${PLINKDIR}/Ascaris.65.pruned.bed .

    cp ${PLINKDIR}/Ascaris.65.pruned.bim .

    cp ${PLINKDIR}/Ascaris.65.pruned.fam .

    # Randomly thin the LD-pruned dataset to 100,000 SNPs for ADMIXTURE.

    plink2 \
      --bfile Ascaris.65.pruned \
      --thin-count 100000 \
      --make-bed \
      --out Ascaris.65.pruned.100kSNPs

    cp Ascaris.65.pruned.100kSNPs.bim Ascaris.65.pruned.100kSNPs.original.bim

    # ADMIXTURE expects numeric chromosome codes. This rewrites chromosome names in the BIM file using

    # the order defined in chrom_order.txt while preserving variant IDs, genetic distance, position and

    # alleles.

    awk 'BEGIN{FS="[[:space:]]+"; OFS="\\t"} NR==FNR {map[$1]=NR; next} {$1=map[$1]; print $1,$2,$3,$4,$5,$6}' \
      ${CHROM_FILE} \
      Ascaris.65.pruned.100kSNPs.original.bim > Ascaris.65.pruned.100kSNPs.bim

else

    echo "100k SNP dataset already exists — skipping creation"

fi

echo "Starting replicated ADMIXTURE runs..."

date

# Run K = 1 to 10, with 10 replicates per K. Replicates are useful because ADMIXTURE can converge on

# different local optima.

for K in {1..10}

do

    for REP in {1..10}

    do

        echo "Running ADMIXTURE K=${K}, replicate=${REP}"

        # --cv writes cross-validation error to the log, which should be parsed later to choose the best K.

        admixture --cv -j16 Ascaris.65.pruned.100kSNPs.bed ${K} \
          | tee log100kSNPs_K${K}_rep${REP}.out

        mv Ascaris.65.pruned.100kSNPs.${K}.Q Ascaris.65.pruned.100kSNPs.K${K}.rep${REP}.Q

        mv Ascaris.65.pruned.100kSNPs.${K}.P Ascaris.65.pruned.100kSNPs.K${K}.rep${REP}.P

    done

done

echo "ADMIXTURE analysis finished"

date

# Analysis in RStudio

# Purpose: read selected ADMIXTURE Q matrices, join them to metadata, and produce stacked ancestry plots

# with country and host annotation bars.

setwd("S:/Research/Lauren_PopGen/ADMIXTURE")

library(tidyverse)

library(patchwork)

# ============================================================

# ADMIXTURE CV error summary

# ============================================================

# ADMIXTURE was run with 10 replicates per K. CV errors were

# extracted to cv_values_100k.txt. These values are summarised

# here to identify the best-supported K and to check whether

# replicate runs produced stable CV errors.

cv <- read.table(

  "cv_values_100k.txt",

  header = FALSE,

  col.names = c("K", "replicate", "cv_error")

)

cv_summary <- cv %>%

  group_by(K) %>%

  summarise(

    n_replicates = n(),

    mean_cv = mean(cv_error),

    sd_cv = sd(cv_error),

    min_cv = min(cv_error),

    max_cv = max(cv_error),

    .groups = "drop"

  ) %>%

  arrange(K)

print(cv_summary)

write.csv(

  cv_summary,

  "ADMIXTURE_CV_summary_100k.csv",

  row.names = FALSE

)

p_cv <- ggplot(cv, aes(x = K, y = cv_error)) +

  geom_point(size = 2, alpha = 0.7) +

  stat_summary(fun = mean, geom = "line") +

  scale_x_continuous(breaks = sort(unique(cv$K))) +

  theme_bw() +

  labs(

    x = "K",

    y = "Cross-validation error",

    title = "ADMIXTURE cross-validation error"

  )

p_cv

ggsave(

  "ADMIXTURE_CV_error_100k.png",

  p_cv,

  width = 7,

  height = 5,

  dpi = 600

)

# Identify best-supported K

best_K <- cv_summary %>%

  filter(mean_cv == min(mean_cv))

print(best_K)

# In this dataset, the lowest mean CV error supports the best-fitting K.

# ADMIXTURE cluster colours

cluster_cols <- c(

  Cluster1 = "#03045e",

  Cluster2 = "#0077b6",

  Cluster3 = "#00b4d8",

  Cluster4 = "#90e0ef"

)

# country colours for annotation bar

country_cols <- c(

  Ethiopia    = "#8c510a",

  Kenya       = "#bf812d",

  Uganda      = "#dfc27d",

  India       = "#9FE2BF",

  Philippines = "#40826D",

  Thailand    = "#004643",

  Germany     = "#470047",

  Ireland     = "#E0B0FF",

  Hungary     = "#A300A3"

)

# host colours for annotation bar

host_cols <- c(

  Human = "black",

  Pig   = "white"

)

# Function to combine ADMIXTURE Q proportions with PLINK sample order and metadata.

prepare_admixture_data <- function(q_file, fam_file, meta_file) {

  q <- read.table(q_file, header = FALSE)

  fam <- read.table(fam_file, header = FALSE)

  meta <- read.csv(meta_file, stringsAsFactors = FALSE)

  

  meta$sample_id <- trimws(meta$sample_id)

  meta$country   <- trimws(meta$country)

  meta$region    <- trimws(meta$region)

  meta$host      <- trimws(meta$host)

  

  fam$sample_id <- trimws(fam$V2)

  q$sample_id <- fam$sample_id

  

  k <- ncol(q) - 1

  colnames(q) <- c(paste0("Cluster", 1:k), "sample_id")

  

  admix <- left_join(q, meta, by = "sample_id") %>%

    arrange(region, country, sample_id)

  

  admix$sample_id <- factor(admix$sample_id, levels = admix$sample_id)

  

  admix

}

make_admixture_plot <- function(admix, K_label, show_legend = FALSE) {

  admix_long <- admix %>%

    pivot_longer(

      cols = starts_with("Cluster"),

      names_to = "Cluster",

      values_to = "Ancestry"

    )

  

  p <- ggplot(admix_long, aes(x = sample_id, y = Ancestry, fill = Cluster)) +

    geom_col(width = 1) +

    scale_fill_manual(values = cluster_cols, drop = FALSE) +

    theme_bw() +

    labs(

      title = K_label,

      x = NULL,

      y = "Ancestry proportion"

    ) +

    theme(

      axis.text.x = element_blank(),

      axis.ticks.x = element_blank(),

      strip.background = element_rect(fill = "grey90"),

      strip.text = element_text(size = 10),

      legend.position = if (show_legend) "right" else "none",

      panel.spacing.x = unit(0.3, "lines")

    )

  

  p <- p + facet_grid(\~region, scales = "free_x", space = "free_x")

  return(p)

}

make_country_bar <- function(admix) {

  p <- ggplot(admix, aes(x = sample_id, y = 1, fill = country)) +

    geom_tile() +

    scale_fill_manual(

      name = "Country",

      values = country_cols,

      drop = FALSE

    ) +

    guides(fill = guide_legend(ncol = 1)) +

    theme_void() +

    theme(

      legend.position = "right",

      legend.title = element_text(size = 10),

      legend.text = element_text(size = 9),

      panel.spacing = unit(0, "lines"),

      panel.spacing.x = unit(0.3, "lines"),

      strip.text = element_blank(),

      strip.background = element_blank()

    ) +

    facet_grid(\~region, scales = "free_x", space = "free_x")

  

  return(p)

}

make_host_bar <- function(admix) {

  p <- ggplot(admix, aes(x = sample_id, y = 1, fill = host)) +

    geom_tile(color = "black", linewidth = 0.25) +

    scale_fill_manual(

      name = "Host",

      values = host_cols,

      drop = FALSE

    ) +

    guides(fill = guide_legend(nrow = 1)) +

    theme_void() +

    theme(

      legend.position = "right",

      legend.title = element_text(size = 10),

      legend.text = element_text(size = 9),

      panel.spacing = unit(0, "lines"),

      panel.spacing.x = unit(0.3, "lines"),

      strip.text = element_blank(),

      strip.background = element_blank(),

      plot.margin = margin(5.5, 20, 5.5, 5.5)

    ) +

    facet_grid(\~region, scales = "free_x", space = "free_x")

  

  return(p)

}

# input files

fam_file <- "Ascaris.65.pruned.100kSNPs.fam"

meta_file <- "S:/Research/Lauren_PopGen/sample_QC_list.csv"

# prepare datasets

# ADMIXTURE was run with 10 replicates per K. CV errors were extracted to cv_values_100k.txt

# and were identical or near-identical across replicates for each K, indicating stable

# convergence. Replicate 1 was therefore used for plotting K = 2, K = 3 and K = 4.

admix_k2 <- prepare_admixture_data(

  "Ascaris.65.pruned.100kSNPs.K2.rep1.Q",

  fam_file,

  meta_file

)

admix_k3 <- prepare_admixture_data(

  "Ascaris.65.pruned.100kSNPs.K3.rep1.Q",

  fam_file,

  meta_file

)

admix_k4 <- prepare_admixture_data(

  "Ascaris.65.pruned.100kSNPs.K4.rep1.Q",

  fam_file,

  meta_file

)

# make plots

p2_main <- make_admixture_plot(admix_k2, "K = 2", show_legend = FALSE)

p2_bar  <- make_country_bar(admix_k2)

p2_host <- make_host_bar(admix_k2)

p3      <- make_admixture_plot(admix_k3, "K = 3", show_legend = FALSE)

p4      <- make_admixture_plot(admix_k4, "K = 4", show_legend = TRUE)

# assemble everything

combined_plot <- wrap_plots(

  p2_main, p2_bar, p2_host, p3, p4,

  ncol = 1,

  heights = c(10, 0.7, 0.7, 10, 10)

)

 combined_plot

ggsave(

  "ADMIXTURE_K2_K3_K4_combined_clean.png",

  combined_plot,

  width = 14,

  height = 12,

  dpi = 600

)
```



---

## Outputs and QC checks

PLINK pruned data, n=65 VCF, PCA, IBS, PERMANOVA, NJ/ML trees, ADMIXTURE. Consult recorded logs and intermediate tables to distinguish expected settings from actual execution success.

## Provenance and limitations

- Original population structure used LD-pruning r²=0.15; balanced n=65 construction inserted.
- Sample identifiers and host/country metadata were manually curated in `sample_QC_list.csv`; it is not generated by an algorithm in these files.
- Commands, versions, seeds and paths are retained where present; annotations are not evidence that all jobs have been freshly rerun.
- For the discovery and correction of the PIXY all-sites input issue, see [`99_PIXY_allsites_correction.md`](99_PIXY_allsites_correction.md).
