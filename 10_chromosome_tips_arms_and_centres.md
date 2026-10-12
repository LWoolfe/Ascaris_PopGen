# 10 — Chromosomal tip, arm and centre analysis

## Purpose

Chromosomal tip, arm and centre analysis. This file archives the commands and analysis notes from the thesis/paper workflow.

## Inputs

PIXY windows, 25-kb coverage and mappability, chromosome lengths. Source filenames and absolute HPC paths are retained in the original script blocks below.

## Workflow

See the original annotated analysis sequence below. Distinct script blocks are independent jobs/R sessions, **not one executable script**.

## Annotated scripts

> **Historical PIXY-derived analysis:** This uses historical diversity/differentiation outputs, not the corrected PIXY VCF.

```r
# ============================================================
# ASCARIS CHROMOSOME TIP QC ANALYSIS
# ============================================================
#
# Purpose:
# Test whether the M-shaped chromosome-wide pattern in
# nucleotide diversity (pi) and dXY is associated with:
#
#   - sequencing depth
#   - fraction of the window covered
#   - unique mappability
#   - number/proportion of sites contributing to PIXY estimates
#
# Diversity data = 10-kb windows
# Depth/mappability = 25-kb windows
#
# Therefore all datasets are converted to common 50-kb windows.
#
# ============================================================


# ------------------------------------------------------------
# 0. PACKAGES
# ------------------------------------------------------------

required_packages <- c(
  "data.table",
  "ggplot2"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  stop(
    "Please install these packages first: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(data.table)
library(ggplot2)


# ============================================================
# 1. FILE PATHS
# ============================================================

# Folder containing the copied Eureka files
QC_DIR <- "S:/Research/Lauren_PopGen/pi TajimasD/human_pig"


# Existing nucleotide diversity file
PI_FILE <- paste0(
  "S:/Research/Lauren_PopGen/pi TajimasD/",
  "all_146_pi.txt"
)


# Human vs pig dXY file
DXY_FILE <- file.path(
  QC_DIR,
  "human_vs_pig_dxy.txt"
)


# 25-kb depth / coverage file copied from Eureka
COV_FILE <- file.path(
  QC_DIR,
  "New_25kb_cov_stats.tsv"
)


# 25-kb mappability file copied from Eureka
MAP_FILE <- file.path(
  QC_DIR,
  "Ascaris.autosomes19.25kb.mappability.tsv"
)


# Final 146 sample IDs
SAMPLE_FILE <- file.path(
  QC_DIR,
  "final146.samples.txt"
)


# Reference genome chromosome lengths
GENOME_FILE <- file.path(
  QC_DIR,
  "Ascaris_primary_scaffolds_plus_X.genome.txt"
)


# Save results locally to Downloads rather than the S: drive

OUT_DIR <- file.path(
  Sys.getenv("USERPROFILE"),
  "Downloads",
  "Ascaris_tip_QC_results"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

cat(
  "\nOutput directory:\n",
  OUT_DIR,
  "\n"
)

if (!dir.exists(OUT_DIR)) {
  stop("Output directory could not be created.")
}

# ============================================================
# 2. CHECK ALL FILES EXIST
# ============================================================

files_to_check <- c(
  PI_FILE,
  DXY_FILE,
  COV_FILE,
  MAP_FILE,
  SAMPLE_FILE,
  GENOME_FILE
)

file_check <- data.table(
  file = files_to_check,
  exists = file.exists(files_to_check)
)

print(file_check)

if (any(!file_check$exists)) {
  stop(
    "\nOne or more required files cannot be found.\n",
    "Check the paths printed above before continuing."
  )
}


# ============================================================
# 3. ANALYSIS PARAMETERS
# ============================================================

# Common window size
BIN_SIZE <- 50000


# Width used for the normalised chromosome profile
# 0.025 = 2.5% of chromosome length
PROFILE_BIN <- 0.025


# Test three alternative definitions of chromosome tips
TIP_CUTOFFS <- c(
  0.025,
  0.05,
  0.10
)


# Adjacent chromosome arms extend to 25%
# of chromosome length
ARM_BOUNDARY <- 0.25


# QC thresholds

# 90% uniquely mappable
MAP_THRESHOLD <- 0.90

# 90% of bases covered
BREADTH_THRESHOLD <- 0.90

# Require 80% of genomic positions to contribute
# to PIXY statistic
SITE_THRESHOLD <- 0.80


# Rather than imposing an arbitrary depth threshold,
# define low depth as the bottom 10% of the genome-wide
# median-depth distribution.
DEPTH_QUANTILE <- 0.10


# ============================================================
# 4. HELPER FUNCTIONS
# ============================================================

safe_wmean <- function(x, w) {
  
  good <- (
    !is.na(x) &
      !is.na(w) &
      w > 0
  )
  
  if (!any(good)) {
    return(NA_real_)
  }
  
  weighted.mean(
    x[good],
    w[good]
  )
}


safe_spearman <- function(x, y) {
  
  good <- complete.cases(x, y)
  
  if (sum(good) < 4) {
    
    return(
      data.table(
        n = sum(good),
        rho = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  test <- suppressWarnings(
    cor.test(
      x[good],
      y[good],
      method = "spearman",
      exact = FALSE
    )
  )
  
  data.table(
    n = sum(good),
    rho = unname(test$estimate),
    p_value = test$p.value
  )
}


# ============================================================
# 5. IDENTIFY THE 19 AUTOSOMES
# ============================================================
#
# We do NOT need autosomes19.txt.
#
# The mappability dataset is already explicitly restricted
# to the 19 autosomes, so use its chromosome names.
#
# ============================================================

cat(
  "\n========================================\n",
  "IDENTIFYING AUTOSOMES\n",
  "========================================\n"
)


map_chroms <- fread(
  MAP_FILE,
  select = "chrom"
)

autosomes <- unique(
  map_chroms$chrom
)

rm(map_chroms)


cat(
  "\nAutosomes identified from mappability file:",
  length(autosomes),
  "\n\n"
)

print(autosomes)


if (length(autosomes) != 19) {
  stop(
    "Expected 19 autosomes in the mappability file, but found ",
    length(autosomes),
    "."
  )
}


# ============================================================
# 6. READ REFERENCE CHROMOSOME LENGTHS
# ============================================================

cat(
  "\n========================================\n",
  "READING CHROMOSOME LENGTHS\n",
  "========================================\n"
)


genome <- fread(
  GENOME_FILE,
  header = FALSE
)


if (ncol(genome) < 2) {
  stop(
    "Genome file does not contain at least two columns."
  )
}


# Only use first two columns:
# chromosome and chromosome length

genome <- genome[, 1:2]


setnames(
  genome,
  c(
    "chromosome",
    "chrom_length"
  )
)


# Keep only the 19 autosomes

genome <- genome[
  chromosome %chin% autosomes
]


cat(
  "\nAutosome lengths found:",
  nrow(genome),
  "\n\n"
)


print(genome)


if (nrow(genome) != 19) {
  stop(
    "Expected reference lengths for 19 autosomes but found ",
    nrow(genome),
    "."
  )
}


# ============================================================
# 7. READ FINAL 146 SAMPLE LIST
# ============================================================

cat(
  "\n========================================\n",
  "READING FINAL SAMPLE LIST\n",
  "========================================\n"
)


samples146 <- scan(
  SAMPLE_FILE,
  what = "character",
  quiet = TRUE
)


samples146 <- unique(
  samples146
)


cat(
  "\nFinal sample IDs:",
  length(samples146),
  "\n"
)


if (length(samples146) != 146) {
  warning(
    "Expected 146 sample IDs but found ",
    length(samples146),
    "."
  )
}


# ============================================================
# PART A
# BUILD 50-KB DEPTH + COVERAGE TABLE
# ============================================================

cat(
  "\n========================================\n",
  "READING DEPTH / COVERAGE DATA\n",
  "========================================\n"
)


cov <- fread(
  COV_FILE,
  select = c(
    "chrom",
    "start",
    "end",
    "bases_covered",
    "window_length",
    "fraction_covered",
    "mean_depth",
    "sample"
  )
)


cat(
  "\nSamples in raw coverage file:",
  uniqueN(cov$sample),
  "\n"
)


# ------------------------------------------------------------
# Keep the 19 autosomes
# and exact final 146 samples
# ------------------------------------------------------------

cov <- cov[
  chrom %chin% autosomes &
    sample %chin% samples146
]


cat(
  "Samples retained after filtering:",
  uniqueN(cov$sample),
  "\n"
)


# Check that all 146 final samples are present

missing_samples <- setdiff(
  samples146,
  unique(cov$sample)
)


if (length(missing_samples) > 0) {
  
  stop(
    "\nThe following final samples were not found in ",
    "New_25kb_cov_stats.tsv:\n\n",
    paste(
      missing_samples,
      collapse = "\n"
    )
  )
}


# ------------------------------------------------------------
# Assign each 25-kb window to an exact 50-kb bin
# ------------------------------------------------------------

cov[
  ,
  bin50_start :=
    floor(
      start / BIN_SIZE
    ) *
    BIN_SIZE
]


# ------------------------------------------------------------
# First combine the two 25-kb windows WITHIN each sample
# ------------------------------------------------------------

cov50_sample <- cov[
  ,
  .(
    
    mean_depth_50kb =
      safe_wmean(
        mean_depth,
        window_length
      ),
    
    fraction_covered_50kb =
      sum(
        bases_covered,
        na.rm = TRUE
      ) /
      sum(
        window_length,
        na.rm = TRUE
      ),
    
    contributing_length =
      sum(
        window_length,
        na.rm = TRUE
      )
    
  ),
  by = .(
    sample,
    chrom,
    bin50_start
  )
]


# We no longer need huge raw table

rm(cov)

gc()


# ------------------------------------------------------------
# Now summarise across the final 146 samples
# ------------------------------------------------------------

cov50 <- cov50_sample[
  ,
  .(
    
    median_depth =
      median(
        mean_depth_50kb,
        na.rm = TRUE
      ),
    
    mean_depth =
      mean(
        mean_depth_50kb,
        na.rm = TRUE
      ),
    
    median_fraction_covered =
      median(
        fraction_covered_50kb,
        na.rm = TRUE
      ),
    
    mean_fraction_covered =
      mean(
        fraction_covered_50kb,
        na.rm = TRUE
      ),
    
    n_samples = .N
    
  ),
  by = .(
    chrom,
    bin50_start
  )
]


rm(cov50_sample)

gc()


cat(
  "\nNumber of samples contributing to each 50-kb window:\n"
)

print(
  summary(
    cov50$n_samples
  )
)


# ============================================================
# PART B
# BUILD 50-KB MAPPABILITY TABLE
# ============================================================

cat(
  "\n========================================\n",
  "READING MAPPABILITY DATA\n",
  "========================================\n"
)


map <- fread(
  MAP_FILE
)


map <- map[
  chrom %chin% autosomes
]


# Assign each 25-kb window to its 50-kb bin

map[
  ,
  bin50_start :=
    floor(
      start / BIN_SIZE
    ) *
    BIN_SIZE
]


# Sum unique bases across the two 25-kb windows

map50 <- map[
  ,
  .(
    
    unique_bases =
      sum(
        unique_bases,
        na.rm = TRUE
      ),
    
    mappability_window_length =
      sum(
        window_length,
        na.rm = TRUE
      )
    
  ),
  by = .(
    chrom,
    bin50_start
  )
]


map50[
  ,
  unique_fraction :=
    unique_bases /
    mappability_window_length
]


rm(map)

gc()


# ============================================================
# PART C
# CREATE MASTER TECHNICAL QC TABLE
# ============================================================

qc50 <- merge(
  cov50,
  map50,
  by = c(
    "chrom",
    "bin50_start"
  ),
  all = TRUE
)


rm(cov50)
rm(map50)

gc()


# Rename chromosome column

setnames(
  qc50,
  "chrom",
  "chromosome"
)


# Add true reference chromosome length

qc50 <- merge(
  qc50,
  genome,
  by = "chromosome",
  all.x = TRUE
)


# ------------------------------------------------------------
# Correct the final bin if a chromosome does not end exactly
# on a 50-kb boundary
# ------------------------------------------------------------

qc50[
  ,
  bin50_end :=
    pmin(
      bin50_start + BIN_SIZE,
      chrom_length
    )
]


qc50[
  ,
  bin_length :=
    bin50_end -
    bin50_start
]


# Remove impossible / zero-length bins

qc50 <- qc50[
  bin_length > 0
]


# Midpoint of each bin

qc50[
  ,
  midpoint :=
    bin50_start +
    bin_length / 2
]


# Position from 0 to 1 along chromosome

qc50[
  ,
  relative_position :=
    midpoint /
    chrom_length
]


# Distance to nearest chromosome end in bases

qc50[
  ,
  distance_from_end :=
    pmin(
      midpoint,
      chrom_length - midpoint
    )
]


setorder(
  qc50,
  chromosome,
  bin50_start
)


# ------------------------------------------------------------
# Determine empirical low-depth threshold
# ------------------------------------------------------------

DEPTH_THRESHOLD <- as.numeric(
  quantile(
    qc50$median_depth,
    probs = DEPTH_QUANTILE,
    na.rm = TRUE
  )
)


cat(
  "\n========================================\n",
  "QC THRESHOLDS\n",
  "========================================\n",
  "Unique mappability >= ",
  MAP_THRESHOLD,
  "\n",
  "Median breadth >= ",
  BREADTH_THRESHOLD,
  "\n",
  "PIXY site fraction >= ",
  SITE_THRESHOLD,
  "\n",
  "Bottom 10% median-depth threshold = ",
  round(
    DEPTH_THRESHOLD,
    3
  ),
  "\n",
  sep = ""
)


# Save technical QC table

fwrite(
  qc50,
  file.path(
    OUT_DIR,
    "QC_50kb_master.tsv"
  ),
  sep = "\t"
)


# ============================================================
# PART D
# AGGREGATE HUMAN VS PIG DXY TO 50 KB
# ============================================================

cat(
  "\n========================================\n",
  "READING HUMAN VS PIG DXY\n",
  "========================================\n"
)


dxy <- fread(
  DXY_FILE
)


# Keep only autosomes and intended comparison

dxy <- dxy[
  chromosome %chin% autosomes &
    pop1 == "Human" &
    pop2 == "Pig"
]


cat(
  "\n10-kb dXY rows retained:",
  nrow(dxy),
  "\n"
)


# ------------------------------------------------------------
# PIXY coordinates are 1-based.
#
# Windows:
#       1-10000
#   10001-20000
#
# Convert to 0-based 50-kb bins.
# ------------------------------------------------------------

dxy[
  ,
  bin50_start :=
    floor(
      (window_pos_1 - 1) /
        BIN_SIZE
    ) *
    BIN_SIZE
]


# ------------------------------------------------------------
# Aggregate correctly by summing PIXY raw components.
#
# avg_dxy for the new 50-kb window is:
#
# sum(count_diffs) / sum(count_comparisons)
#
# ------------------------------------------------------------

dxy50 <- dxy[
  ,
  .(
    
    dxy_diffs =
      sum(
        count_diffs,
        na.rm = TRUE
      ),
    
    dxy_comparisons =
      sum(
        count_comparisons,
        na.rm = TRUE
      ),
    
    dxy_sites =
      sum(
        no_sites,
        na.rm = TRUE
      )
    
  ),
  by = .(
    chromosome,
    bin50_start
  )
]


dxy50[
  ,
  avg_dxy_50kb :=
    fifelse(
      dxy_comparisons > 0,
      dxy_diffs /
        dxy_comparisons,
      NA_real_
    )
]


rm(dxy)

gc()


# ------------------------------------------------------------
# Add technical QC information
# ------------------------------------------------------------

master_dxy <- merge(
  qc50,
  dxy50,
  by = c(
    "chromosome",
    "bin50_start"
  ),
  all.x = TRUE
)


rm(dxy50)

gc()


# Fraction of genomic positions contributing to PIXY dXY

master_dxy[
  ,
  dxy_site_fraction :=
    dxy_sites /
    bin_length
]


# ------------------------------------------------------------
# QC flags
# ------------------------------------------------------------

master_dxy[
  ,
  high_map :=
    !is.na(unique_fraction) &
    unique_fraction >= MAP_THRESHOLD
]


master_dxy[
  ,
  high_depth :=
    !is.na(median_depth) &
    median_depth >= DEPTH_THRESHOLD
]


master_dxy[
  ,
  high_breadth :=
    !is.na(median_fraction_covered) &
    median_fraction_covered >=
    BREADTH_THRESHOLD
]


master_dxy[
  ,
  high_dxy_sites :=
    !is.na(dxy_site_fraction) &
    dxy_site_fraction >=
    SITE_THRESHOLD
]


master_dxy[
  ,
  strict_qc :=
    high_map &
    high_depth &
    high_breadth &
    high_dxy_sites
]


fwrite(
  master_dxy,
  file.path(
    OUT_DIR,
    "DXY_50kb_master.tsv"
  ),
  sep = "\t"
)


# ============================================================
# PART E
# AGGREGATE NUCLEOTIDE DIVERSITY TO 50 KB
# ============================================================

cat(
  "\n========================================\n",
  "READING NUCLEOTIDE DIVERSITY\n",
  "========================================\n"
)


pi <- fread(
  PI_FILE
)


pi <- pi[
  chromosome %chin% autosomes
]


cat(
  "\nPopulations in nucleotide diversity file:\n"
)

print(
  unique(pi$pop)
)


cat(
  "\n10-kb pi rows retained:",
  nrow(pi),
  "\n"
)


# Assign 10-kb windows to 50-kb bins

pi[
  ,
  bin50_start :=
    floor(
      (window_pos_1 - 1) /
        BIN_SIZE
    ) *
    BIN_SIZE
]


# ------------------------------------------------------------
# Reconstruct combined pi from PIXY raw components:
#
# sum(count_diffs) / sum(count_comparisons)
# ------------------------------------------------------------

pi50 <- pi[
  ,
  .(
    
    pi_diffs =
      sum(
        count_diffs,
        na.rm = TRUE
      ),
    
    pi_comparisons =
      sum(
        count_comparisons,
        na.rm = TRUE
      ),
    
    pi_sites =
      sum(
        no_sites,
        na.rm = TRUE
      )
    
  ),
  by = .(
    pop,
    chromosome,
    bin50_start
  )
]


pi50[
  ,
  avg_pi_50kb :=
    fifelse(
      pi_comparisons > 0,
      pi_diffs /
        pi_comparisons,
      NA_real_
    )
]


rm(pi)

gc()


# Add technical QC information

master_pi <- merge(
  pi50,
  qc50,
  by = c(
    "chromosome",
    "bin50_start"
  ),
  all.x = TRUE
)


rm(pi50)

gc()


# Fraction of positions represented in PIXY

master_pi[
  ,
  pi_site_fraction :=
    pi_sites /
    bin_length
]


# QC flags

master_pi[
  ,
  high_map :=
    !is.na(unique_fraction) &
    unique_fraction >=
    MAP_THRESHOLD
]


master_pi[
  ,
  high_depth :=
    !is.na(median_depth) &
    median_depth >=
    DEPTH_THRESHOLD
]


master_pi[
  ,
  high_breadth :=
    !is.na(median_fraction_covered) &
    median_fraction_covered >=
    BREADTH_THRESHOLD
]


master_pi[
  ,
  high_pi_sites :=
    !is.na(pi_site_fraction) &
    pi_site_fraction >=
    SITE_THRESHOLD
]


master_pi[
  ,
  strict_qc :=
    high_map &
    high_depth &
    high_breadth &
    high_pi_sites
]


fwrite(
  master_pi,
  file.path(
    OUT_DIR,
    "PI_50kb_master.tsv"
  ),
  sep = "\t"
)


# ============================================================
# PART F
# BASIC QC SUMMARY
# ============================================================

cat(
  "\n\n========================================\n",
  "50-KB QC SUMMARY\n",
  "========================================\n"
)


cat(
  "\nUnique mappability:\n"
)

print(
  summary(
    qc50$unique_fraction
  )
)


cat(
  "\nMedian depth:\n"
)

print(
  summary(
    qc50$median_depth
  )
)


cat(
  "\nMedian coverage breadth:\n"
)

print(
  summary(
    qc50$median_fraction_covered
  )
)


cat(
  "\ndXY site fraction:\n"
)

print(
  summary(
    master_dxy$dxy_site_fraction
  )
)


cat(
  "\nDXY windows total:",
  nrow(master_dxy),
  "\n"
)


cat(
  "DXY windows >=90% mappability:",
  sum(
    master_dxy$high_map,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "DXY windows passing strict QC:",
  sum(
    master_dxy$strict_qc,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# PART G
# NORMALISED CHROMOSOME POSITION
# ============================================================

# ------------------------------------------------------------
# Bin relative chromosome position in 2.5% increments
# ------------------------------------------------------------

master_dxy[
  ,
  rel_bin :=
    floor(
      relative_position /
        PROFILE_BIN
    ) *
    PROFILE_BIN +
    PROFILE_BIN / 2
]


master_pi[
  ,
  rel_bin :=
    floor(
      relative_position /
        PROFILE_BIN
    ) *
    PROFILE_BIN +
    PROFILE_BIN / 2
]


qc50[
  ,
  rel_bin :=
    floor(
      relative_position /
        PROFILE_BIN
    ) *
    PROFILE_BIN +
    PROFILE_BIN / 2
]


# ============================================================
# PART H
# DXY NORMALISED CHROMOSOME PROFILES
# ============================================================

make_dxy_profile <- function(
    dat,
    label
) {
  
  # Calculate one estimate for each chromosome
  # within each relative-position bin
  
  chr_profile <- dat[
    !is.na(dxy_comparisons) &
      dxy_comparisons > 0,
    .(
      
      value =
        sum(
          dxy_diffs,
          na.rm = TRUE
        ) /
        sum(
          dxy_comparisons,
          na.rm = TRUE
        )
      
    ),
    by = .(
      chromosome,
      rel_bin
    )
  ]
  
  
  # Then summarise ACROSS chromosomes
  
  out <- chr_profile[
    ,
    .(
      
      median =
        median(
          value,
          na.rm = TRUE
        ),
      
      q25 =
        quantile(
          value,
          0.25,
          na.rm = TRUE
        ),
      
      q75 =
        quantile(
          value,
          0.75,
          na.rm = TRUE
        ),
      
      n_chr =
        uniqueN(chromosome)
      
    ),
    by = rel_bin
  ]
  
  
  out[
    ,
    QC := label
  ]
  
  
  out
}


# Produce profiles for:
#
# 1. all windows
# 2. >=90% uniquely mappable
# 3. strict QC

dxy_profile <- rbindlist(
  list(
    
    make_dxy_profile(
      master_dxy,
      "All windows"
    ),
    
    make_dxy_profile(
      master_dxy[
        high_map == TRUE
      ],
      "Mappability >= 90%"
    ),
    
    make_dxy_profile(
      master_dxy[
        strict_qc == TRUE
      ],
      "Strict QC"
    )
    
  )
)


fwrite(
  dxy_profile,
  file.path(
    OUT_DIR,
    "DXY_relative_position_profile.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# Plot dXY
# ------------------------------------------------------------

p_dxy <- ggplot(
  dxy_profile,
  aes(
    x = rel_bin,
    y = median
  )
) +
  
  geom_ribbon(
    aes(
      ymin = q25,
      ymax = q75
    ),
    alpha = 0.20
  ) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  facet_wrap(
    ~ QC,
    ncol = 1
  ) +
  
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(
      0,
      1,
      by = 0.1
    )
  ) +
  
  labs(
    x = "Relative chromosome position",
    y = expression(d[XY]),
    title = "Human vs pig dXY across normalised chromosomes"
  ) +
  
  theme_classic()


ggsave(
  file.path(
    OUT_DIR,
    "01_DXY_relative_position_QC.png"
  ),
  p_dxy,
  width = 8,
  height = 9,
  dpi = 300
)


# ============================================================
# PART I
# PI NORMALISED CHROMOSOME PROFILES
# ============================================================

make_pi_profile <- function(
    dat,
    label
) {
  
  # One value per population/chromosome/position bin
  
  chr_profile <- dat[
    !is.na(pi_comparisons) &
      pi_comparisons > 0,
    .(
      
      value =
        sum(
          pi_diffs,
          na.rm = TRUE
        ) /
        sum(
          pi_comparisons,
          na.rm = TRUE
        )
      
    ),
    by = .(
      pop,
      chromosome,
      rel_bin
    )
  ]
  
  
  # Median across chromosomes
  
  out <- chr_profile[
    ,
    .(
      
      median =
        median(
          value,
          na.rm = TRUE
        ),
      
      q25 =
        quantile(
          value,
          0.25,
          na.rm = TRUE
        ),
      
      q75 =
        quantile(
          value,
          0.75,
          na.rm = TRUE
        ),
      
      n_chr =
        uniqueN(chromosome)
      
    ),
    by = .(
      pop,
      rel_bin
    )
  ]
  
  
  out[
    ,
    QC := label
  ]
  
  
  out
}


pi_profile <- rbindlist(
  list(
    
    make_pi_profile(
      master_pi,
      "All windows"
    ),
    
    make_pi_profile(
      master_pi[
        high_map == TRUE
      ],
      "Mappability >= 90%"
    ),
    
    make_pi_profile(
      master_pi[
        strict_qc == TRUE
      ],
      "Strict QC"
    )
    
  )
)


fwrite(
  pi_profile,
  file.path(
    OUT_DIR,
    "PI_relative_position_profile.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# Plot pi
# ------------------------------------------------------------

p_pi <- ggplot(
  pi_profile,
  aes(
    x = rel_bin,
    y = median,
    group = pop,
    linetype = pop
  )
) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  facet_wrap(
    ~ QC,
    ncol = 1
  ) +
  
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(
      0,
      1,
      by = 0.1
    )
  ) +
  
  labs(
    x = "Relative chromosome position",
    y = expression(pi),
    linetype = "Population",
    title = "Nucleotide diversity across normalised chromosomes"
  ) +
  
  theme_classic()


ggsave(
  file.path(
    OUT_DIR,
    "02_PI_relative_position_QC.png"
  ),
  p_pi,
  width = 8,
  height = 9,
  dpi = 300
)


# ============================================================
# PART J
# TECHNICAL QC PROFILES
# ============================================================

# Turn QC metrics into long format

qc_long <- melt(
  qc50,
  id.vars = c(
    "chromosome",
    "rel_bin"
  ),
  measure.vars = c(
    "median_depth",
    "median_fraction_covered",
    "unique_fraction"
  ),
  variable.name = "metric",
  value.name = "value"
)


# Add PIXY dXY site fraction

dxy_sites_qc <- master_dxy[
  ,
  .(
    chromosome,
    rel_bin,
    metric = "dxy_site_fraction",
    value = dxy_site_fraction
  )
]


qc_long <- rbindlist(
  list(
    qc_long,
    dxy_sites_qc
  ),
  use.names = TRUE
)


# First get one value per chromosome

qc_chr_profile <- qc_long[
  !is.na(value),
  .(
    value =
      median(
        value,
        na.rm = TRUE
      )
  ),
  by = .(
    chromosome,
    rel_bin,
    metric
  )
]


# Then median and IQR across chromosomes

qc_profile <- qc_chr_profile[
  ,
  .(
    
    median =
      median(
        value,
        na.rm = TRUE
      ),
    
    q25 =
      quantile(
        value,
        0.25,
        na.rm = TRUE
      ),
    
    q75 =
      quantile(
        value,
        0.75,
        na.rm = TRUE
      )
    
  ),
  by = .(
    rel_bin,
    metric
  )
]


fwrite(
  qc_profile,
  file.path(
    OUT_DIR,
    "QC_relative_position_profile.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# Plot QC metrics
# ------------------------------------------------------------

p_qc <- ggplot(
  qc_profile,
  aes(
    x = rel_bin,
    y = median
  )
) +
  
  geom_ribbon(
    aes(
      ymin = q25,
      ymax = q75
    ),
    alpha = 0.20
  ) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  facet_wrap(
    ~ metric,
    scales = "free_y",
    ncol = 1
  ) +
  
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(
      0,
      1,
      by = 0.1
    )
  ) +
  
  labs(
    x = "Relative chromosome position",
    y = NULL,
    title = "Technical metrics across normalised chromosomes"
  ) +
  
  theme_classic()


ggsave(
  file.path(
    OUT_DIR,
    "03_QC_relative_position.png"
  ),
  p_qc,
  width = 8,
  height = 11,
  dpi = 300
)


# ============================================================
# PART K
# DEFINE CHROMOSOME TIPS AND ADJACENT ARMS
# ============================================================
#
# Example using 5%:
#
# LEFT
#
# 0-5%       = tip
# 5-25%      = adjacent arm
#
# RIGHT
#
# 75-95%     = adjacent arm
# 95-100%    = tip
#
# Centre is ignored for this particular comparison.
#
# ============================================================

assign_tip_arm <- function(
    position,
    cutoff
) {
  
  fcase(
    
    # Left tip
    position < cutoff,
    "left_tip",
    
    # Left adjacent arm
    position >= cutoff &
      position < ARM_BOUNDARY,
    "left_arm",
    
    # Right tip
    position > (1 - cutoff),
    "right_tip",
    
    # Right adjacent arm
    position > (1 - ARM_BOUNDARY) &
      position <= (1 - cutoff),
    "right_arm",
    
    default = NA_character_
  )
}


# ============================================================
# PART L
# DXY TIP VS ADJACENT ARM
# ============================================================

summarise_dxy_tips <- function(
    dat,
    cutoff,
    qc_name
) {
  
  x <- copy(dat)
  
  
  x[
    ,
    region :=
      assign_tip_arm(
        relative_position,
        cutoff
      )
  ]
  
  
  # Ignore chromosome centre
  
  x <- x[
    !is.na(region)
  ]
  
  
  # Left or right chromosome end
  
  x[
    ,
    end :=
      fifelse(
        grepl(
          "^left",
          region
        ),
        "Left",
        "Right"
      )
  ]
  
  
  # Tip or adjacent arm
  
  x[
    ,
    zone :=
      fifelse(
        grepl(
          "tip$",
          region
        ),
        "Tip",
        "Arm"
      )
  ]
  
  
  # ----------------------------------------------------------
  # Summarise each chromosome end independently
  # ----------------------------------------------------------
  
  region_summary <- x[
    ,
    .(
      
      dxy =
        ifelse(
          sum(
            dxy_comparisons,
            na.rm = TRUE
          ) > 0,
          
          sum(
            dxy_diffs,
            na.rm = TRUE
          ) /
            sum(
              dxy_comparisons,
              na.rm = TRUE
            ),
          
          NA_real_
        ),
      
      site_fraction =
        sum(
          dxy_sites,
          na.rm = TRUE
        ) /
        sum(
          bin_length,
          na.rm = TRUE
        ),
      
      depth =
        safe_wmean(
          median_depth,
          bin_length
        ),
      
      breadth =
        safe_wmean(
          median_fraction_covered,
          bin_length
        ),
      
      mappability =
        safe_wmean(
          unique_fraction,
          bin_length
        ),
      
      n_windows = .N
      
    ),
    by = .(
      chromosome,
      end,
      zone
    )
  ]
  
  
  # Separate tip and arm rows
  
  tip <- region_summary[
    zone == "Tip"
  ]
  
  
  arm <- region_summary[
    zone == "Arm"
  ]
  
  
  # Rename variables before joining
  
  setnames(
    tip,
    c(
      "dxy",
      "site_fraction",
      "depth",
      "breadth",
      "mappability",
      "n_windows"
    ),
    paste0(
      c(
        "dxy",
        "site_fraction",
        "depth",
        "breadth",
        "mappability",
        "n_windows"
      ),
      "_tip"
    )
  )
  
  
  setnames(
    arm,
    c(
      "dxy",
      "site_fraction",
      "depth",
      "breadth",
      "mappability",
      "n_windows"
    ),
    paste0(
      c(
        "dxy",
        "site_fraction",
        "depth",
        "breadth",
        "mappability",
        "n_windows"
      ),
      "_arm"
    )
  )
  
  
  # Join tip to its adjacent arm
  
  out <- merge(
    tip,
    arm,
    by = c(
      "chromosome",
      "end"
    )
  )
  
  
  # Calculate tip minus arm
  
  out[
    ,
    `:=`(
      
      delta_dxy =
        dxy_tip -
        dxy_arm,
      
      delta_site_fraction =
        site_fraction_tip -
        site_fraction_arm,
      
      delta_depth =
        depth_tip -
        depth_arm,
      
      delta_breadth =
        breadth_tip -
        breadth_arm,
      
      delta_mappability =
        mappability_tip -
        mappability_arm,
      
      cutoff =
        cutoff,
      
      QC =
        qc_name
    )
  ]
  
  
  out
}


# ------------------------------------------------------------
# Run for 2.5%, 5%, 10%
# and alternative QC datasets
# ------------------------------------------------------------

dxy_tip_results <- rbindlist(
  
  lapply(
    TIP_CUTOFFS,
    
    function(cutoff) {
      
      rbindlist(
        list(
          
          summarise_dxy_tips(
            master_dxy,
            cutoff,
            "All windows"
          ),
          
          summarise_dxy_tips(
            master_dxy[
              high_map == TRUE
            ],
            cutoff,
            "Mappability >= 90%"
          ),
          
          summarise_dxy_tips(
            master_dxy[
              strict_qc == TRUE
            ],
            cutoff,
            "Strict QC"
          )
          
        ),
        fill = TRUE
      )
    }
  ),
  
  fill = TRUE
)


fwrite(
  dxy_tip_results,
  file.path(
    OUT_DIR,
    "DXY_tip_vs_arm_by_chromosome_end.tsv"
  ),
  sep = "\t"
)


# ============================================================
# PART M
# PI TIP VS ADJACENT ARM
# ============================================================

summarise_pi_tips <- function(
    dat,
    cutoff,
    qc_name
) {
  
  x <- copy(dat)
  
  
  x[
    ,
    region :=
      assign_tip_arm(
        relative_position,
        cutoff
      )
  ]
  
  
  x <- x[
    !is.na(region)
  ]
  
  
  x[
    ,
    end :=
      fifelse(
        grepl(
          "^left",
          region
        ),
        "Left",
        "Right"
      )
  ]
  
  
  x[
    ,
    zone :=
      fifelse(
        grepl(
          "tip$",
          region
        ),
        "Tip",
        "Arm"
      )
  ]
  
  
  region_summary <- x[
    ,
    .(
      
      pi =
        ifelse(
          sum(
            pi_comparisons,
            na.rm = TRUE
          ) > 0,
          
          sum(
            pi_diffs,
            na.rm = TRUE
          ) /
            sum(
              pi_comparisons,
              na.rm = TRUE
            ),
          
          NA_real_
        ),
      
      site_fraction =
        sum(
          pi_sites,
          na.rm = TRUE
        ) /
        sum(
          bin_length,
          na.rm = TRUE
        ),
      
      depth =
        safe_wmean(
          median_depth,
          bin_length
        ),
      
      breadth =
        safe_wmean(
          median_fraction_covered,
          bin_length
        ),
      
      mappability =
        safe_wmean(
          unique_fraction,
          bin_length
        ),
      
      n_windows =
        .N
      
    ),
    by = .(
      pop,
      chromosome,
      end,
      zone
    )
  ]
  
  
  tip <- region_summary[
    zone == "Tip"
  ]
  
  
  arm <- region_summary[
    zone == "Arm"
  ]
  
  
  setnames(
    tip,
    c(
      "pi",
      "site_fraction",
      "depth",
      "breadth",
      "mappability",
      "n_windows"
    ),
    paste0(
      c(
        "pi",
        "site_fraction",
        "depth",
        "breadth",
        "mappability",
        "n_windows"
      ),
      "_tip"
    )
  )
  
  
  setnames(
    arm,
    c(
      "pi",
      "site_fraction",
      "depth",
      "breadth",
      "mappability",
      "n_windows"
    ),
    paste0(
      c(
        "pi",
        "site_fraction",
        "depth",
        "breadth",
        "mappability",
        "n_windows"
      ),
      "_arm"
    )
  )
  
  
  out <- merge(
    tip,
    arm,
    by = c(
      "pop",
      "chromosome",
      "end"
    )
  )
  
  
  out[
    ,
    `:=`(
      
      delta_pi =
        pi_tip -
        pi_arm,
      
      delta_site_fraction =
        site_fraction_tip -
        site_fraction_arm,
      
      delta_depth =
        depth_tip -
        depth_arm,
      
      delta_breadth =
        breadth_tip -
        breadth_arm,
      
      delta_mappability =
        mappability_tip -
        mappability_arm,
      
      cutoff =
        cutoff,
      
      QC =
        qc_name
    )
  ]
  
  
  out
}


pi_tip_results <- rbindlist(
  
  lapply(
    TIP_CUTOFFS,
    
    function(cutoff) {
      
      rbindlist(
        list(
          
          summarise_pi_tips(
            master_pi,
            cutoff,
            "All windows"
          ),
          
          summarise_pi_tips(
            master_pi[
              high_map == TRUE
            ],
            cutoff,
            "Mappability >= 90%"
          ),
          
          summarise_pi_tips(
            master_pi[
              strict_qc == TRUE
            ],
            cutoff,
            "Strict QC"
          )
          
        ),
        fill = TRUE
      )
    }
  ),
  
  fill = TRUE
)


fwrite(
  pi_tip_results,
  file.path(
    OUT_DIR,
    "PI_tip_vs_arm_by_chromosome_end.tsv"
  ),
  sep = "\t"
)


# ============================================================
# PART N
# STATISTICAL TEST:
# IS DXY LOWER AT TIPS THAN ADJACENT ARMS?
# ============================================================
#
# Negative delta:
#
#       tip dXY < adjacent-arm dXY
#
# We average left and right tip effects for each chromosome,
# then use chromosome as the inferential unit.
#
# ============================================================

test_dxy_tip_decline <- function(dat) {
  
  
  chromosome_values <- dat[
    !is.na(delta_dxy),
    .(
      
      delta =
        mean(
          delta_dxy,
          na.rm = TRUE
        )
      
    ),
    by = chromosome
  ]
  
  
  if (nrow(chromosome_values) < 4) {
    
    return(
      data.table(
        
        n_chromosomes =
          nrow(chromosome_values),
        
        median_delta =
          NA_real_,
        
        chromosomes_negative =
          NA_integer_,
        
        proportion_negative =
          NA_real_,
        
        wilcoxon_p =
          NA_real_
      )
    )
  }
  
  
  wt <- suppressWarnings(
    wilcox.test(
      chromosome_values$delta,
      mu = 0,
      exact = FALSE
    )
  )
  
  
  data.table(
    
    n_chromosomes =
      nrow(chromosome_values),
    
    median_delta =
      median(
        chromosome_values$delta,
        na.rm = TRUE
      ),
    
    chromosomes_negative =
      sum(
        chromosome_values$delta < 0,
        na.rm = TRUE
      ),
    
    proportion_negative =
      mean(
        chromosome_values$delta < 0,
        na.rm = TRUE
      ),
    
    wilcoxon_p =
      wt$p.value
  )
}


dxy_tests <- dxy_tip_results[
  ,
  test_dxy_tip_decline(.SD),
  by = .(
    cutoff,
    QC
  )
]


fwrite(
  dxy_tests,
  file.path(
    OUT_DIR,
    "DXY_tip_decline_tests.tsv"
  ),
  sep = "\t"
)


# ============================================================
# PART O
# STATISTICAL TEST:
# IS PI LOWER AT TIPS THAN ADJACENT ARMS?
# ============================================================

test_pi_tip_decline <- function(dat) {
  
  
  chromosome_values <- dat[
    !is.na(delta_pi),
    .(
      
      delta =
        mean(
          delta_pi,
          na.rm = TRUE
        )
      
    ),
    by = chromosome
  ]
  
  
  if (nrow(chromosome_values) < 4) {
    
    return(
      data.table(
        
        n_chromosomes =
          nrow(chromosome_values),
        
        median_delta =
          NA_real_,
        
        chromosomes_negative =
          NA_integer_,
        
        proportion_negative =
          NA_real_,
        
        wilcoxon_p =
          NA_real_
      )
    )
  }
  
  
  wt <- suppressWarnings(
    wilcox.test(
      chromosome_values$delta,
      mu = 0,
      exact = FALSE
    )
  )
  
  
  data.table(
    
    n_chromosomes =
      nrow(chromosome_values),
    
    median_delta =
      median(
        chromosome_values$delta,
        na.rm = TRUE
      ),
    
    chromosomes_negative =
      sum(
        chromosome_values$delta < 0,
        na.rm = TRUE
      ),
    
    proportion_negative =
      mean(
        chromosome_values$delta < 0,
        na.rm = TRUE
      ),
    
    wilcoxon_p =
      wt$p.value
  )
}


pi_tests <- pi_tip_results[
  ,
  test_pi_tip_decline(.SD),
  by = .(
    pop,
    cutoff,
    QC
  )
]


fwrite(
  pi_tests,
  file.path(
    OUT_DIR,
    "PI_tip_decline_tests.tsv"
  ),
  sep = "\t"
)


# ============================================================
# PART P
# DOES THE TECHNICAL DROP TRACK THE DXY DROP?
# ============================================================

correlation_variables <- c(
  "delta_depth",
  "delta_breadth",
  "delta_mappability",
  "delta_site_fraction"
)


dxy_correlations <- rbindlist(
  
  lapply(
    TIP_CUTOFFS,
    
    function(cutoff_value) {
      
      
      temp <- dxy_tip_results[
        cutoff == cutoff_value &
          QC == "All windows"
      ]
      
      
      rbindlist(
        
        lapply(
          correlation_variables,
          
          function(v) {
            
            
            ans <- safe_spearman(
              temp$delta_dxy,
              temp[[v]]
            )
            
            
            ans[
              ,
              `:=`(
                
                cutoff =
                  cutoff_value,
                
                predictor =
                  v
              )
            ]
            
            
            ans
          }
        )
      )
    }
  )
)


fwrite(
  dxy_correlations,
  file.path(
    OUT_DIR,
    "DXY_tip_delta_QC_correlations.tsv"
  ),
  sep = "\t"
)


# ============================================================
# PART Q
# PI DROP VS TECHNICAL DROP
# ============================================================

pi_correlations <- rbindlist(
  
  lapply(
    unique(
      pi_tip_results$pop
    ),
    
    function(pop_value) {
      
      
      rbindlist(
        
        lapply(
          TIP_CUTOFFS,
          
          function(cutoff_value) {
            
            
            temp <- pi_tip_results[
              pop == pop_value &
                cutoff == cutoff_value &
                QC == "All windows"
            ]
            
            
            rbindlist(
              
              lapply(
                correlation_variables,
                
                function(v) {
                  
                  
                  ans <- safe_spearman(
                    temp$delta_pi,
                    temp[[v]]
                  )
                  
                  
                  ans[
                    ,
                    `:=`(
                      
                      pop =
                        pop_value,
                      
                      cutoff =
                        cutoff_value,
                      
                      predictor =
                        v
                    )
                  ]
                  
                  
                  ans
                }
              )
            )
          }
        )
      )
    }
  )
)


fwrite(
  pi_correlations,
  file.path(
    OUT_DIR,
    "PI_tip_delta_QC_correlations.tsv"
  ),
  sep = "\t"
)


# ============================================================
# PART R
# HOW MANY TIP WINDOWS SURVIVE QC?
# ============================================================

master_dxy[
  ,
  region_5pct :=
    assign_tip_arm(
      relative_position,
      0.05
    )
]


retention_5pct <- master_dxy[
  !is.na(region_5pct),
  .(
    
    total_windows =
      .N,
    
    map90_windows =
      sum(
        high_map,
        na.rm = TRUE
      ),
    
    strict_windows =
      sum(
        strict_qc,
        na.rm = TRUE
      )
    
  ),
  by = region_5pct
]


retention_5pct[
  ,
  `:=`(
    
    proportion_map90 =
      map90_windows /
      total_windows,
    
    proportion_strict =
      strict_windows /
      total_windows
  )
]


fwrite(
  retention_5pct,
  file.path(
    OUT_DIR,
    "DXY_5pct_tip_QC_retention.tsv"
  ),
  sep = "\t"
)


# ============================================================
# PART S
# DXY TIP-ARM DELTA PLOT
# ============================================================

dxy_tip_plot_data <- copy(
  dxy_tip_results
)


dxy_tip_plot_data[
  ,
  cutoff_label :=
    paste0(
      cutoff * 100,
      "% tips"
    )
]


p_dxy_delta <- ggplot(
  dxy_tip_plot_data,
  aes(
    x = cutoff_label,
    y = delta_dxy
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  
  geom_boxplot(
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.12,
    alpha = 0.6
  ) +
  
  facet_wrap(
    ~ QC,
    ncol = 1
  ) +
  
  labs(
    x = "Tip definition",
    y = expression(
      Delta*d[XY]~"(tip - adjacent arm)"
    ),
    title = "Human vs pig dXY at chromosome tips"
  ) +
  
  theme_classic()


ggsave(
  file.path(
    OUT_DIR,
    "04_DXY_tip_arm_deltas.png"
  ),
  p_dxy_delta,
  width = 8,
  height = 9,
  dpi = 300
)


# ============================================================
# PART T
# PI TIP-ARM DELTA PLOT
# ============================================================

pi_tip_plot_data <- copy(
  pi_tip_results
)


pi_tip_plot_data[
  ,
  cutoff_label :=
    paste0(
      cutoff * 100,
      "% tips"
    )
]


p_pi_delta <- ggplot(
  pi_tip_plot_data,
  aes(
    x = cutoff_label,
    y = delta_pi
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  
  geom_boxplot(
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.12,
    alpha = 0.6
  ) +
  
  facet_grid(
    pop ~ QC
  ) +
  
  labs(
    x = "Tip definition",
    y = expression(
      Delta*pi~"(tip - adjacent arm)"
    ),
    title = "Nucleotide diversity at chromosome tips"
  ) +
  
  theme_classic()


ggsave(
  file.path(
    OUT_DIR,
    "05_PI_tip_arm_deltas.png"
  ),
  p_pi_delta,
  width = 12,
  height = 10,
  dpi = 300
)


# ============================================================
# PART U
# PRINT THE KEY RESULTS
# ============================================================

cat(
  "\n\n========================================\n",
  "DXY TIP DECLINE TESTS\n",
  "========================================\n"
)


print(
  dxy_tests
)


cat(
  "\n\n========================================\n",
  "PI TIP DECLINE TESTS\n",
  "========================================\n"
)


print(
  pi_tests
)


cat(
  "\n\n========================================\n",
  "5% TIP QC RETENTION\n",
  "========================================\n"
)


print(
  retention_5pct
)


cat(
  "\n\n========================================\n",
  "DXY TIP DROP VS QC METRICS\n",
  "========================================\n"
)


print(
  dxy_correlations
)


cat(
  "\n\n========================================\n",
  "PI TIP DROP VS QC METRICS\n",
  "========================================\n"
)


print(
  pi_correlations
)


cat(
  "\n\n========================================\n",
  "ANALYSIS COMPLETE\n",
  "========================================\n",
  "\nResults saved to:\n",
  OUT_DIR,
  "\n"
)

# ============================================================
# TERMINAL GRADIENT ANALYSIS
# ============================================================
#
# Purpose:
#
# Examine the actual transition from chromosome tip -> arm.
#
# Distance from the nearest chromosome end is divided into:
#
#   0-2.5%
#   2.5-5%
#   5-7.5%
#   7.5-10%
#   10-15%
#   15-20%
#   20-25%
#
# Diversity/divergence is examined under STRICT QC.
#
# Technical QC metrics are examined using ALL WINDOWS so that
# we can see where depth/mappability/coverage actually recover.
#
# For formal testing, each terminal band from 0-15% is compared
# with an internal arm reference spanning 15-25%.
#
# Left and right chromosome ends are first compared separately,
# then averaged within chromosome so the inferential unit
# remains the chromosome (maximum n = 19).
#
# ============================================================


# ------------------------------------------------------------
# 1. TERMINAL BAND DEFINITIONS
# ------------------------------------------------------------

terminal_breaks <- c(
  0,
  0.025,
  0.05,
  0.075,
  0.10,
  0.15,
  0.20,
  0.2500001
)


terminal_labels <- c(
  "0-2.5%",
  "2.5-5%",
  "5-7.5%",
  "7.5-10%",
  "10-15%",
  "15-20%",
  "20-25%"
)


# Ordered factor so plots appear in the correct order

terminal_levels <- terminal_labels


# ============================================================
# 2. ADD DISTANCE-FROM-END VARIABLES
# ============================================================

# ------------------------------------------------------------
# dXY
# ------------------------------------------------------------

master_dxy[
  ,
  terminal_distance :=
    distance_from_end /
    chrom_length
]


master_dxy[
  ,
  chromosome_end :=
    fifelse(
      midpoint <= chrom_length / 2,
      "Left",
      "Right"
    )
]


master_dxy[
  ,
  terminal_band :=
    cut(
      terminal_distance,
      breaks = terminal_breaks,
      labels = terminal_labels,
      include.lowest = TRUE,
      right = FALSE
    )
]


master_dxy[
  ,
  terminal_band :=
    factor(
      terminal_band,
      levels = terminal_levels
    )
]


# ------------------------------------------------------------
# pi
# ------------------------------------------------------------

master_pi[
  ,
  terminal_distance :=
    distance_from_end /
    chrom_length
]


master_pi[
  ,
  chromosome_end :=
    fifelse(
      midpoint <= chrom_length / 2,
      "Left",
      "Right"
    )
]


master_pi[
  ,
  terminal_band :=
    cut(
      terminal_distance,
      breaks = terminal_breaks,
      labels = terminal_labels,
      include.lowest = TRUE,
      right = FALSE
    )
]


master_pi[
  ,
  terminal_band :=
    factor(
      terminal_band,
      levels = terminal_levels
    )
]


# ============================================================
# 3. STRICT-QC DXY BY TERMINAL BAND
# ============================================================

# First calculate one dXY estimate for each:
#
# chromosome x chromosome end x terminal band

dxy_terminal_end <- master_dxy[
  strict_qc == TRUE &
    !is.na(terminal_band) &
    dxy_comparisons > 0,
  .(
    
    dxy =
      sum(
        dxy_diffs,
        na.rm = TRUE
      ) /
      sum(
        dxy_comparisons,
        na.rm = TRUE
      ),
    
    n_windows = .N
    
  ),
  by = .(
    chromosome,
    chromosome_end,
    terminal_band
  )
]


# ------------------------------------------------------------
# Average left and right ends within each chromosome
# ------------------------------------------------------------

dxy_terminal_chr <- dxy_terminal_end[
  ,
  .(
    
    dxy =
      mean(
        dxy,
        na.rm = TRUE
      ),
    
    n_ends =
      uniqueN(
        chromosome_end
      )
    
  ),
  by = .(
    chromosome,
    terminal_band
  )
]


# ------------------------------------------------------------
# Median profile across chromosomes
# ------------------------------------------------------------

dxy_terminal_profile <- dxy_terminal_chr[
  ,
  .(
    
    median =
      median(
        dxy,
        na.rm = TRUE
      ),
    
    q25 =
      quantile(
        dxy,
        0.25,
        na.rm = TRUE
      ),
    
    q75 =
      quantile(
        dxy,
        0.75,
        na.rm = TRUE
      ),
    
    n_chr =
      uniqueN(
        chromosome
      )
    
  ),
  by = terminal_band
]


print(
  dxy_terminal_profile
)


# ============================================================
# 4. STRICT-QC PI BY TERMINAL BAND
# ============================================================

# One estimate for each:
#
# population x chromosome x chromosome end x terminal band

pi_terminal_end <- master_pi[
  strict_qc == TRUE &
    !is.na(terminal_band) &
    pi_comparisons > 0,
  .(
    
    pi =
      sum(
        pi_diffs,
        na.rm = TRUE
      ) /
      sum(
        pi_comparisons,
        na.rm = TRUE
      ),
    
    n_windows = .N
    
  ),
  by = .(
    pop,
    chromosome,
    chromosome_end,
    terminal_band
  )
]


# Average left and right ends within chromosome

pi_terminal_chr <- pi_terminal_end[
  ,
  .(
    
    pi =
      mean(
        pi,
        na.rm = TRUE
      ),
    
    n_ends =
      uniqueN(
        chromosome_end
      )
    
  ),
  by = .(
    pop,
    chromosome,
    terminal_band
  )
]


# Median across chromosomes

pi_terminal_profile <- pi_terminal_chr[
  ,
  .(
    
    median =
      median(
        pi,
        na.rm = TRUE
      ),
    
    q25 =
      quantile(
        pi,
        0.25,
        na.rm = TRUE
      ),
    
    q75 =
      quantile(
        pi,
        0.75,
        na.rm = TRUE
      ),
    
    n_chr =
      uniqueN(
        chromosome
      )
    
  ),
  by = .(
    pop,
    terminal_band
  )
]


print(
  pi_terminal_profile
)


# ============================================================
# 5. TECHNICAL QC ACROSS THE SAME TERMINAL BANDS
# ============================================================
#
# IMPORTANT:
#
# Use ALL WINDOWS here rather than strict-QC windows.
#
# Otherwise we would artificially flatten the technical
# profiles by removing exactly the poor-quality sequence
# we are trying to examine.
#
# ============================================================

qc_terminal <- copy(
  master_dxy
)


qc_terminal <- qc_terminal[
  !is.na(terminal_band)
]


# First one technical value per chromosome end / band

qc_terminal_end <- qc_terminal[
  ,
  .(
    
    median_depth =
      median(
        median_depth,
        na.rm = TRUE
      ),
    
    coverage_breadth =
      median(
        median_fraction_covered,
        na.rm = TRUE
      ),
    
    mappability =
      median(
        unique_fraction,
        na.rm = TRUE
      ),
    
    dxy_site_fraction =
      median(
        dxy_site_fraction,
        na.rm = TRUE
      ),
    
    n_windows = .N
    
  ),
  by = .(
    chromosome,
    chromosome_end,
    terminal_band
  )
]


# Average both chromosome ends within each chromosome

qc_terminal_chr <- qc_terminal_end[
  ,
  .(
    
    median_depth =
      mean(
        median_depth,
        na.rm = TRUE
      ),
    
    coverage_breadth =
      mean(
        coverage_breadth,
        na.rm = TRUE
      ),
    
    mappability =
      mean(
        mappability,
        na.rm = TRUE
      ),
    
    dxy_site_fraction =
      mean(
        dxy_site_fraction,
        na.rm = TRUE
      )
    
  ),
  by = .(
    chromosome,
    terminal_band
  )
]


# Median across chromosomes

qc_terminal_profile <- qc_terminal_chr[
  ,
  .(
    
    median_depth =
      median(
        median_depth,
        na.rm = TRUE
      ),
    
    coverage_breadth =
      median(
        coverage_breadth,
        na.rm = TRUE
      ),
    
    mappability =
      median(
        mappability,
        na.rm = TRUE
      ),
    
    dxy_site_fraction =
      median(
        dxy_site_fraction,
        na.rm = TRUE
      ),
    
    n_chr =
      uniqueN(
        chromosome
      )
    
  ),
  by = terminal_band
]


print(
  qc_terminal_profile
)


# ============================================================
# 6. FORMAL DXY COMPARISON AGAINST ARM REFERENCE
# ============================================================
#
# Reference:
#
#        15-25% from chromosome end
#
# This combines:
#
#        15-20%
#        20-25%
#
# For every chromosome END independently.
#
# Then:
#
# delta = terminal band - arm reference
#
# Negative delta means the terminal band has LOWER dXY than
# the established arm region.
#
# ============================================================


# ------------------------------------------------------------
# Build 15-25% arm reference for each chromosome end
# ------------------------------------------------------------

dxy_arm_reference <- master_dxy[
  strict_qc == TRUE &
    terminal_distance >= 0.15 &
    terminal_distance < 0.25 &
    dxy_comparisons > 0,
  .(
    
    arm_dxy =
      sum(
        dxy_diffs,
        na.rm = TRUE
      ) /
      sum(
        dxy_comparisons,
        na.rm = TRUE
      )
    
  ),
  by = .(
    chromosome,
    chromosome_end
  )
]


# Only test bands BEFORE the arm reference

dxy_test_bands <- c(
  "0-2.5%",
  "2.5-5%",
  "5-7.5%",
  "7.5-10%",
  "10-15%"
)


dxy_band_tests_end <- dxy_terminal_end[
  terminal_band %in% dxy_test_bands
]


dxy_band_tests_end <- merge(
  dxy_band_tests_end,
  dxy_arm_reference,
  by = c(
    "chromosome",
    "chromosome_end"
  )
)


dxy_band_tests_end[
  ,
  delta_from_arm :=
    dxy -
    arm_dxy
]


# ------------------------------------------------------------
# Average left/right effects within chromosome
# ------------------------------------------------------------

dxy_band_tests_chr <- dxy_band_tests_end[
  ,
  .(
    
    delta =
      mean(
        delta_from_arm,
        na.rm = TRUE
      ),
    
    n_ends =
      uniqueN(
        chromosome_end
      )
    
  ),
  by = .(
    chromosome,
    terminal_band
  )
]


# ------------------------------------------------------------
# Wilcoxon test for each terminal band
# ------------------------------------------------------------

dxy_gradient_tests <- dxy_band_tests_chr[
  ,
  {
    
    x <- delta[
      !is.na(delta)
    ]
    
    if (length(x) >= 4) {
      
      wt <- suppressWarnings(
        wilcox.test(
          x,
          mu = 0,
          exact = FALSE
        )
      )
      
      list(
        
        n_chromosomes =
          length(x),
        
        median_delta =
          median(
            x,
            na.rm = TRUE
          ),
        
        chromosomes_lower =
          sum(
            x < 0,
            na.rm = TRUE
          ),
        
        proportion_lower =
          mean(
            x < 0,
            na.rm = TRUE
          ),
        
        p_value =
          wt$p.value
      )
      
    } else {
      
      list(
        
        n_chromosomes =
          length(x),
        
        median_delta =
          NA_real_,
        
        chromosomes_lower =
          NA_integer_,
        
        proportion_lower =
          NA_real_,
        
        p_value =
          NA_real_
      )
    }
    
  },
  by = terminal_band
]


# Correct the five tests using Benjamini-Hochberg

dxy_gradient_tests[
  ,
  p_BH :=
    p.adjust(
      p_value,
      method = "BH"
    )
]


print(
  dxy_gradient_tests
)


# ============================================================
# 7. FORMAL PI COMPARISON AGAINST ARM REFERENCE
# ============================================================

pi_arm_reference <- master_pi[
  strict_qc == TRUE &
    terminal_distance >= 0.15 &
    terminal_distance < 0.25 &
    pi_comparisons > 0,
  .(
    
    arm_pi =
      sum(
        pi_diffs,
        na.rm = TRUE
      ) /
      sum(
        pi_comparisons,
        na.rm = TRUE
      )
    
  ),
  by = .(
    pop,
    chromosome,
    chromosome_end
  )
]


pi_band_tests_end <- pi_terminal_end[
  terminal_band %in% dxy_test_bands
]


pi_band_tests_end <- merge(
  pi_band_tests_end,
  pi_arm_reference,
  by = c(
    "pop",
    "chromosome",
    "chromosome_end"
  )
)


pi_band_tests_end[
  ,
  delta_from_arm :=
    pi -
    arm_pi
]


# Average left/right within each chromosome

pi_band_tests_chr <- pi_band_tests_end[
  ,
  .(
    
    delta =
      mean(
        delta_from_arm,
        na.rm = TRUE
      ),
    
    n_ends =
      uniqueN(
        chromosome_end
      )
    
  ),
  by = .(
    pop,
    chromosome,
    terminal_band
  )
]


# Test each population separately

pi_gradient_tests <- pi_band_tests_chr[
  ,
  {
    
    x <- delta[
      !is.na(delta)
    ]
    
    if (length(x) >= 4) {
      
      wt <- suppressWarnings(
        wilcox.test(
          x,
          mu = 0,
          exact = FALSE
        )
      )
      
      list(
        
        n_chromosomes =
          length(x),
        
        median_delta =
          median(
            x,
            na.rm = TRUE
          ),
        
        chromosomes_lower =
          sum(
            x < 0,
            na.rm = TRUE
          ),
        
        proportion_lower =
          mean(
            x < 0,
            na.rm = TRUE
          ),
        
        p_value =
          wt$p.value
      )
      
    } else {
      
      list(
        
        n_chromosomes =
          length(x),
        
        median_delta =
          NA_real_,
        
        chromosomes_lower =
          NA_integer_,
        
        proportion_lower =
          NA_real_,
        
        p_value =
          NA_real_
      )
    }
    
  },
  by = .(
    pop,
    terminal_band
  )
]


# BH correction separately within population

pi_gradient_tests[
  ,
  p_BH :=
    p.adjust(
      p_value,
      method = "BH"
    ),
  by = pop
]


print(
  pi_gradient_tests
)


# ============================================================
# 8. PLOT STRICT-QC DXY TERMINAL GRADIENT
# ============================================================

p_dxy_terminal <- ggplot(
  dxy_terminal_profile,
  aes(
    x = terminal_band,
    y = median,
    group = 1
  )
) +
  
  geom_ribbon(
    aes(
      ymin = q25,
      ymax = q75,
      group = 1
    ),
    alpha = 0.20
  ) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 2
  ) +
  
  labs(
    x = "Distance from nearest chromosome end",
    y = expression(d[XY]),
    title = "Human vs pig dXY from chromosome termini into the arms",
    subtitle = "Strict-QC windows only"
  ) +
  
  theme_classic()


ggsave(
  file.path(
    OUT_DIR,
    "06_DXY_terminal_gradient_strictQC.png"
  ),
  p_dxy_terminal,
  width = 8,
  height = 5,
  dpi = 300
)


# ============================================================
# 9. PLOT STRICT-QC PI TERMINAL GRADIENT
# ============================================================

p_pi_terminal <- ggplot(
  pi_terminal_profile,
  aes(
    x = terminal_band,
    y = median,
    group = pop,
    linetype = pop
  )
) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 2
  ) +
  
  labs(
    x = "Distance from nearest chromosome end",
    y = expression(pi),
    linetype = "Population",
    title = "Nucleotide diversity from chromosome termini into the arms",
    subtitle = "Strict-QC windows only"
  ) +
  
  theme_classic()


ggsave(
  file.path(
    OUT_DIR,
    "07_PI_terminal_gradient_strictQC.png"
  ),
  p_pi_terminal,
  width = 8,
  height = 5,
  dpi = 300
)


# ============================================================
# 10. TECHNICAL-QC TERMINAL GRADIENT PLOT
# ============================================================

qc_terminal_long <- melt(
  qc_terminal_profile,
  id.vars = c(
    "terminal_band",
    "n_chr"
  ),
  measure.vars = c(
    "median_depth",
    "coverage_breadth",
    "mappability",
    "dxy_site_fraction"
  ),
  variable.name = "metric",
  value.name = "value"
)


p_qc_terminal <- ggplot(
  qc_terminal_long,
  aes(
    x = terminal_band,
    y = value,
    group = 1
  )
) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 2
  ) +
  
  facet_wrap(
    ~ metric,
    scales = "free_y",
    ncol = 1
  ) +
  
  labs(
    x = "Distance from nearest chromosome end",
    y = NULL,
    title = "Technical QC from chromosome termini into the arms",
    subtitle = "All windows"
  ) +
  
  theme_classic()


ggsave(
  file.path(
    OUT_DIR,
    "08_QC_terminal_gradient.png"
  ),
  p_qc_terminal,
  width = 8,
  height = 10,
  dpi = 300
)


# ============================================================
# 11. SAVE TABLES
# ============================================================

fwrite(
  dxy_terminal_profile,
  file.path(
    OUT_DIR,
    "DXY_terminal_gradient_profile.tsv"
  ),
  sep = "\t"
)


fwrite(
  pi_terminal_profile,
  file.path(
    OUT_DIR,
    "PI_terminal_gradient_profile.tsv"
  ),
  sep = "\t"
)


fwrite(
  qc_terminal_profile,
  file.path(
    OUT_DIR,
    "QC_terminal_gradient_profile.tsv"
  ),
  sep = "\t"
)


fwrite(
  dxy_gradient_tests,
  file.path(
    OUT_DIR,
    "DXY_terminal_gradient_tests.tsv"
  ),
  sep = "\t"
)


fwrite(
  pi_gradient_tests,
  file.path(
    OUT_DIR,
    "PI_terminal_gradient_tests.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 12. PRINT THE IMPORTANT OUTPUT
# ============================================================

cat(
  "\n\n========================================\n",
  "STRICT-QC DXY TERMINAL PROFILE\n",
  "========================================\n"
)

print(
  dxy_terminal_profile
)


cat(
  "\n\n========================================\n",
  "STRICT-QC PI TERMINAL PROFILE\n",
  "========================================\n"
)

print(
  pi_terminal_profile
)


cat(
  "\n\n========================================\n",
  "TECHNICAL QC TERMINAL PROFILE\n",
  "========================================\n"
)

print(
  qc_terminal_profile
)


cat(
  "\n\n========================================\n",
  "DXY VS 15-25% ARM REFERENCE\n",
  "========================================\n"
)

print(
  dxy_gradient_tests
)


cat(
  "\n\n========================================\n",
  "PI VS 15-25% ARM REFERENCE\n",
  "========================================\n"
)

print(
  pi_gradient_tests
)

install.packages("patchwork")   # only if not already installed

library(data.table)
library(ggplot2)
library(patchwork)

# ============================================================
# FINAL COMBINED SUPPLEMENTARY FIGURE
# ============================================================

library(ggplot2)
library(data.table)
library(patchwork)


# ------------------------------------------------------------
# A. FULL NORMALISED dXY PROFILE - STRICT QC
# ------------------------------------------------------------

plot_dxy_full <- dxy_profile[
  QC == "Strict QC"
]

pA <- ggplot(
  plot_dxy_full,
  aes(
    x = rel_bin,
    y = median
  )
) +
  geom_ribbon(
    aes(
      ymin = q25,
      ymax = q75
    ),
    alpha = 0.20
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    x = "Relative chromosome position",
    y = expression(d[XY])
  ) +
  theme_classic(base_size = 11)


# ------------------------------------------------------------
# B. FULL NORMALISED PI PROFILE - STRICT QC
# ------------------------------------------------------------

plot_pi_full <- pi_profile[
  QC == "Strict QC"
]

pB <- ggplot(
  plot_pi_full,
  aes(
    x = rel_bin,
    y = median,
    group = pop,
    linetype = pop
  )
) +
  geom_line(
    linewidth = 0.9
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    x = "Relative chromosome position",
    y = expression(pi),
    linetype = "Region"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  )


# ------------------------------------------------------------
# C. dXY TERMINAL GRADIENT - STRICT QC
# ------------------------------------------------------------

pC <- ggplot(
  dxy_terminal_profile,
  aes(
    x = terminal_band,
    y = median,
    group = 1
  )
) +
  geom_ribbon(
    aes(
      ymin = q25,
      ymax = q75,
      group = 1
    ),
    alpha = 0.20
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 2
  ) +
  labs(
    x = "Distance from nearest chromosome end",
    y = expression(d[XY])
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


# ------------------------------------------------------------
# D. PI TERMINAL GRADIENT - STRICT QC
#
# Duplicate Region legend removed.
# ------------------------------------------------------------

pD <- ggplot(
  pi_terminal_profile,
  aes(
    x = terminal_band,
    y = median,
    group = pop,
    linetype = pop
  )
) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 2
  ) +
  labs(
    x = "Distance from nearest chromosome end",
    y = expression(pi)
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


# ------------------------------------------------------------
# E. TECHNICAL QC TERMINAL GRADIENT
#
# Uses ALL windows, not just strict-QC windows.
# ------------------------------------------------------------

qc_terminal_long <- melt(
  qc_terminal_profile,
  id.vars = c(
    "terminal_band",
    "n_chr"
  ),
  measure.vars = c(
    "median_depth",
    "coverage_breadth",
    "mappability",
    "dxy_site_fraction"
  ),
  variable.name = "metric",
  value.name = "value"
)


# Rename facet labels
qc_terminal_long[
  ,
  metric := factor(
    metric,
    levels = c(
      "median_depth",
      "coverage_breadth",
      "mappability",
      "dxy_site_fraction"
    ),
    labels = c(
      "Median depth",
      "Coverage breadth",
      "Unique mappability",
      "dXY site fraction"
    )
  )
]


pE <- ggplot(
  qc_terminal_long,
  aes(
    x = terminal_band,
    y = value,
    group = 1
  )
) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 1.8
  ) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    nrow = 1
  ) +
  labs(
    x = "Distance from nearest chromosome end",
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold"
    )
  )


# ============================================================
# COMBINE PANELS
# ============================================================

supp_figure <- (
  pA | pB
) /
  (
    pC | pD
  ) /
  pE +
  plot_annotation(
    tag_levels = "A"
  ) +
  plot_layout(
    heights = c(
      1,
      1,
      1.15
    )
  )


supp_figure


# ============================================================
# SAVE
# ============================================================

ggsave(
  file.path(
    OUT_DIR,
    "Figure_S_terminal_chromosome_QC_combined.png"
  ),
  supp_figure,
  width = 12,
  height = 12,
  dpi = 600
)


ggsave(
  file.path(
    OUT_DIR,
    "Figure_S_terminal_chromosome_QC_combined.pdf"
  ),
  supp_figure,
  width = 12,
  height = 12
)

# ============================================================
# SIMPLIFIED SUPPLEMENTARY TABLE
# ============================================================

bands_keep <- c(
  "0-2.5%",
  "2.5-5%",
  "5-7.5%"
)


# ------------------------------------------------------------
# dXY
# ------------------------------------------------------------

dxy_supp <- copy(
  dxy_gradient_tests[
    terminal_band %in% bands_keep
  ]
)

dxy_supp[
  ,
  `:=`(
    Metric = "dXY",
    Region = "Human vs Pig"
  )
]


# ------------------------------------------------------------
# PI
# ------------------------------------------------------------

pi_supp <- copy(
  pi_gradient_tests[
    terminal_band %in% bands_keep
  ]
)

pi_supp[
  ,
  `:=`(
    Metric = "Nucleotide diversity (pi)",
    Region = pop
  )
]


# ------------------------------------------------------------
# Put them into same format
# ------------------------------------------------------------

dxy_supp_final <- dxy_supp[
  ,
  .(
    Metric,
    Region,
    `Distance from chromosome end` =
      as.character(terminal_band),
    
    `n chromosomes` =
      n_chromosomes,
    
    `Median difference from arm` =
      median_delta,
    
    `Chromosomes with lower value` =
      chromosomes_lower,
    
    `Proportion lower` =
      proportion_lower,
    
    `BH-adjusted p` =
      p_BH
  )
]


pi_supp_final <- pi_supp[
  ,
  .(
    Metric,
    Region,
    `Distance from chromosome end` =
      as.character(terminal_band),
    
    `n chromosomes` =
      n_chromosomes,
    
    `Median difference from arm` =
      median_delta,
    
    `Chromosomes with lower value` =
      chromosomes_lower,
    
    `Proportion lower` =
      proportion_lower,
    
    `BH-adjusted p` =
      p_BH
  )
]


supp_table <- rbindlist(
  list(
    dxy_supp_final,
    pi_supp_final
  ),
  use.names = TRUE
)


# ------------------------------------------------------------
# Make ordering sensible
# ------------------------------------------------------------

supp_table[
  ,
  `Distance from chromosome end` :=
    factor(
      `Distance from chromosome end`,
      levels = bands_keep
    )
]


setorder(
  supp_table,
  Metric,
  Region,
  `Distance from chromosome end`
)


# ------------------------------------------------------------
# Round for publication
# ------------------------------------------------------------

supp_table[
  ,
  `Median difference from arm` :=
    round(
      `Median difference from arm`,
      5
    )
]


supp_table[
  ,
  `Proportion lower` :=
    round(
      `Proportion lower`,
      3
    )
]


supp_table[
  ,
  `BH-adjusted p` :=
    signif(
      `BH-adjusted p`,
      3
    )
]


print(
  supp_table
)


# Save
fwrite(
  supp_table,
  file.path(
    OUT_DIR,
    "Table_S_terminal_diversity_tests.tsv"
  ),
  sep = "\t"
)

write.csv(
  supp_table,
  file.path(
    OUT_DIR,
    "Table_S_terminal_diversity_tests.csv"
  ),
  row.names = FALSE
)

# ============================================================
# CREATE SUPPLEMENTARY TERMINAL QC + DIVERSITY TABLE
# ============================================================

library(data.table)

# ------------------------------------------------------------
# 1. PI: convert Africa / Asia / Europe from long to wide
# ------------------------------------------------------------

pi_table <- dcast(
  pi_terminal_profile,
  terminal_band ~ pop,
  value.var = "median"
)

# Rename population columns so it is obvious these are pi values
setnames(
  pi_table,
  old = c("Africa", "Asia", "Europe"),
  new = c(
    "Africa_pi",
    "Asia_pi",
    "Europe_pi"
  )
)


# ------------------------------------------------------------
# 2. dXY
# ------------------------------------------------------------

dxy_table <- dxy_terminal_profile[
  ,
  .(
    terminal_band,
    n_chromosomes = n_chr,
    Human_Pig_dXY = median
  )
]


# ------------------------------------------------------------
# 3. TECHNICAL QC
# ------------------------------------------------------------
#
# These are the ALL-WINDOW technical summaries,
# not the strict-QC subset.
# ------------------------------------------------------------

qc_table <- qc_terminal_profile[
  ,
  .(
    terminal_band,
    median_depth,
    coverage_breadth,
    unique_mappability = mappability,
    dXY_site_fraction = dxy_site_fraction
  )
]


# ------------------------------------------------------------
# 4. MERGE EVERYTHING
# ------------------------------------------------------------

supp_terminal_table <- merge(
  dxy_table,
  pi_table,
  by = "terminal_band",
  all = TRUE
)

supp_terminal_table <- merge(
  supp_terminal_table,
  qc_table,
  by = "terminal_band",
  all = TRUE
)


# ------------------------------------------------------------
# 5. PUT BANDS IN CORRECT BIOLOGICAL ORDER
# ------------------------------------------------------------

band_order <- c(
  "0-2.5%",
  "2.5-5%",
  "5-7.5%",
  "7.5-10%",
  "10-15%",
  "15-20%",
  "20-25%"
)

supp_terminal_table[
  ,
  terminal_band :=
    factor(
      terminal_band,
      levels = band_order
    )
]

setorder(
  supp_terminal_table,
  terminal_band
)


# ------------------------------------------------------------
# 6. CONVERT PROPORTIONS TO PERCENTAGES
# ------------------------------------------------------------
#
# Keep depth as x coverage.
# Convert breadth, mappability and site fraction to %.
# ------------------------------------------------------------

supp_terminal_table[
  ,
  coverage_breadth_percent :=
    coverage_breadth * 100
]

supp_terminal_table[
  ,
  unique_mappability_percent :=
    unique_mappability * 100
]

supp_terminal_table[
  ,
  dXY_site_fraction_percent :=
    dXY_site_fraction * 100
]


# Remove original proportion columns now percentages exist

supp_terminal_table[
  ,
  c(
    "coverage_breadth",
    "unique_mappability",
    "dXY_site_fraction"
  ) := NULL
]


# ------------------------------------------------------------
# 7. REORDER COLUMNS
# ------------------------------------------------------------

setcolorder(
  supp_terminal_table,
  c(
    "terminal_band",
    "n_chromosomes",
    "Africa_pi",
    "Asia_pi",
    "Europe_pi",
    "Human_Pig_dXY",
    "median_depth",
    "coverage_breadth_percent",
    "unique_mappability_percent",
    "dXY_site_fraction_percent"
  )
)


# ------------------------------------------------------------
# 8. ROUND FOR REPORTING
# ------------------------------------------------------------

supp_terminal_table[
  ,
  `:=`(
    Africa_pi =
      round(Africa_pi, 5),
    
    Asia_pi =
      round(Asia_pi, 5),
    
    Europe_pi =
      round(Europe_pi, 5),
    
    Human_Pig_dXY =
      round(Human_Pig_dXY, 5),
    
    median_depth =
      round(median_depth, 2),
    
    coverage_breadth_percent =
      round(coverage_breadth_percent, 1),
    
    unique_mappability_percent =
      round(unique_mappability_percent, 1),
    
    dXY_site_fraction_percent =
      round(dXY_site_fraction_percent, 1)
  )
]


# ------------------------------------------------------------
# 9. GIVE COLUMNS PUBLICATION-FRIENDLY NAMES
# ------------------------------------------------------------

setnames(
  supp_terminal_table,
  old = c(
    "terminal_band",
    "n_chromosomes",
    "Africa_pi",
    "Asia_pi",
    "Europe_pi",
    "Human_Pig_dXY",
    "median_depth",
    "coverage_breadth_percent",
    "unique_mappability_percent",
    "dXY_site_fraction_percent"
  ),
  new = c(
    "Distance from chromosome end",
    "n chromosomes",
    "Africa pi",
    "Asia pi",
    "Europe pi",
    "Human-pig dXY",
    "Median depth (x)",
    "Coverage breadth (%)",
    "Unique mappability (%)",
    "dXY site fraction (%)"
  )
)


# ------------------------------------------------------------
# 10. VIEW FINAL TABLE
# ------------------------------------------------------------

print(
  supp_terminal_table
)


# ------------------------------------------------------------
# 11. SAVE
# ------------------------------------------------------------

fwrite(
  supp_terminal_table,
  file.path(
    OUT_DIR,
    "Table_S_terminal_diversity_QC.tsv"
  ),
  sep = "\t"
)

fwrite(
  supp_terminal_table,
  file.path(
    OUT_DIR,
    "Table_S_terminal_diversity_QC.csv"
  )
)
```

---

## Outputs and QC checks

50-kb joint QC summaries and spatial comparisons. Consult recorded logs and intermediate tables to distinguish expected settings from actual execution success.

## Provenance and limitations

- Relies on the historical PIXY results and should be revisited after future PIXY rerun.
- Sample identifiers and host/country metadata were manually curated in `sample_QC_list.csv`; it is not generated by an algorithm in these files.
- Commands, versions, seeds and paths are retained where present; annotations are not evidence that all jobs have been freshly rerun.
- For the discovery and correction of the PIXY all-sites input issue, see [`99_PIXY_allsites_correction.md`](99_PIXY_allsites_correction.md).
