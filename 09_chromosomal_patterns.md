# 09 — Chromosome-scale pi and Tajima’s D

## Purpose

Chromosome-scale pi and Tajima’s D. This file archives the commands and analysis notes from the thesis/paper workflow.

## Inputs

10-kb regional PIXY output windows, chromosome lengths. Source filenames and absolute HPC paths are retained in the original script blocks below.

## Workflow

See the original annotated analysis sequence below. Distinct script blocks are independent jobs/R sessions, **not one executable script**.

## Annotated scripts

> **Historical PIXY-derived analysis:** This uses historical diversity/differentiation outputs, not the corrected PIXY VCF.

```r
#========================================================
# CHROMOSOME-LEVEL PI, TAJIMA'S D, AND
# CHROMOSOME-SIZE CORRELATION
#========================================================
#
# Purpose:
# 1. Read 10 kb windowed pi and Tajima's D files for
#    Africa, Asia and Europe.
# 2. Check for invalid rows, unexpected chromosomes,
#    missing chromosomes, duplicated windows and
#    overlapping windows.
# 3. Summarise each statistic by chromosome and region.
# 4. Calculate ordinary and window-length-weighted
#    chromosome-level means.
# 5. Join the true lengths of the 19 autosomes.
# 6. Test whether chromosome-level weighted mean pi is
#    associated with chromosome length using Spearman's
#    rank correlation.
# 7. Produce regional descriptive summaries, Table S11,
#    and thesis-ready plots.
#
#========================================================

library(tidyverse)

#========================================================
# 1. FILE PATHS
#========================================================

base_dir <- "S:/Research/Lauren_PopGen/pi TajimasD/text files"

outdir <- file.path(
  base_dir,
  "chromosome_level_outputs"
)

dir.create(
  outdir,
  showWarnings = TRUE,
  recursive = TRUE
)

if (!dir.exists(outdir)) {
  stop(
    "Output directory could not be created:\n",
    outdir
  )
}

#========================================================
# 2. TRUE AUTOSOME LENGTHS AND CHROMOSOME LABELS
#========================================================

chr_lengths <- tibble(
  chr_label = paste0("chr", 1:19),

  chromosome = c(
    "CM024168.1", "CM024169.1", "CM024170.1", "CM024171.1",
    "CM024172.1", "CM024173.1", "CM024174.1", "CM024175.1",
    "CM024176.1", "CM024177.1", "CM024178.1", "CM024179.1",
    "CM024180.1", "CM024181.1", "CM024182.1", "CM024183.1",
    "CM024184.1", "CM024185.1", "CM024186.1"
  ),

  chr_length_bp = c(
    18835767, 14205397, 14128337, 13918991, 13427824,
    13248486, 13137961, 12190940, 11532163, 11047307,
    10268763, 9833456, 8681520, 8487591, 8253480,
    8141563, 7348195, 7084669, 7082418
  )
)

if (
  nrow(chr_lengths) != 19 ||
  anyDuplicated(chr_lengths$chromosome) > 0 ||
  anyDuplicated(chr_lengths$chr_label) > 0
) {
  stop(
    "The chromosome metadata must contain 19 unique ",
    "autosomes and 19 unique chromosome labels."
  )
}

#========================================================
# 3. HELPER FUNCTIONS
#========================================================

# Read a windowed statistic file and stop if the file is
# missing or contains no rows.

read_window_table <- function(file) {

  if (!file.exists(file)) {
    stop(
      "Missing input file:\n",
      file
    )
  }

  dat <- read_tsv(
    file,
    show_col_types = FALSE,
    progress = FALSE
  )

  if (nrow(dat) == 0) {
    stop(
      "Input file contains no rows:\n",
      file
    )
  }

  dat
}


# Standardise the region-specific files into a common
# structure.
#
# Rows lacking finite coordinates or statistic values are
# removed before chromosome-level calculations.
#
# Duplicate or overlapping windows cause the script to stop
# so that they are not silently included or discarded.

standardise_windows <- function(
    dat,
    population,
    value_col,
    statistic_name
) {

  required_cols <- c(
    "chromosome",
    "window_pos_1",
    "window_pos_2",
    value_col
  )

  missing_cols <- setdiff(
    required_cols,
    colnames(dat)
  )

  if (length(missing_cols) > 0) {
    stop(
      "Missing columns for ",
      statistic_name,
      " in ",
      population,
      ": ",
      paste(missing_cols, collapse = ", "),
      "\nAvailable columns: ",
      paste(colnames(dat), collapse = ", ")
    )
  }

  dat <- dat %>%
    transmute(
      population = population,
      statistic = statistic_name,

      chromosome = trimws(
        as.character(chromosome)
      ),

      window_pos_1 = suppressWarnings(
        as.numeric(window_pos_1)
      ),

      window_pos_2 = suppressWarnings(
        as.numeric(window_pos_2)
      ),

      value = suppressWarnings(
        as.numeric(.data[[value_col]])
      ),

      window_length_bp =
        window_pos_2 - window_pos_1 + 1
    )

  #------------------------------------------------------
  # Remove invalid rows
  #------------------------------------------------------

  invalid_rows <- dat %>%
    filter(
      is.na(chromosome) |
        chromosome == "" |
        !is.finite(window_pos_1) |
        !is.finite(window_pos_2) |
        !is.finite(value) |
        !is.finite(window_length_bp) |
        window_pos_1 < 1 |
        window_pos_2 < window_pos_1 |
        window_length_bp <= 0
    )

  cat(
    population,
    statistic_name,
    "- rows read:",
    nrow(dat),
    "- invalid rows removed:",
    nrow(invalid_rows),
    "\n"
  )

  dat <- dat %>%
    filter(
      !is.na(chromosome),
      chromosome != "",
      is.finite(window_pos_1),
      is.finite(window_pos_2),
      is.finite(value),
      is.finite(window_length_bp),
      window_pos_1 >= 1,
      window_pos_2 >= window_pos_1,
      window_length_bp > 0
    )

  #------------------------------------------------------
  # Check chromosome names
  #------------------------------------------------------

  unexpected_chromosomes <- setdiff(
    unique(dat$chromosome),
    chr_lengths$chromosome
  )

  if (length(unexpected_chromosomes) > 0) {
    stop(
      population,
      " ",
      statistic_name,
      " contains unexpected chromosomes: ",
      paste(unexpected_chromosomes, collapse = ", ")
    )
  }

  missing_chromosomes <- setdiff(
    chr_lengths$chromosome,
    unique(dat$chromosome)
  )

  if (length(missing_chromosomes) > 0) {
    stop(
      population,
      " ",
      statistic_name,
      " is missing chromosomes: ",
      paste(missing_chromosomes, collapse = ", ")
    )
  }

  #------------------------------------------------------
  # Check duplicated window coordinates
  #------------------------------------------------------

  duplicate_windows <- dat %>%
    count(
      chromosome,
      window_pos_1,
      window_pos_2,
      name = "n"
    ) %>%
    filter(
      n > 1
    )

  if (nrow(duplicate_windows) > 0) {

    cat(
      "\nDuplicated windows found in ",
      population,
      " ",
      statistic_name,
      ":\n",
      sep = ""
    )

    print(
      duplicate_windows,
      n = Inf
    )

    stop(
      population,
      " ",
      statistic_name,
      " contains duplicated window coordinates."
    )
  }

  #------------------------------------------------------
  # Check overlapping windows
  #------------------------------------------------------

  overlapping_windows <- dat %>%
    arrange(
      chromosome,
      window_pos_1,
      window_pos_2
    ) %>%
    group_by(
      chromosome
    ) %>%
    mutate(
      maximum_previous_end = lag(
        cummax(window_pos_2)
      )
    ) %>%
    ungroup() %>%
    filter(
      !is.na(maximum_previous_end),
      window_pos_1 <= maximum_previous_end
    )

  if (nrow(overlapping_windows) > 0) {

    cat(
      "\nOverlapping windows found in ",
      population,
      " ",
      statistic_name,
      ":\n",
      sep = ""
    )

    print(
      overlapping_windows,
      n = 20
    )

    stop(
      population,
      " ",
      statistic_name,
      " contains overlapping windows."
    )
  }

  cat(
    population,
    statistic_name,
    "- valid windows:",
    nrow(dat),
    "- chromosomes:",
    n_distinct(dat$chromosome),
    "- duplicates: 0",
    "- overlaps: 0\n"
  )

  dat
}


# For each population, statistic and chromosome:
#
# - The ordinary mean treats every retained window equally.
# - The weighted mean accounts for the represented length of
#   each window. This mainly affects shorter terminal windows.
#
# No weighting by chromosome size is used in the subsequent
# Spearman correlation.

summarise_by_chromosome <- function(dat) {

  dat %>%
    left_join(
      chr_lengths,
      by = "chromosome"
    ) %>%
    mutate(
      chr_label = factor(
        chr_label,
        levels = paste0("chr", 1:19)
      )
    ) %>%
    group_by(
      population,
      statistic,
      chromosome,
      chr_label,
      chr_length_bp
    ) %>%
    summarise(
      n_windows = n(),

      total_window_bp = sum(
        window_length_bp,
        na.rm = TRUE
      ),

      mean = mean(
        value,
        na.rm = TRUE
      ),

      weighted_mean = weighted.mean(
        value,
        w = window_length_bp,
        na.rm = TRUE
      ),

      median = median(
        value,
        na.rm = TRUE
      ),

      sd = sd(
        value,
        na.rm = TRUE
      ),

      min = min(
        value,
        na.rm = TRUE
      ),

      max = max(
        value,
        na.rm = TRUE
      ),

      .groups = "drop"
    ) %>%
    arrange(
      population,
      statistic,
      chr_label
    )
}


# Write a CSV and provide a clearer error if the output file
# is open in Excel or the directory is not writable.

safe_write_csv <- function(dat, filename) {

  output_file <- file.path(
    outdir,
    filename
  )

  tryCatch(
    {
      write_csv(
        dat,
        output_file
      )
    },
    error = function(e) {
      stop(
        "Could not write:\n",
        output_file,
        "\nThe file may be open in Excel, or the output ",
        "directory may not be writable.\n",
        conditionMessage(e)
      )
    }
  )

  message(
    "Written: ",
    output_file
  )
}

#========================================================
# 4. LOAD WINDOWED PI FILES
#========================================================

pi_windows <- bind_rows(

  standardise_windows(
    read_window_table(
      file.path(
        base_dir,
        "africa_pi.txt"
      )
    ),
    population = "Africa",
    value_col = "avg_pi",
    statistic_name = "pi"
  ),

  standardise_windows(
    read_window_table(
      file.path(
        base_dir,
        "asia_pi.txt"
      )
    ),
    population = "Asia",
    value_col = "avg_pi",
    statistic_name = "pi"
  ),

  standardise_windows(
    read_window_table(
      file.path(
        base_dir,
        "europe_pi.txt"
      )
    ),
    population = "Europe",
    value_col = "avg_pi",
    statistic_name = "pi"
  )
)

#========================================================
# 5. LOAD WINDOWED TAJIMA'S D FILES
#========================================================

tajd_windows <- bind_rows(

  standardise_windows(
    read_window_table(
      file.path(
        base_dir,
        "africa_tajima_d.txt"
      )
    ),
    population = "Africa",
    value_col = "tajima_d",
    statistic_name = "tajimas_d"
  ),

  standardise_windows(
    read_window_table(
      file.path(
        base_dir,
        "asia_tajima_d.txt"
      )
    ),
    population = "Asia",
    value_col = "tajima_d",
    statistic_name = "tajimas_d"
  ),

  standardise_windows(
    read_window_table(
      file.path(
        base_dir,
        "europe_tajima_d.txt"
      )
    ),
    population = "Europe",
    value_col = "tajima_d",
    statistic_name = "tajimas_d"
  )
)

all_windows <- bind_rows(
  pi_windows,
  tajd_windows
)

cat("\n==============================\n")
cat("WINDOW COUNTS AFTER FILTERING\n")
cat("==============================\n")

print(
  all_windows %>%
    count(
      population,
      statistic
    ),
  n = Inf
)

#========================================================
# 6. CHROMOSOME-LEVEL SUMMARIES
#========================================================

chromosome_summary_long <- summarise_by_chromosome(
  all_windows
)

pi_chr_summary <- chromosome_summary_long %>%
  filter(
    statistic == "pi"
  )

tajd_chr_summary <- chromosome_summary_long %>%
  filter(
    statistic == "tajimas_d"
  )

#--------------------------------------------------------
# Check that each population/statistic produced exactly
# 19 chromosome-level summaries
#--------------------------------------------------------

chromosome_count_check <- chromosome_summary_long %>%
  count(
    population,
    statistic,
    name = "n_chromosomes"
  )

cat("\n==============================\n")
cat("CHROMOSOME COUNT CHECK\n")
cat("==============================\n")

print(
  chromosome_count_check,
  n = Inf
)

if (
  any(
    chromosome_count_check$n_chromosomes != 19
  )
) {
  stop(
    "At least one population/statistic does not have ",
    "exactly 19 chromosome summaries."
  )
}

#--------------------------------------------------------
# Check for duplicated chromosome summaries
#--------------------------------------------------------

duplicate_chr_summaries <- chromosome_summary_long %>%
  count(
    population,
    statistic,
    chromosome,
    name = "n"
  ) %>%
  filter(
    n > 1
  )

if (nrow(duplicate_chr_summaries) > 0) {

  print(
    duplicate_chr_summaries,
    n = Inf
  )

  stop(
    "Duplicate chromosome-level summaries were produced."
  )
}

#--------------------------------------------------------
# Check that chromosome lengths and weighted means are
# finite before running the correlation
#--------------------------------------------------------

invalid_chr_summaries <- chromosome_summary_long %>%
  filter(
    !is.finite(chr_length_bp) |
      !is.finite(weighted_mean)
  )

if (nrow(invalid_chr_summaries) > 0) {

  print(
    invalid_chr_summaries,
    n = Inf
  )

  stop(
    "Non-finite chromosome-level summary values were produced."
  )
}

#========================================================
# 7. SPEARMAN CORRELATION
#========================================================
#
# Each chromosome contributes one observation:
#
# 1. Chromosome length
# 2. Chromosome-level window-length-weighted mean pi
#
# Chromosomes are not additionally weighted in Spearman's
# rank correlation. Each of the 19 chromosomes contributes
# one paired observation.
#
#========================================================

spearman_pi_chr_length <- pi_chr_summary %>%
  group_by(
    population
  ) %>%
  summarise(
    n_chromosomes = n(),

    n_complete_chromosomes = sum(
      is.finite(chr_length_bp) &
        is.finite(weighted_mean)
    ),

    spearman_rho = cor(
      chr_length_bp,
      weighted_mean,
      method = "spearman",
      use = "complete.obs"
    ),

    spearman_p = cor.test(
      chr_length_bp,
      weighted_mean,
      method = "spearman",
      exact = FALSE
    )$p.value,

    .groups = "drop"
  )

if (
  any(
    spearman_pi_chr_length$n_complete_chromosomes != 19
  )
) {
  stop(
    "At least one population has fewer than 19 complete ",
    "chromosome-level observations for the Spearman test."
  )
}

cat("\n==============================\n")
cat("SPEARMAN: CHROMOSOME LENGTH VS WEIGHTED MEAN PI\n")
cat("==============================\n")

print(
  spearman_pi_chr_length,
  n = Inf,
  width = Inf
)

#========================================================
# 8. REGION-LEVEL DESCRIPTIVE SUMMARIES
#========================================================
#
# mean_across_chromosomes:
# Arithmetic mean of the 19 chromosome-level weighted means.
# This is the descriptive value reported in the final row
# of Table S11.
#
# genome_wide_weighted_mean:
# Chromosome-level weighted means combined according to the
# total valid window length represented on each chromosome.
#
# The genome-wide weighted mean therefore gives greater
# contribution to chromosomes containing more represented
# sequence and is separate from the descriptive mean across
# chromosomes.
#
#========================================================

region_summary <- chromosome_summary_long %>%
  group_by(
    population,
    statistic
  ) %>%
  summarise(
    n_chromosomes = n(),

    mean_across_chromosomes = mean(
      weighted_mean,
      na.rm = TRUE
    ),

    genome_wide_weighted_mean = weighted.mean(
      weighted_mean,
      w = total_window_bp,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

cat("\n==============================\n")
cat("REGION-LEVEL SUMMARY VALUES\n")
cat("==============================\n")

print(
  region_summary,
  n = Inf,
  width = Inf
)

#========================================================
# 9. GENERATE TABLE S11
#========================================================
#
# The chromosome rows report the mean, weighted mean,
# median and standard deviation of the retained 10 kb
# windows for each chromosome.
#
# The final row for each region is the arithmetic mean of
# the 19 chromosome-level values in each column. It is a
# descriptive chromosome-level summary rather than a
# genome-wide weighted estimate.
#
#========================================================

table_s11 <- chromosome_summary_long %>%
  select(
    population,
    chr_label,
    chromosome,
    statistic,
    mean,
    weighted_mean,
    median,
    sd
  ) %>%
  pivot_wider(
    names_from = statistic,

    values_from = c(
      mean,
      weighted_mean,
      median,
      sd
    ),

    names_glue = "{statistic}_{.value}"
  ) %>%
  select(
    population,
    chr_label,
    chromosome,

    pi_mean,
    pi_weighted_mean,
    pi_median,
    pi_sd,

    tajimas_d_mean,
    tajimas_d_weighted_mean,
    tajimas_d_median,
    tajimas_d_sd
  ) %>%
  arrange(
    population,
    chr_label
  )

table_s11_region_means <- table_s11 %>%
  group_by(
    population
  ) %>%
  summarise(
    chr_label = "Mean across chromosomes",
    chromosome = NA_character_,

    pi_mean = mean(
      pi_mean,
      na.rm = TRUE
    ),

    pi_weighted_mean = mean(
      pi_weighted_mean,
      na.rm = TRUE
    ),

    pi_median = mean(
      pi_median,
      na.rm = TRUE
    ),

    pi_sd = mean(
      pi_sd,
      na.rm = TRUE
    ),

    tajimas_d_mean = mean(
      tajimas_d_mean,
      na.rm = TRUE
    ),

    tajimas_d_weighted_mean = mean(
      tajimas_d_weighted_mean,
      na.rm = TRUE
    ),

    tajimas_d_median = mean(
      tajimas_d_median,
      na.rm = TRUE
    ),

    tajimas_d_sd = mean(
      tajimas_d_sd,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

table_s11_final <- bind_rows(
  table_s11,
  table_s11_region_means
) %>%
  mutate(
    chr_order = match(
      as.character(chr_label),
      c(
        paste0("chr", 1:19),
        "Mean across chromosomes"
      )
    )
  ) %>%
  arrange(
    population,
    chr_order
  ) %>%
  select(
    -chr_order
  )

cat("\n==============================\n")
cat("TABLE S11\n")
cat("==============================\n")

print(
  table_s11_final,
  n = Inf,
  width = Inf
)

#========================================================
# 10. SAVE OUTPUT TABLES
#========================================================

safe_write_csv(
  chromosome_summary_long,
  "chromosome_level_pi_tajimasd_summary.csv"
)

safe_write_csv(
  spearman_pi_chr_length,
  "chromosome_length_vs_weighted_pi_spearman.csv"
)

safe_write_csv(
  region_summary,
  "region_level_pi_tajimasd_summary.csv"
)

safe_write_csv(
  table_s11_final,
  "Table_S11_chromosome_level_pi_tajimasD.csv"
)

#========================================================
# 11. PLOT: CHROMOSOME SIZE VS WEIGHTED MEAN PI
#========================================================
#
# The fitted linear regression and 95% confidence band are
# included only to visualise the direction of the trend.
#
# The formal association reported in the thesis is the
# Spearman rank correlation calculated above.
#
#========================================================

p_pi_chr_length <- ggplot(
  pi_chr_summary,
  aes(
    x = chr_length_bp / 1e6,
    y = weighted_mean
  )
) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.8,
    colour = "#3366FF",
    fill = "grey70",
    alpha = 0.55
  ) +
  geom_point(
    size = 1.6,
    colour = "black"
  ) +
  facet_wrap(
    ~ population,
    nrow = 1,
    scales = "fixed"
  ) +
  scale_x_continuous(
    breaks = c(
      10,
      15
    ),
    limits = c(
      6.5,
      19.5
    ),
    expand = expansion(
      mult = c(
        0.01,
        0.01
      )
    )
  ) +
  scale_y_continuous(
    breaks = seq(
      0.006,
      0.014,
      by = 0.002
    ),
    limits = c(
      0.0058,
      0.0144
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  labs(
    title = "Relationship between chromosome size and nucleotide diversity",
    x = "Chromosome size (Mb)",
    y = expression(
      "Weighted mean nucleotide diversity (" * pi * ")"
    )
  ) +
  theme_bw(
    base_size = 14
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      face = "plain",
      hjust = 0
    ),

    axis.title = element_text(
      size = 16
    ),

    axis.text = element_text(
      size = 13,
      colour = "grey30"
    ),

    strip.background = element_rect(
      fill = "grey85",
      colour = "grey30",
      linewidth = 0.7
    ),

    strip.text = element_text(
      size = 15,
      face = "plain"
    ),

    panel.grid.major = element_line(
      colour = "grey88",
      linewidth = 0.6
    ),

    panel.grid.minor = element_line(
      colour = "grey93",
      linewidth = 0.4
    ),

    panel.spacing = grid::unit(
      0.15,
      "cm"
    ),

    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  )

print(
  p_pi_chr_length
)

ggsave(
  filename = file.path(
    outdir,
    "chromosome_length_vs_weighted_pi.png"
  ),
  plot = p_pi_chr_length,
  width = 14,
  height = 8,
  dpi = 600
)

ggsave(
  filename = file.path(
    outdir,
    "chromosome_length_vs_weighted_pi.svg"
  ),
  plot = p_pi_chr_length,
  width = 14,
  height = 8
)

#========================================================
# 12. PLOT: WEIGHTED PI BY CHROMOSOME
#========================================================

p_pi_by_chromosome <- ggplot(
  pi_chr_summary,
  aes(
    x = chr_label,
    y = weighted_mean
  )
) +
  geom_col() +
  facet_wrap(
    ~ population,
    scales = "free_y"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    title = "Chromosome-level weighted mean nucleotide diversity",
    x = "Chromosome",
    y = expression(
      "Weighted mean " * pi
    )
  )

print(
  p_pi_by_chromosome
)

ggsave(
  filename = file.path(
    outdir,
    "chromosome_level_weighted_pi_by_region.png"
  ),
  plot = p_pi_by_chromosome,
  width = 10,
  height = 5.5,
  dpi = 600
)

#========================================================
# 13. SCRIPT COMPLETE
#========================================================

cat("\n==============================\n")
cat("SCRIPT COMPLETE\n")
cat("==============================\n")
cat("Outputs written to:\n")
cat(outdir, "\n")
```

---

## Outputs and QC checks

Chromosome summaries, chromosome-size Spearman correlation, tables. Consult recorded logs and intermediate tables to distinguish expected settings from actual execution success.

## Provenance and limitations

- These summarize historical PIXY-derived values; affected branch requires reanalysis for publication.
- Sample identifiers and host/country metadata were manually curated in `sample_QC_list.csv`; it is not generated by an algorithm in these files.
- Commands, versions, seeds and paths are retained where present; annotations are not evidence that all jobs have been freshly rerun.
- For the discovery and correction of the PIXY all-sites input issue, see [`99_PIXY_allsites_correction.md`](99_PIXY_allsites_correction.md).
