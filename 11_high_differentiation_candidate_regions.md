# 11 — Empirical genomic differentiation outliers

## Purpose

Empirical genomic differentiation outliers. This file archives the commands and analysis notes from the thesis/paper workflow.

## Inputs

10-kb PIXY FST/dXY windows, gene annotations. Source filenames and absolute HPC paths are retained in the original script blocks below.

## Workflow

See the original annotated analysis sequence below. Distinct script blocks are independent jobs/R sessions, **not one executable script**.

## Annotated scripts

> **Rerun note:** The first candidate-region block calls `file.path(outdir, ...)` before defining `outdir`. In the historical RStudio session this object may already have existed. To rerun the block independently, the analyst must set/verify `outdir` first, without treating the directory value as an original method. All candidate estimates are downstream of the historical PIXY input.


This annotated section documents the exploratory candidate-region workflow. It is separated from the genome-wide summary section because the candidate-region scripts use empirical outlier thresholds and focus on two focal comparisons:

1. Africa vs non-Africa
2. Human vs pig

The aim is to identify windows and regions showing high relative differentiation (FST), high absolute divergence (dXY), or both, and then connect those regions to nearby or overlapping genes.

## Overall rationale

- High FST alone can reflect reduced diversity in one group as well as true divergence.
- High dXY helps identify regions with elevated absolute divergence.
- Overlap between high FST and high dXY is therefore more conservative than using FST alone.
- Nearby windows are also considered because linked or regulatory variation may not fall in the exact same 10 kb window.
- Candidate genes are interpreted cautiously because host and geography are confounded in the sampling design.



## Africa vs non-Africa candidate-region and gene annotation code


```r

#========================================================
# PURPOSE OF THIS SCRIPT
#========================================================
# This script:
# 1. Reads in windowed FST and dXY results
# 2. Keeps only the focal comparison: Africa vs non-Africa
# 3. Removes duplicated windows
# 4. Calculates bootstrap confidence intervals (CIs) for
#    the genome-wide mean FST and genome-wide mean dXY
# 5. Identifies top 1% and top 5% outlier windows
# 6. Tests whether FST and dXY outliers overlap exactly
# 7. Tests whether FST and dXY outliers cluster on the
#    same chromosomes
# 8. Tests whether FST and dXY outliers occur nearby
#    on the same chromosome even if they are not the
#    exact same window
# 9. Merges shared outlier windows into candidate regions
# 10. Intersects candidate regions with annotated genes
# 11. Uses tx61820 transcript links to recover additional
#     transcript-level annotation where possible
# 12. Builds a final cleaned candidate gene summary table
#     for figures and tables
#========================================================

library(tidyverse)
library(readxl)

set.seed(123)

#========================================================
# 1. READ INPUT FILES
#========================================================

fst <- read_tsv("africa_vs_nonafrica_fst.txt", show_col_types = FALSE)
dxy <- read_tsv("africa_vs_nonafrica_dxy.txt", show_col_types = FALSE)

#========================================================
# 2. INSPECT WHICH COMPARISONS ARE PRESENT
#========================================================

cat("\n==============================\n")
cat("UNIQUE FST COMPARISONS\n")
cat("==============================\n")
print(fst %>% distinct(pop1, pop2), n = Inf)

cat("\n==============================\n")
cat("UNIQUE dXY COMPARISONS\n")
cat("==============================\n")
print(dxy %>% distinct(pop1, pop2), n = Inf)

#========================================================
# 3. DEFINE THE VALUE COLUMNS
#========================================================

fst_col <- "avg_wc_fst"
dxy_col <- "avg_dxy"

#========================================================
# 4. KEEP ONLY THE FOCAL COMPARISON
#========================================================

fst <- fst %>%
  filter(pop1 == "Africa", pop2 == "non-Africa")

dxy <- dxy %>%
  filter(pop1 == "Africa", pop2 == "non-Africa")

#========================================================
# 5. REMOVE ANY DUPLICATED WINDOWS
#========================================================

fst <- fst %>%
  distinct(chromosome, window_pos_1, window_pos_2, .keep_all = TRUE)

dxy <- dxy %>%
  distinct(chromosome, window_pos_1, window_pos_2, .keep_all = TRUE)

#========================================================
# 6. BASIC SANITY CHECKS
#========================================================

cat("\n==============================\n")
cat("CHECKS AFTER FILTERING\n")
cat("==============================\n")
cat("FST rows:", nrow(fst), "\n")
cat("dXY rows:", nrow(dxy), "\n")

cat("\nUnique windows:\n")
cat("FST unique windows:", n_distinct(fst$chromosome, fst$window_pos_1, fst$window_pos_2), "\n")
cat("dXY unique windows:", n_distinct(dxy$chromosome, dxy$window_pos_1, dxy$window_pos_2), "\n")

#========================================================
# 7. HELPER FUNCTION: BOOTSTRAP CI FOR MEAN
#========================================================

boot_mean_ci <- function(x, nboot = 1000) {
  x <- x[is.finite(x)]
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    stop("No valid values found after filtering.")
  }
  
  boots <- replicate(nboot, {
    mean(sample(x, size = length(x), replace = TRUE))
  })
  
  tibble(
    n_windows = length(x),
    mean      = mean(x),
    median    = median(x),
    sd        = sd(x),
    min       = min(x),
    max       = max(x),
    ci_lower  = unname(quantile(boots, 0.025, na.rm = TRUE)),
    ci_upper  = unname(quantile(boots, 0.975, na.rm = TRUE))
  )
}

#========================================================
# 8. RUN BOOTSTRAP FOR GENOME-WIDE FST AND dXY
#========================================================

fst_summary <- boot_mean_ci(fst[[fst_col]], nboot = 1000) %>%
  mutate(statistic = "FST")

dxy_summary <- boot_mean_ci(dxy[[dxy_col]], nboot = 1000) %>%
  mutate(statistic = "dXY")

summary_table <- bind_rows(fst_summary, dxy_summary) %>%
  select(statistic, n_windows, mean, median, sd, min, max, ci_lower, ci_upper)

cat("\n==============================\n")
cat("BOOTSTRAP SUMMARY TABLE\n")
cat("==============================\n")
print(as.data.frame(summary_table))

ci_table <- summary_table %>%
  select(statistic, mean, ci_lower, ci_upper)

cat("\n==============================\n")
cat("BOOTSTRAP CI TABLE ONLY\n")
cat("==============================\n")
print(as.data.frame(ci_table))

write_csv(summary_table, file.path(outdir, "africa_vs_nonafrica_bootstrap_summary.csv"))
write_csv(ci_table, file.path(outdir, "africa_vs_nonafrica_bootstrap_ci_only.csv"))

#========================================================
# 9. PREPARE WINDOWED DATA FOR OUTLIER ANALYSIS
#========================================================

fst2 <- fst %>%
  filter(is.finite(.data[[fst_col]]), !is.na(.data[[fst_col]])) %>%
  mutate(midpoint = (window_pos_1 + window_pos_2) / 2)

dxy2 <- dxy %>%
  filter(is.finite(.data[[dxy_col]]), !is.na(.data[[dxy_col]])) %>%
  mutate(midpoint = (window_pos_1 + window_pos_2) / 2)

#========================================================
# 10. TOP 1% OUTLIERS
#========================================================

fst_cutoff_1 <- quantile(fst2[[fst_col]], 0.99, na.rm = TRUE)
dxy_cutoff_1 <- quantile(dxy2[[dxy_col]], 0.99, na.rm = TRUE)

cat("\n==============================\n")
cat("TOP 1% CUTOFFS\n")
cat("==============================\n")
cat("FST cutoff:", fst_cutoff_1, "\n")
cat("dXY cutoff:", dxy_cutoff_1, "\n")

fst_outliers_1 <- fst2 %>%
  filter(.data[[fst_col]] >= fst_cutoff_1) %>%
  arrange(desc(.data[[fst_col]])) %>%
  mutate(comparison = "Africa_vs_nonAfrica") %>%
  select(comparison, chromosome, window_pos_1, window_pos_2, midpoint, all_of(fst_col))

dxy_outliers_1 <- dxy2 %>%
  filter(.data[[dxy_col]] >= dxy_cutoff_1) %>%
  arrange(desc(.data[[dxy_col]])) %>%
  mutate(comparison = "Africa_vs_nonAfrica") %>%
  select(comparison, chromosome, window_pos_1, window_pos_2, midpoint, all_of(dxy_col))

#========================================================
# 11. EXACT OVERLAP OF TOP 1% OUTLIERS
#========================================================

shared_outliers_1 <- fst_outliers_1 %>%
  inner_join(
    dxy_outliers_1,
    by = c("comparison", "chromosome", "window_pos_1", "window_pos_2", "midpoint")
  ) %>%
  arrange(desc(avg_wc_fst), desc(avg_dxy))

cat("\n==============================\n")
cat("TOP 1% FST OUTLIERS\n")
cat("==============================\n")
print(fst_outliers_1, n = 20)

cat("\n==============================\n")
cat("TOP 1% dXY OUTLIERS\n")
cat("==============================\n")
print(dxy_outliers_1, n = 20)

cat("\n==============================\n")
cat("SHARED TOP 1% WINDOWS\n")
cat("==============================\n")
print(shared_outliers_1, n = 50)
cat("Number of shared top 1% windows:", nrow(shared_outliers_1), "\n")

write_csv(fst_outliers_1, file.path(outdir, "africa_vs_nonafrica_fst_top1pct_windows.csv"))
write_csv(dxy_outliers_1, file.path(outdir, "africa_vs_nonafrica_dxy_top1pct_windows.csv"))
write_csv(shared_outliers_1, file.path(outdir, "africa_vs_nonafrica_shared_top1pct_windows.csv"))

#========================================================
# 12. TOP 5% OUTLIERS
#========================================================

fst_cutoff_5 <- quantile(fst2[[fst_col]], 0.95, na.rm = TRUE)
dxy_cutoff_5 <- quantile(dxy2[[dxy_col]], 0.95, na.rm = TRUE)

cat("\n==============================\n")
cat("TOP 5% CUTOFFS\n")
cat("==============================\n")
cat("FST cutoff:", fst_cutoff_5, "\n")
cat("dXY cutoff:", dxy_cutoff_5, "\n")

fst_outliers_5 <- fst2 %>%
  filter(.data[[fst_col]] >= fst_cutoff_5) %>%
  arrange(desc(.data[[fst_col]])) %>%
  mutate(comparison = "Africa_vs_nonAfrica") %>%
  select(comparison, chromosome, window_pos_1, window_pos_2, midpoint, all_of(fst_col))

dxy_outliers_5 <- dxy2 %>%
  filter(.data[[dxy_col]] >= dxy_cutoff_5) %>%
  arrange(desc(.data[[dxy_col]])) %>%
  mutate(comparison = "Africa_vs_nonAfrica") %>%
  select(comparison, chromosome, window_pos_1, window_pos_2, midpoint, all_of(dxy_col))

#========================================================
# 13. EXACT OVERLAP OF TOP 5% OUTLIERS
#========================================================

shared_outliers_5 <- fst_outliers_5 %>%
  inner_join(
    dxy_outliers_5,
    by = c("comparison", "chromosome", "window_pos_1", "window_pos_2", "midpoint")
  ) %>%
  arrange(desc(avg_wc_fst), desc(avg_dxy))

cat("\n==============================\n")
cat("SHARED TOP 5% WINDOWS\n")
cat("==============================\n")
print(shared_outliers_5, n = 50)
cat("Number of shared top 5% windows:", nrow(shared_outliers_5), "\n")

write_csv(fst_outliers_5, file.path(outdir, "africa_vs_nonafrica_fst_top5pct_windows.csv"))
write_csv(dxy_outliers_5, file.path(outdir, "africa_vs_nonafrica_dxy_top5pct_windows.csv"))
write_csv(shared_outliers_5, file.path(outdir, "africa_vs_nonafrica_shared_top5pct_windows.csv"))

#========================================================
# 14. CHROMOSOME-LEVEL CLUSTERING: TOP 1%
#========================================================

fst_chr_counts_1 <- fst_outliers_1 %>%
  count(chromosome, sort = TRUE) %>%
  rename(fst_top1_count = n)

dxy_chr_counts_1 <- dxy_outliers_1 %>%
  count(chromosome, sort = TRUE) %>%
  rename(dxy_top1_count = n)

chr_overlap_1 <- full_join(fst_chr_counts_1, dxy_chr_counts_1, by = "chromosome") %>%
  mutate(
    fst_top1_count = replace_na(fst_top1_count, 0),
    dxy_top1_count = replace_na(dxy_top1_count, 0),
    shared_chromosome = fst_top1_count > 0 & dxy_top1_count > 0
  ) %>%
  arrange(desc(shared_chromosome), desc(fst_top1_count + dxy_top1_count))

cat("\n==============================\n")
cat("CHROMOSOME OVERLAP FOR TOP 1% OUTLIERS\n")
cat("==============================\n")
print(chr_overlap_1, n = Inf)

write_csv(fst_chr_counts_1, file.path(outdir, "africa_vs_nonafrica_fst_top1pct_counts_by_chromosome.csv"))
write_csv(dxy_chr_counts_1, file.path(outdir, "africa_vs_nonafrica_dxy_top1pct_counts_by_chromosome.csv"))
write_csv(chr_overlap_1, file.path(outdir, "africa_vs_nonafrica_top1pct_chromosome_overlap.csv"))

#========================================================
# 15. CHROMOSOME-LEVEL CLUSTERING: TOP 5%
#========================================================

fst_chr_counts_5 <- fst_outliers_5 %>%
  count(chromosome, sort = TRUE) %>%
  rename(fst_top5_count = n)

dxy_chr_counts_5 <- dxy_outliers_5 %>%
  count(chromosome, sort = TRUE) %>%
  rename(dxy_top5_count = n)

chr_overlap_5 <- full_join(fst_chr_counts_5, dxy_chr_counts_5, by = "chromosome") %>%
  mutate(
    fst_top5_count = replace_na(fst_top5_count, 0),
    dxy_top5_count = replace_na(dxy_top5_count, 0),
    shared_chromosome = fst_top5_count > 0 & dxy_top5_count > 0
  ) %>%
  arrange(desc(shared_chromosome), desc(fst_top5_count + dxy_top5_count))

cat("\n==============================\n")
cat("CHROMOSOME OVERLAP FOR TOP 5% OUTLIERS\n")
cat("==============================\n")
print(chr_overlap_5, n = Inf)

write_csv(fst_chr_counts_5, file.path(outdir, "africa_vs_nonafrica_fst_top5pct_counts_by_chromosome.csv"))
write_csv(dxy_chr_counts_5, file.path(outdir, "africa_vs_nonafrica_dxy_top5pct_counts_by_chromosome.csv"))
write_csv(chr_overlap_5, file.path(outdir, "africa_vs_nonafrica_top5pct_chromosome_overlap.csv"))

#========================================================
# 16. FUNCTION TO TEST NEARBY OVERLAPS
#========================================================

find_nearby_overlaps <- function(fst_df, dxy_df, max_dist = 10000) {
  inner_join(
    fst_df %>% rename(fst_midpoint = midpoint, fst_value = avg_wc_fst),
    dxy_df %>% rename(dxy_midpoint = midpoint, dxy_value = avg_dxy),
    by = c("comparison", "chromosome"),
    relationship = "many-to-many"
  ) %>%
    mutate(distance = abs(fst_midpoint - dxy_midpoint)) %>%
    filter(distance <= max_dist) %>%
    arrange(chromosome, distance, desc(fst_value), desc(dxy_value))
}

#========================================================
# 17. NEARBY OVERLAP: TOP 1%
#========================================================

nearby_overlap_1_10kb <- find_nearby_overlaps(fst_outliers_1, dxy_outliers_1, max_dist = 10000)
nearby_overlap_1_50kb <- find_nearby_overlaps(fst_outliers_1, dxy_outliers_1, max_dist = 50000)

cat("\n==============================\n")
cat("NEARBY OVERLAP FOR TOP 1% OUTLIERS (WITHIN 10 kb)\n")
cat("==============================\n")
print(nearby_overlap_1_10kb, n = 50)
cat("Number of nearby top 1% overlaps within 10 kb:", nrow(nearby_overlap_1_10kb), "\n")

cat("\n==============================\n")
cat("NEARBY OVERLAP FOR TOP 1% OUTLIERS (WITHIN 50 kb)\n")
cat("==============================\n")
print(nearby_overlap_1_50kb, n = 50)
cat("Number of nearby top 1% overlaps within 50 kb:", nrow(nearby_overlap_1_50kb), "\n")

write_csv(nearby_overlap_1_10kb, file.path(outdir, "africa_vs_nonafrica_nearby_overlap_top1pct_10kb.csv"))
write_csv(nearby_overlap_1_50kb, file.path(outdir, "africa_vs_nonafrica_nearby_overlap_top1pct_50kb.csv"))

#========================================================
# 18. NEARBY OVERLAP: TOP 5%
#========================================================

nearby_overlap_5_10kb <- find_nearby_overlaps(fst_outliers_5, dxy_outliers_5, max_dist = 10000)
nearby_overlap_5_50kb <- find_nearby_overlaps(fst_outliers_5, dxy_outliers_5, max_dist = 50000)

cat("\n==============================\n")
cat("NEARBY OVERLAP FOR TOP 5% OUTLIERS (WITHIN 10 kb)\n")
cat("==============================\n")
print(nearby_overlap_5_10kb, n = 50)
cat("Number of nearby top 5% overlaps within 10 kb:", nrow(nearby_overlap_5_10kb), "\n")

cat("\n==============================\n")
cat("NEARBY OVERLAP FOR TOP 5% OUTLIERS (WITHIN 50 kb)\n")
cat("==============================\n")
print(nearby_overlap_5_50kb, n = 50)
cat("Number of nearby top 5% overlaps within 50 kb:", nrow(nearby_overlap_5_50kb), "\n")

write_csv(nearby_overlap_5_10kb, file.path(outdir, "africa_vs_nonafrica_nearby_overlap_top5pct_10kb.csv"))
write_csv(nearby_overlap_5_50kb, file.path(outdir, "africa_vs_nonafrica_nearby_overlap_top5pct_50kb.csv"))

#========================================================
# 19. OPTIONAL: COUNT UNIQUE WINDOWS WITH NEARBY MATCHES
#========================================================

nearby_fst_top1_10kb <- nearby_overlap_1_10kb %>%
  distinct(comparison, chromosome, window_pos_1.x, window_pos_2.x, fst_midpoint, fst_value)

nearby_dxy_top1_10kb <- nearby_overlap_1_10kb %>%
  distinct(comparison, chromosome, window_pos_1.y, window_pos_2.y, dxy_midpoint, dxy_value)

cat("\n==============================\n")
cat("UNIQUE TOP 1% WINDOWS WITH A NEARBY MATCH WITHIN 10 kb\n")
cat("==============================\n")
cat("Unique FST windows with a nearby dXY peak:", nrow(nearby_fst_top1_10kb), "\n")
cat("Unique dXY windows with a nearby FST peak:", nrow(nearby_dxy_top1_10kb), "\n")

write_csv(nearby_fst_top1_10kb, file.path(outdir, "africa_vs_nonafrica_unique_fst_top1pct_with_nearby_dxy_10kb.csv"))
write_csv(nearby_dxy_top1_10kb, file.path(outdir, "africa_vs_nonafrica_unique_dxy_top1pct_with_nearby_fst_10kb.csv"))

#========================================================
# 20. RANK EXACT SHARED TOP 5% WINDOWS
#========================================================

shared_top5_ranked <- shared_outliers_5 %>%
  mutate(
    fst_rank = min_rank(desc(avg_wc_fst)),
    dxy_rank = min_rank(desc(avg_dxy)),
    combined_rank_score = fst_rank + dxy_rank
  ) %>%
  arrange(combined_rank_score, desc(avg_wc_fst), desc(avg_dxy))

cat("\n==============================\n")
cat("RANKED SHARED TOP 5% WINDOWS\n")
cat("==============================\n")
print(shared_top5_ranked, n = 50)

write_csv(shared_top5_ranked, file.path(outdir, "africa_vs_nonafrica_ranked_shared_top5pct_windows.csv"))

#========================================================
# 21. RANK NEARBY TOP 1% OVERLAPS WITHIN 10 kb
#========================================================

nearby_top1_ranked_10kb <- nearby_overlap_1_10kb %>%
  mutate(
    fst_rank = min_rank(desc(fst_value)),
    dxy_rank = min_rank(desc(dxy_value)),
    dist_rank = min_rank(distance),
    combined_rank_score = fst_rank + dxy_rank + dist_rank
  ) %>%
  arrange(combined_rank_score, distance, desc(fst_value), desc(dxy_value))

cat("\n==============================\n")
cat("RANKED NEARBY TOP 1% OVERLAPS (10 kb)\n")
cat("==============================\n")
print(nearby_top1_ranked_10kb, n = 50)

write_csv(nearby_top1_ranked_10kb, file.path(outdir, "africa_vs_nonafrica_ranked_nearby_top1pct_10kb.csv"))

#========================================================
# 22. FUNCTION TO MERGE NEARBY WINDOWS INTO REGIONS
#========================================================

merge_windows_to_regions <- function(df, gap_allowance = 10000) {
  df %>%
    arrange(chromosome, window_pos_1, window_pos_2) %>%
    group_by(chromosome) %>%
    mutate(
      prev_end = lag(window_pos_2),
      new_region = if_else(
        is.na(prev_end) | window_pos_1 > (prev_end + gap_allowance),
        1L, 0L
      ),
      region_id = cumsum(new_region)
    ) %>%
    group_by(chromosome, region_id) %>%
    summarise(
      region_start = min(window_pos_1),
      region_end = max(window_pos_2),
      n_windows = n(),
      mean_fst = mean(avg_wc_fst, na.rm = TRUE),
      max_fst = max(avg_wc_fst, na.rm = TRUE),
      mean_dxy = mean(avg_dxy, na.rm = TRUE),
      max_dxy = max(avg_dxy, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(region_width_bp = region_end - region_start + 1)
}

#========================================================
# 23. CANDIDATE REGIONS FROM SHARED TOP 5% WINDOWS
#========================================================

shared_top5_regions <- merge_windows_to_regions(shared_outliers_5, gap_allowance = 10000) %>%
  arrange(desc(max_fst), desc(max_dxy), desc(n_windows))

cat("\n==============================\n")
cat("CANDIDATE REGIONS FROM SHARED TOP 5% WINDOWS\n")
cat("==============================\n")
print(shared_top5_regions, n = 50)

write_csv(shared_top5_regions, file.path(outdir, "africa_vs_nonafrica_candidate_regions_shared_top5pct.csv"))

#========================================================
# 24. BUILD A SHORTLIST OF TOP CANDIDATE REGIONS
#========================================================

priority_candidates <- shared_top5_regions %>%
  mutate(
    priority_score = rank(-max_fst, ties.method = "min") +
      rank(-max_dxy, ties.method = "min") +
      rank(-n_windows, ties.method = "min")
  ) %>%
  arrange(priority_score, desc(max_fst), desc(max_dxy), desc(n_windows))

cat("\n==============================\n")
cat("PRIORITY CANDIDATE REGIONS\n")
cat("==============================\n")
print(priority_candidates, n = 20)

write_csv(priority_candidates, file.path(outdir, "africa_vs_nonafrica_priority_candidate_regions.csv"))

#========================================================
# 25. LOAD GENE ANNOTATION FILE
#========================================================

genes <- read_excel("Ascaris_genes.15714.annotation_expression.xlsx")

cat("\n==============================\n")
cat("GENE FILE COLUMNS\n")
cat("==============================\n")
print(colnames(genes))

#========================================================
# 26. BUILD CHROMOSOME LOOKUP TABLE
#========================================================

chr_key <- tibble(
  gene_chr = c(
    "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
    "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
    "chrX1","chrX2","chrX3","chrX4","chrX5","chrMT"
  ),
  chromosome = c(
    "CM024168.1","CM024169.1","CM024170.1","CM024171.1","CM024172.1",
    "CM024173.1","CM024174.1","CM024175.1","CM024176.1","CM024177.1",
    "CM024178.1","CM024179.1","CM024180.1","CM024181.1","CM024182.1",
    "CM024183.1","CM024184.1","CM024185.1","CM024186.1",
    "CM024187.1","CM024188.1","CM024189.1","CM024190.1","CM024191.1",
    "CM024213.1"
  )
)

#========================================================
# 27. KEEP AND RENAME RELEVANT GENE COLUMNS
#========================================================

genes2 <- genes %>%
  select(
    chr,
    start,
    end,
    name,
    `Original gene description`,
    `Protein function group name`,
    `Curated annotation`,
    `AL's updated description`,
    `AL's comments`,
    `AL's functional class`,
    ASC10_names,
    tx61820_names,
    `ISO-seq_names`,
    expression,
    enriched_genes,
    average,
    max
  ) %>%
  rename(
    gene_chr   = chr,
    gene_start = start,
    gene_end   = end,
    gene_name  = name
  ) %>%
  left_join(chr_key, by = "gene_chr")

cat("\n==============================\n")
cat("CHROMOSOME NAMING CHECK\n")
cat("==============================\n")
cat("Example chromosomes in genes file:\n")
print(head(unique(genes2$gene_chr), 20))

cat("\nMapped chromosomes in genes file:\n")
print(head(unique(genes2$chromosome), 20))

cat("\nExample chromosomes in candidate regions:\n")
print(head(unique(priority_candidates$chromosome), 20))

cat("\nNumber of genes with missing chromosome mapping:\n")
print(sum(is.na(genes2$chromosome)))

unmapped_gene_chromosomes <- genes2 %>%
  filter(is.na(chromosome)) %>%
  distinct(gene_chr)

cat("\nUnmapped gene chromosome names (if any):\n")
print(unmapped_gene_chromosomes, n = Inf)

#========================================================
# 28. GENES OVERLAPPING PRIORITY CANDIDATE REGIONS
#========================================================

genes_overlapping_priority <- genes2 %>%
  filter(!is.na(chromosome)) %>%
  inner_join(priority_candidates, by = "chromosome") %>%
  filter(
    gene_start <= region_end,
    gene_end >= region_start
  ) %>%
  mutate(
    overlap_type = "overlapping",
    overlap_start = pmax(gene_start, region_start),
    overlap_end   = pmin(gene_end, region_end),
    overlap_bp    = overlap_end - overlap_start + 1
  ) %>%
  arrange(priority_score, chromosome, region_start, gene_start)

cat("\n==============================\n")
cat("GENES OVERLAPPING PRIORITY CANDIDATE REGIONS\n")
cat("==============================\n")
print(genes_overlapping_priority, n = 50)

write_csv(
  genes_overlapping_priority,
  file.path(outdir, "africa_vs_nonafrica_genes_overlapping_priority_regions.csv")
)

#========================================================
# 29. GENES NEAR PRIORITY CANDIDATE REGIONS (WITHIN 10 kb)
#========================================================

genes_near_priority_10kb <- genes2 %>%
  filter(!is.na(chromosome)) %>%
  inner_join(priority_candidates, by = "chromosome") %>%
  mutate(
    distance_to_region = case_when(
      gene_end < region_start ~ region_start - gene_end,
      gene_start > region_end ~ gene_start - region_end,
      TRUE ~ 0
    )
  ) %>%
  filter(distance_to_region <= 10000) %>%
  mutate(
    overlap_type = case_when(
      distance_to_region == 0 ~ "overlapping",
      TRUE ~ "nearby_10kb"
    )
  ) %>%
  arrange(priority_score, distance_to_region, chromosome, region_start, gene_start)

cat("\n==============================\n")
cat("GENES OVERLAPPING OR NEAR PRIORITY REGIONS (<=10 kb)\n")
cat("==============================\n")
print(genes_near_priority_10kb, n = 50)

write_csv(
  genes_near_priority_10kb,
  file.path(outdir, "africa_vs_nonafrica_genes_overlapping_or_near_priority_regions_10kb.csv")
)

#========================================================
# 30. CREATE ONE SUMMARY ROW PER GENE
#========================================================

priority_gene_summary <- genes_near_priority_10kb %>%
  group_by(
    chromosome,
    gene_chr,
    gene_start,
    gene_end,
    gene_name,
    `Original gene description`,
    `Protein function group name`,
    `Curated annotation`,
    `AL's updated description`,
    `AL's comments`,
    `AL's functional class`,
    ASC10_names,
    tx61820_names,
    `ISO-seq_names`,
    expression,
    enriched_genes,
    average,
    max
  ) %>%
  summarise(
    n_matching_regions = n(),
    best_priority_score = min(priority_score, na.rm = TRUE),
    best_max_fst = max(max_fst, na.rm = TRUE),
    best_max_dxy = max(max_dxy, na.rm = TRUE),
    total_windows_in_matched_regions = sum(n_windows, na.rm = TRUE),
    min_distance_to_region = min(distance_to_region, na.rm = TRUE),
    match_type = case_when(
      min(distance_to_region, na.rm = TRUE) == 0 ~ "overlapping",
      TRUE ~ "nearby_10kb"
    ),
    .groups = "drop"
  ) %>%
  arrange(best_priority_score, min_distance_to_region, desc(best_max_fst), desc(best_max_dxy))

cat("\n==============================\n")
cat("PRIORITISED CANDIDATE GENE SUMMARY\n")
cat("==============================\n")
print(priority_gene_summary, n = 50)

write_csv(
  priority_gene_summary,
  file.path(outdir, "africa_vs_nonafrica_priority_candidate_genes_summary.csv")
)

#========================================================
# 31. USE ALL CANDIDATE GENES
#========================================================
# From this point onwards, all plots and summaries use the full
# priority_gene_summary dataset; no candidate-gene ranking cutoff is applied.

all_candidate_genes <- priority_gene_summary %>%
  arrange(best_priority_score, min_distance_to_region,
          desc(best_max_fst), desc(best_max_dxy))

cat("\n==============================\n")
cat("ALL CANDIDATE GENES\n")
cat("==============================\n")
print(all_candidate_genes, n = Inf)

write_csv(
  all_candidate_genes,
  file.path(outdir, "africa_vs_nonafrica_all_candidate_genes.csv")
)

#========================================================
# 32. LOAD TRANSCRIPT ANNOTATION FILE
#========================================================

tx_annot_raw <- read_excel("Ascaris_transcripts_2017_version.xlsx")

cat("\nTRANSCRIPT FILE COLUMNS:\n")
print(colnames(tx_annot_raw))

tx_annot <- tx_annot_raw %>%
  rename_with(~tolower(.)) %>%
  rename(
    transcript_name = transcript_name,
    tx_chromosome = chromosome,
    tx_start = start,
    tx_end = end,
    tx_annotation = annotation
  ) %>%
  mutate(
    transcript_name = trimws(transcript_name)
  )

#========================================================
# 33. EXPAND tx61820 LINKS FOR ALL CANDIDATE GENES
#========================================================

candidate_tx_links <- all_candidate_genes %>%
  select(gene_name, tx61820_names) %>%
  separate_rows(tx61820_names, sep = ",\\s*") %>%
  mutate(tx61820_names = trimws(tx61820_names)) %>%
  filter(!is.na(tx61820_names),
         tx61820_names != "",
         tx61820_names != "NA")

cat("\n==============================\n")
cat("EXPANDED TRANSCRIPT LINKS FOR ALL CANDIDATE GENES\n")
cat("==============================\n")
print(candidate_tx_links %>% count(gene_name), n = Inf)

#========================================================
# 34. JOIN TO TRANSCRIPT ANNOTATION
#========================================================

candidate_tx_annot <- candidate_tx_links %>%
  left_join(
    tx_annot,
    by = c("tx61820_names" = "transcript_name"),
    keep = TRUE
  ) %>%
  mutate(
    tx_gene_id = sub("_t[0-9]+$", "", tx61820_names)
  )

cat("\n==============================\n")
cat("TRANSCRIPT MATCHING SUMMARY")
cat("\n==============================\n")
print(
  candidate_tx_annot %>%
    summarise(
      total_rows = n(),
      matched_rows = sum(!is.na(tx_annotation)),
      unmatched_rows = sum(is.na(tx_annotation)),
      matched_fraction = matched_rows / total_rows
    )
)

write_csv(
  candidate_tx_annot,
  file.path(outdir, "all_candidate_gene_transcript_annotation_long.csv")
)

#========================================================
# 35. SUMMARISE TRANSCRIPT ANNOTATION PER GENE
#========================================================

transcript_annotation_summary <- candidate_tx_annot %>%
  group_by(gene_name) %>%
  summarise(
    n_tx61820_links = n_distinct(tx61820_names),
    n_matched_transcripts = n_distinct(transcript_name[!is.na(transcript_name)]),
    recovered_transcript_names = paste(
      unique(transcript_name[!is.na(transcript_name)]),
      collapse = "; "
    ),
    recovered_annotations = tx_annotation[!is.na(tx_annotation) & tx_annotation != ""] %>%
      unique() %>%
      .[!grepl("hypothetical|uncharacterized|unknown", ., ignore.case = TRUE)] %>%
      { if (length(.) == 0) "hypothetical protein" else paste(., collapse = "; ") },
    .groups = "drop"
  )

write_csv(
  transcript_annotation_summary,
  file.path(outdir, "all_candidate_gene_transcript_annotation_summary.csv")
)

#========================================================
# 36. MERGE TRANSCRIPT ANNOTATION BACK INTO ALL CANDIDATES
#========================================================

all_candidate_genes_enriched <- all_candidate_genes %>%
  left_join(transcript_annotation_summary, by = "gene_name") %>%
  mutate(
    preferred_annotation = case_when(
      !is.na(`Curated annotation`) &
        `Curated annotation` != "" &
        `Curated annotation` != "NA" ~ `Curated annotation`,
      
      !is.na(`AL's updated description`) &
        `AL's updated description` != "" &
        `AL's updated description` != "NA" ~ `AL's updated description`,
      
      !is.na(recovered_annotations) &
        recovered_annotations != "" &
        recovered_annotations != "NA" ~ recovered_annotations,
      
      !is.na(`Original gene description`) &
        `Original gene description` != "" &
        `Original gene description` != "NA" ~ `Original gene description`,
      
      TRUE ~ "Uncharacterised"
    )
  )

cat("\n==============================\n")
cat("ENRICHED ALL CANDIDATE GENES\n")
cat("==============================\n")
print(all_candidate_genes_enriched, n = Inf)

write_csv(
  all_candidate_genes_enriched,
  file.path(outdir, "all_candidate_genes_enriched.csv")
)

#========================================================
# 37. FUNCTIONAL CLASS SUMMARY FOR ALL CANDIDATE GENES
#========================================================

all_candidate_annotation_summary <- all_candidate_genes_enriched %>%
  mutate(
    annotation_class = case_when(
      is.na(preferred_annotation) ~ "Uncharacterised",
      grepl("hypothetical|uncharacterized|unknown", preferred_annotation, ignore.case = TRUE) ~ "Uncharacterised",
      grepl("secreted", `AL's functional class`, ignore.case = TRUE) ~ "Secreted",
      TRUE ~ "Annotated"
    )
  ) %>%
  count(annotation_class)

cat("\n==============================\n")
cat("ANNOTATION SUMMARY FOR ALL CANDIDATE GENES\n")
cat("==============================\n")
print(all_candidate_annotation_summary)

cat("\nTotal genes counted: ", sum(all_candidate_annotation_summary$n), "\n")

write_csv(
  all_candidate_annotation_summary,
  file.path(outdir, "all_candidate_gene_annotation_summary.csv")
)

#========================================================
# 38. PLOT: CHROMOSOME DISTRIBUTION FOR ALL CANDIDATE GENES
#========================================================

p_chr_distribution <- ggplot(
  all_candidate_genes_enriched %>% count(gene_chr),
  aes(x = reorder(gene_chr, n), y = n)
) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Distribution of all candidate genes across chromosomes",
    x = "Chromosome",
    y = "Number of genes"
  )

print(p_chr_distribution)

ggsave(
  filename = file.path(outdir, "all_candidate_genes_chr_distribution.png"),
  plot = p_chr_distribution,
  width = 8,
  height = 5,
  dpi = 600
)

ggsave(
  filename = file.path(outdir, "all_candidate_genes_chr_distribution.svg"),
  plot = p_chr_distribution,
  width = 8,
  height = 5
)

#========================================================
# 39. PLOT: DISTANCE TO REGION FOR ALL CANDIDATE GENES
#========================================================

distance_plot_df <- all_candidate_genes_enriched %>%
  mutate(
    distance_category = case_when(
      min_distance_to_region == 0 ~ "Overlapping",
      min_distance_to_region > 0 & min_distance_to_region <= 1000 ~ "0–1 kb",
      min_distance_to_region > 1000 & min_distance_to_region <= 5000 ~ "1–5 kb",
      min_distance_to_region > 5000 & min_distance_to_region <= 10000 ~ "5–10 kb",
      TRUE ~ NA_character_
    ),
    distance_category = factor(
      distance_category,
      levels = c("Overlapping", "0–1 kb", "1–5 kb", "5–10 kb")
    )
  ) %>%
  count(distance_category, .drop = FALSE)

p_distance <- ggplot(distance_plot_df, aes(x = distance_category, y = n)) +
  geom_col(fill = "grey30", width = 0.6) +
  geom_text(aes(label = n), vjust = -0.4, size = 4) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Distance of all candidate genes to prioritised differentiated regions",
    x = "Distance to nearest region",
    y = "Number of genes"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

print(p_distance)

ggsave(
  filename = file.path(outdir, "all_candidate_genes_distance_to_region.png"),
  plot = p_distance,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  filename = file.path(outdir, "all_candidate_genes_distance_to_region.svg"),
  plot = p_distance,
  width = 7,
  height = 5
)

#========================================================
# 40. LOAD GROUPED FUNCTIONAL METADATA
#========================================================
# This file should contain:
# Gene_name and Group

functional_annotation_grouped <- read_excel(
  file.path(base_dir, "functional_annotation_grouped.xlsx")
)

functional_annotation_clean <- functional_annotation_grouped %>%
  group_by(Gene_name) %>%
  summarise(
    Group = case_when(
      any(Group != "Unknown / poorly characterised") ~ first(Group[Group != "Unknown / poorly characterised"]),
      TRUE ~ "Unknown / poorly characterised"
    ),
    .groups = "drop"
  )

#========================================================
# 41. DEFINE SHARED FUNCTIONAL CATEGORY COLOURS
#========================================================
# Defining this once here prevents the error:
# Error: object 'custom_cols' not found

custom_cols <- c(
  "Poorly characterised" = "#fd91a0",
  "Unknown / poorly characterised" = "#fd91a0",
  "Host-parasite interaction / secreted / interface" = "#fdc170",
  "Host–parasite interaction / secreted / interface" = "#fdc170",
  "Metabolism, mitochondrial function, and redox homeostasis" = "#ffec95",
  "Membrane transport and surface-associated proteins" = "#a2daff",
  "Gene expression and genome regulation" = "#c7fa97",
  "Protein synthesis and turnover" = "#dfabff",
  "Signalling and sensory/membrane receptors" = "#8fd3c7",
  "Structural / cell division" = "#ad63ba"
)

custom_cols_chrplot <- custom_cols

#========================================================
# 42. BUILD FINAL CLEAN SUMMARY TABLE FOR ALL CANDIDATE GENES
#========================================================

candidate_gene_summary_table_for_figure <- all_candidate_genes_enriched %>%
  left_join(functional_annotation_clean, by = c("gene_name" = "Gene_name")) %>%
  transmute(
    `Gene Name` = gene_name,
    gene_chr,
    chromosome,
    gene_start,
    gene_end,
    `Max fst` = round(best_max_fst, 4),
    `Max dxy` = round(best_max_dxy, 4),
    `Protein Annotation` = preferred_annotation,
    `Functional Class` = replace_na(`AL's functional class`, "Unknown"),
    Group = replace_na(Group, "Unknown / poorly characterised"),
    best_priority_score,
    min_distance_to_region,
    match_type
  ) %>%
  mutate(
    `Protein Annotation` = case_when(
      `Gene Name` == "ascaris_gene_10238" ~ "Bloom syndrome protein",
      `Gene Name` == "ascaris_gene_04439" ~ "mRNA-decapping enzyme 2",
      `Gene Name` == "ascaris_gene_06316" ~ "AATF protein",
      `Gene Name` == "ascaris_gene_06314" ~ "Pseudouridine-metabolizing protein",
      `Gene Name` == "ascaris_gene_03967" ~ "39S ribosomal protein",
      `Gene Name` == "ascaris_gene_04441" ~ "Uncharacterized protein C28H8.8",
      `Gene Name` %in% c("ascaris_gene_08354", "ascaris_gene_08356") ~ "Serpin-like protein",
      `Gene Name` == "ascaris_gene_10811" ~ "C-type lectin",
      `Gene Name` == "ascaris_gene_06317" ~ "Domain-containing protein",
      `Gene Name` == "ascaris_gene_10841" ~ "Venom allergen 5",
      `Gene Name` == "ascaris_gene_10634" ~ "Lipopolysaccharide-binding protein",
      `Gene Name` == "ascaris_gene_10237" ~ "Ectonucleotide pyrophosphatase/phosphodiesterase",
      `Gene Name` == "ascaris_gene_10546" ~ "Transmembrane protein",
      `Gene Name` == "ascaris_gene_10547" ~ "Transmembrane protein 144",
      `Gene Name` == "ascaris_gene_03968" ~ "Ligand-gated ion channel 50",
      `Gene Name` == "ascaris_gene_03948" ~ "F-box WD repeat-containing protein 7",
      `Gene Name` == "ascaris_gene_03947" ~ "Ubiquinol-cytochrome-c reductase complex assembly factor 2",
      `Gene Name` == "ascaris_gene_11191" ~ "Maltase-glucoamylase-like protein",
      `Gene Name` == "ascaris_gene_00906" ~ "PTCD3 protein",
      `Gene Name` %in% c("ascaris_gene_10636", "ascaris_gene_10637") ~ "Peroxidasin",
      `Gene Name` == "ascaris_gene_06212" ~ "Serpentine receptor class X45",
      `Gene Name` == "ascaris_gene_10548" ~ "Serpentine receptor class delta-19",
      `Gene Name` == "ascaris_gene_03966" ~ "Serine threonine-protein kinase",
      `Gene Name` == "ascaris_gene_11235" ~ "Abnormal spindle protein",
      `Gene Name` == "ascaris_gene_10236" ~ "Prion-like Q/N-rich domain-bearing protein 25",
      `Gene Name` == "ascaris_gene_10635" ~ "GNAT family protein",
      str_detect(str_to_lower(`Protein Annotation`), "hypothetical") ~ "Hypothetical protein",
      TRUE ~ `Protein Annotation`
    ),
    `Functional Class` = case_when(
      `Gene Name` == "ascaris_gene_10238" ~ "Nuclear Regulation",
      `Gene Name` %in% c("ascaris_gene_04439", "ascaris_gene_06316") ~ "Transcription Machinery",
      `Gene Name` == "ascaris_gene_06314" ~ "Unknown",
      `Gene Name` == "ascaris_gene_03967" ~ "Protein Synthesis",
      `Gene Name` == "ascaris_gene_04441" ~ "Secreted / Unknown",
      `Gene Name` %in% c("ascaris_gene_08354", "ascaris_gene_08356") ~ "Secreted / Protease Inhibitor / GPI",
      `Gene Name` == "ascaris_gene_10811" ~ "Immunologically Relevant",
      `Gene Name` == "ascaris_gene_06317" ~ "Secreted",
      `Gene Name` == "ascaris_gene_10841" ~ "Secreted / Unknown",
      `Gene Name` == "ascaris_gene_10634" ~ "Protein Export / GPI",
      `Gene Name` == "ascaris_gene_10237" ~ "Unknown / membrane",
      `Gene Name` == "ascaris_gene_10546" ~ "Transporters / membrane / GPI",
      `Gene Name` %in% c("ascaris_gene_10547", "ascaris_gene_03968") ~ "Transporters / membrane",
      `Gene Name` == "ascaris_gene_03948" ~ "Proteasome",
      `Gene Name` %in% c("ascaris_gene_03947", "ascaris_gene_00906") ~ "Energy Metabolism",
      `Gene Name` %in% c("ascaris_gene_10636", "ascaris_gene_10637") ~ "Detoxification / Oxidation",
      `Gene Name` %in% c("ascaris_gene_06212", "ascaris_gene_10548") ~ "Signal Transduction / membrane",
      `Gene Name` == "ascaris_gene_03966" ~ "Signal Transduction",
      `Gene Name` == "ascaris_gene_11235" ~ "Cytoskeletal",
      str_detect(str_to_lower(`Protein Annotation`), "hypothetical") ~ "Unknown",
      TRUE ~ `Functional Class`
    )
  ) %>%
  arrange(Group, desc(`Max fst`), desc(`Max dxy`))

cat("\n==============================\n")
cat("ALL CANDIDATE GENE SUMMARY TABLE FOR FIGURE\n")
cat("==============================\n")
print(candidate_gene_summary_table_for_figure, n = Inf)

write_csv(
  candidate_gene_summary_table_for_figure,
  file.path(base_dir, "candidate_gene_summary_table_for_figure.csv")
)

#========================================================
# 43. TRUE CHROMOSOME LENGTHS FOR RELATIVE POSITIONING
#========================================================
# These are true assembly chromosome lengths, not estimated
# from the final annotated gene coordinate.

chr_lengths <- tibble(
  gene_chr = c(
    "chr1", "chr2", "chr3", "chr4", "chr5",
    "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19"
  ),
  chromosome = c(
    "CM024168.1", "CM024169.1", "CM024170.1", "CM024171.1", "CM024172.1",
    "CM024173.1", "CM024174.1", "CM024175.1", "CM024176.1", "CM024177.1",
    "CM024178.1", "CM024179.1", "CM024180.1", "CM024181.1", "CM024182.1",
    "CM024183.1", "CM024184.1", "CM024185.1", "CM024186.1"
  ),
  chr_length = c(
    18835767, 14205397, 14128337, 13918991, 13427824,
    13248486, 13137961, 12190940, 11532163, 11047307,
    10268763, 9833456, 8681520, 8487591, 8253480,
    8141563, 7348195, 7084669, 7082418
  )
)

#========================================================
# 44. BUILD GENOMIC PLOTTING TABLE FOR ALL CANDIDATE GENES
#========================================================

candidate_gene_plot_df <- all_candidate_genes_enriched %>%
  left_join(functional_annotation_clean, by = c("gene_name" = "Gene_name")) %>%
  mutate(
    Group = replace_na(Group, "Unknown / poorly characterised"),
    gene_midpoint = (gene_start + gene_end) / 2,
    combined_signal = best_max_fst * best_max_dxy
  ) %>%
  select(
    gene_chr,
    chromosome,
    gene_name,
    gene_start,
    gene_end,
    gene_midpoint,
    Group,
    best_max_fst,
    best_max_dxy,
    combined_signal,
    best_priority_score,
    min_distance_to_region,
    match_type
  ) %>%
  left_join(chr_lengths, by = c("gene_chr", "chromosome")) %>%
  mutate(
    relative_pos = gene_midpoint / chr_length
  )

# Warn if any chromosomes are missing true lengths

missing_chr_lengths <- candidate_gene_plot_df %>%
  filter(is.na(chr_length)) %>%
  distinct(gene_chr, chromosome)

if (nrow(missing_chr_lengths) > 0) {
  cat("\n==============================\n")
  cat("WARNING: SOME CHROMOSOMES ARE MISSING TRUE LENGTHS\n")
  cat("==============================\n")
  print(missing_chr_lengths, n = Inf)
  cat("\nThese genes will have NA relative positions and may not appear in the chromosome plot.\n")
}

chr_levels_present <- c(
  "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
  "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19"
)

candidate_gene_plot_df <- candidate_gene_plot_df %>%
  mutate(
    gene_chr = factor(gene_chr, levels = chr_levels_present)
  ) %>%
  arrange(gene_chr, relative_pos)

write_csv(
  candidate_gene_plot_df,
  file.path(base_dir, "candidate_gene_plot_df.csv")
)

#========================================================
# 45. BUBBLE PLOT: CHROMOSOMAL LOCATION OF ALL CANDIDATE GENES
#========================================================

p_chr_location <- ggplot(
  candidate_gene_plot_df %>% filter(!is.na(relative_pos)),
  aes(
    x = relative_pos,
    y = gene_chr,
    fill = Group,
    size = combined_signal
  )
) +
  geom_point(
    shape = 21,
    colour = "black",
    stroke = 0.35,
    alpha = 1,
    position = position_jitter(width = 0, height = 0.12)
  ) +
  scale_fill_manual(values = custom_cols_chrplot, drop = FALSE) +
  scale_size_continuous(range = c(1, 4)) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0)
  ) +
  labs(
    title = "Chromosomal location of all candidate genes",
    x = "Relative position along chromosome",
    y = "Chromosome",
    fill = "Functional category",
    size = "Combined FST-dXY signal"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p_chr_location)

ggsave(
  filename = file.path(outdir, "all_candidate_genes_chromosomal_location_combined_signal.png"),
  plot = p_chr_location,
  width = 11,
  height = 6.5,
  dpi = 600
)

ggsave(
  filename = file.path(outdir, "all_candidate_genes_chromosomal_location_combined_signal.svg"),
  plot = p_chr_location,
  width = 11,
  height = 6.5
)

#========================================================
# 46. PIE CHART: FUNCTIONAL CATEGORIES OF ALL CANDIDATE GENES
#========================================================

group_counts <- candidate_gene_summary_table_for_figure %>%
  count(Group) %>%
  mutate(
    proportion = n / sum(n),
    percent = round(proportion * 100, 1),
    Group_label = if_else(
      Group == "Unknown / poorly characterised",
      "Poorly characterised",
      Group
    ),
    slice_label = paste0(percent, "%")
  ) %>%
  arrange(desc(n)) %>%
  mutate(Group_label = factor(Group_label, levels = Group_label))

p_functional <- ggplot(group_counts, aes(x = "", y = n, fill = Group_label)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y", start = 0) +
  geom_text(
    aes(label = slice_label),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_fill_manual(values = custom_cols, drop = FALSE) +
  theme_void() +
  labs(
    title = "Functional annotation of all candidate genes",
    fill = "Category"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9)
  )

print(p_functional)

ggsave(
  filename = file.path(outdir, "all_candidate_genes_functional_categories_piechart.png"),
  plot = p_functional,
  width = 9,
  height = 6,
  dpi = 600
)

ggsave(
  filename = file.path(outdir, "all_candidate_genes_functional_categories_piechart.svg"),
  plot = p_functional,
  width = 9,
  height = 6
)

#========================================================
# 47. PIE CHART: OVERLAPPING VS NEARBY GENES FOR ALL CANDIDATES
#========================================================

candidate_location_categories <- all_candidate_genes_enriched %>%
  count(match_type) %>%
  mutate(
    label = paste0(match_type, " (", n, ")")
  )

cat("\n==============================\n")
cat("LOCATION CATEGORIES FOR ALL CANDIDATE GENES\n")
cat("==============================\n")
print(candidate_location_categories)

p_location <- ggplot(candidate_location_categories, aes(x = "", y = n, fill = match_type)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  theme_void() +
  labs(title = "Location of all candidate genes relative to candidate regions")

print(p_location)

ggsave(
  filename = file.path(outdir, "all_candidate_genes_overlap_vs_nearby_piechart.png"),
  plot = p_location,
  width = 7,
  height = 6,
  dpi = 600
)

ggsave(
  filename = file.path(outdir, "all_candidate_genes_overlap_vs_nearby_piechart.svg"),
  plot = p_location,
  width = 7,
  height = 6
)

#========================================================
# 48. BAR CHART: ANNOTATION CATEGORIES FOR ALL CANDIDATES
#========================================================

candidate_gene_categories <- all_candidate_genes_enriched %>%
  mutate(
    annotation_class = case_when(
      !is.na(`AL's functional class`) &
        grepl("secreted", `AL's functional class`, ignore.case = TRUE) ~ "Secreted",
      grepl("hypothetical|uncharacterized|unknown", preferred_annotation, ignore.case = TRUE) ~ "Uncharacterised",
      TRUE ~ "Annotated"
    )
  ) %>%
  count(annotation_class) %>%
  mutate(
    proportion = n / sum(n),
    label = paste0(annotation_class, " (", n, ")")
  )

p_annotation_bar <- ggplot(candidate_gene_categories, aes(x = annotation_class, y = n)) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "Annotation categories of all candidate genes",
    x = "Category",
    y = "Number of genes"
  )

print(p_annotation_bar)

ggsave(
  filename = file.path(outdir, "all_candidate_genes_annotation_categories_barchart.png"),
  plot = p_annotation_bar,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  filename = file.path(outdir, "all_candidate_genes_annotation_categories_barchart.svg"),
  plot = p_annotation_bar,
  width = 7,
  height = 5
)

#========================================================
# 49. LIST ANNOTATED GENES AMONG ALL CANDIDATES
#========================================================

annotated_genes <- all_candidate_genes_enriched %>%
  filter(
    !grepl("hypothetical|uncharacterized|unknown", preferred_annotation, ignore.case = TRUE) &
      !grepl("secreted", `AL's functional class`, ignore.case = TRUE)
  ) %>%
  select(
    gene_name,
    preferred_annotation,
    `AL's functional class`,
    best_max_fst,
    best_max_dxy
  ) %>%
  arrange(desc(best_max_fst))

cat("\n==============================\n")
cat("ANNOTATED GENES AMONG ALL CANDIDATES\n")
cat("==============================\n")
print(annotated_genes, n = Inf)

write_csv(
  annotated_genes,
  file.path(outdir, "all_candidate_genes_annotated_genes.csv")
)

#========================================================
# 50. LIST SECRETED GENES AMONG ALL CANDIDATES
#========================================================

secreted_genes <- all_candidate_genes_enriched %>%
  filter(
    !is.na(`AL's functional class`) &
      grepl("secreted", `AL's functional class`, ignore.case = TRUE)
  ) %>%
  select(
    gene_name,
    preferred_annotation,
    `AL's functional class`,
    best_max_fst,
    best_max_dxy
  ) %>%
  arrange(desc(best_max_fst))

cat("\n==============================\n")
cat("SECRETED GENES AMONG ALL CANDIDATES\n")
cat("==============================\n")
print(secreted_genes, n = Inf)

write_csv(
  secreted_genes,
  file.path(outdir, "all_candidate_genes_secreted_genes.csv")
)

#========================================================
# 51. LIST UNCHARACTERISED GENES AMONG ALL CANDIDATES
#========================================================

uncharacterised_genes <- all_candidate_genes_enriched %>%
  filter(
    is.na(preferred_annotation) |
      grepl("hypothetical|uncharacterized|unknown", preferred_annotation, ignore.case = TRUE)
  ) %>%
  select(
    gene_name,
    preferred_annotation,
    `Original gene description`,
    `AL's functional class`,
    best_max_fst,
    best_max_dxy
  ) %>%
  arrange(desc(best_max_fst))

cat("\n==============================\n")
cat("UNCHARACTERISED GENES AMONG ALL CANDIDATES\n")
cat("==============================\n")
print(uncharacterised_genes, n = Inf)

write_csv(
  uncharacterised_genes,
  file.path(outdir, "all_candidate_genes_uncharacterised_genes.csv")
)

#========================================================
# 52. QC FOR ALL CANDIDATE GENES
#========================================================

cat("\n==============================\n")
cat("QC: OVERLAPPING VS NEARBY")
cat("\n==============================\n")
print(table(genes_near_priority_10kb$overlap_type))

cat("\n==============================\n")
cat("QC: GENE POSITION CHECK")
cat("\n==============================\n")
genes_near_priority_10kb %>%
  select(gene_name, chromosome, gene_start, gene_end,
         region_start, region_end, distance_to_region, overlap_type) %>%
  arrange(chromosome, region_start, gene_start) %>%
  print(n = 20)

cat("\n==============================\n")
cat("QC: OVERLAP LOGIC CHECK")
cat("\n==============================\n")
genes_near_priority_10kb %>%
  mutate(check_overlap = gene_start <= region_end & gene_end >= region_start) %>%
  count(overlap_type, check_overlap) %>%
  print()

cat("\n==============================\n")
cat("QC: CHROMOSOME MAPPING")
cat("\n==============================\n")
genes2 %>%
  count(gene_chr, chromosome, sort = TRUE) %>%
  print(n = 30)

cat("\n==============================\n")
cat("QC: UNMAPPED CHROMOSOMES")
cat("\n==============================\n")
genes2 %>%
  filter(is.na(chromosome)) %>%
  count(gene_chr, sort = TRUE) %>%
  print(n = Inf)

cat("\n==============================\n")
cat("QC: TRANSCRIPT LINKS PER GENE")
cat("\n==============================\n")
candidate_tx_links %>%
  count(gene_name, sort = TRUE) %>%
  print(n = 20)

cat("\n==============================\n")
cat("QC: TRANSCRIPT MATCHING")
cat("\n==============================\n")
candidate_tx_annot %>%
  summarise(
    total_rows = n(),
    matched_rows = sum(!is.na(transcript_name)),
    unmatched_rows = sum(is.na(transcript_name)),
    matched_fraction = matched_rows / total_rows
  ) %>%
  print()

cat("\n==============================\n")
cat("SCRIPT COMPLETE: ALL CANDIDATE GENE ANALYSIS FINISHED")
```


## Human vs pig candidate-region and full candidate-gene set (n = 58)


```r

setwd("S:/Research/Lauren_PopGen/pi TajimasD/human_pig")
outdir <- file.path(getwd(), "candidate_gene_outputs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Use one consistent output directory name throughout
base_dir <- outdir

#========================================================
# PURPOSE OF THIS SCRIPT
#========================================================
# This script:
# 1. Reads in windowed FST and dXY results
# 2. Keeps only the focal comparison: Human vs Pig
# 3. Removes duplicated windows
# 4. Calculates bootstrap confidence intervals (CIs) for
#    the genome-wide mean FST and genome-wide mean dXY
# 5. Identifies top 1% and top 5% outlier windows
# 6. Tests whether FST and dXY outliers overlap exactly
# 7. Tests whether FST and dXY outliers cluster on the
#    same chromosomes
# 8. Tests whether FST and dXY outliers occur nearby
#    on the same chromosome even if they are not the
#    exact same window
# 9. Merges shared outlier windows into candidate regions
# 10. Intersects candidate regions with annotated genes
# 11. Uses tx61820 transcript links to recover additional
#     transcript-level annotation where possible
# 12. Builds a final cleaned candidate gene summary table
#     for figures and tables
#========================================================

library(tidyverse)
library(readxl)

set.seed(123)

#========================================================
# 1. READ INPUT FILES
#========================================================

fst <- read_tsv("human_vs_pig_fst.txt", show_col_types = FALSE)
dxy <- read_tsv("human_vs_pig_dxy.txt", show_col_types = FALSE)

#========================================================
# 2. INSPECT WHICH COMPARISONS ARE PRESENT
#========================================================

cat("\n==============================\n")
cat("UNIQUE FST COMPARISONS\n")
cat("==============================\n")
print(fst %>% distinct(pop1, pop2), n = Inf)

cat("\n==============================\n")
cat("UNIQUE dXY COMPARISONS\n")
cat("==============================\n")
print(dxy %>% distinct(pop1, pop2), n = Inf)

#========================================================
# 3. DEFINE THE VALUE COLUMNS
#========================================================

fst_col <- "avg_wc_fst"
dxy_col <- "avg_dxy"

#========================================================
# 4. KEEP ONLY THE FOCAL COMPARISON
#========================================================

fst <- fst %>%
  filter((pop1 == "Human" & pop2 == "Pig") | (pop1 == "Pig" & pop2 == "Human"))

dxy <- dxy %>%
  filter((pop1 == "Human" & pop2 == "Pig") | (pop1 == "Pig" & pop2 == "Human"))

#========================================================
# 5. REMOVE ANY DUPLICATED WINDOWS
#========================================================

fst <- fst %>%
  distinct(chromosome, window_pos_1, window_pos_2, .keep_all = TRUE)

dxy <- dxy %>%
  distinct(chromosome, window_pos_1, window_pos_2, .keep_all = TRUE)

#========================================================
# 6. BASIC SANITY CHECKS
#========================================================

cat("\n==============================\n")
cat("CHECKS AFTER FILTERING\n")
cat("==============================\n")
cat("FST rows:", nrow(fst), "\n")
cat("dXY rows:", nrow(dxy), "\n")

cat("\nUnique windows:\n")
cat("FST unique windows:", n_distinct(fst$chromosome, fst$window_pos_1, fst$window_pos_2), "\n")
cat("dXY unique windows:", n_distinct(dxy$chromosome, dxy$window_pos_1, dxy$window_pos_2), "\n")

#========================================================
# 7. HELPER FUNCTION: BOOTSTRAP CI FOR MEAN
#========================================================

boot_mean_ci <- function(x, nboot = 1000) {
  x <- x[is.finite(x)]
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    stop("No valid values found after filtering.")
  }
  
  boots <- replicate(nboot, {
    mean(sample(x, size = length(x), replace = TRUE))
  })
  
  tibble(
    n_windows = length(x),
    mean      = mean(x),
    median    = median(x),
    sd        = sd(x),
    min       = min(x),
    max       = max(x),
    ci_lower  = unname(quantile(boots, 0.025, na.rm = TRUE)),
    ci_upper  = unname(quantile(boots, 0.975, na.rm = TRUE))
  )
}

#========================================================
# 8. RUN BOOTSTRAP FOR GENOME-WIDE FST AND dXY
#========================================================

fst_summary <- boot_mean_ci(fst[[fst_col]], nboot = 1000) %>%
  mutate(statistic = "FST")

dxy_summary <- boot_mean_ci(dxy[[dxy_col]], nboot = 1000) %>%
  mutate(statistic = "dXY")

summary_table <- bind_rows(fst_summary, dxy_summary) %>%
  select(statistic, n_windows, mean, median, sd, min, max, ci_lower, ci_upper)

cat("\n==============================\n")
cat("BOOTSTRAP SUMMARY TABLE\n")
cat("==============================\n")
print(as.data.frame(summary_table))

ci_table <- summary_table %>%
  select(statistic, mean, ci_lower, ci_upper)

cat("\n==============================\n")
cat("BOOTSTRAP CI TABLE ONLY\n")
cat("==============================\n")
print(as.data.frame(ci_table))

write_csv(summary_table, file.path(outdir, "human_vs_pig_bootstrap_summary.csv"))
write_csv(ci_table, file.path(outdir, "human_vs_pig_bootstrap_ci_only.csv"))

#========================================================
# 9. PREPARE WINDOWED DATA FOR OUTLIER ANALYSIS
#========================================================

fst2 <- fst %>%
  filter(is.finite(.data[[fst_col]]), !is.na(.data[[fst_col]])) %>%
  mutate(midpoint = (window_pos_1 + window_pos_2) / 2)

dxy2 <- dxy %>%
  filter(is.finite(.data[[dxy_col]]), !is.na(.data[[dxy_col]])) %>%
  mutate(midpoint = (window_pos_1 + window_pos_2) / 2)

#========================================================
# 10. TOP 1% OUTLIERS
#========================================================

fst_cutoff_1 <- quantile(fst2[[fst_col]], 0.99, na.rm = TRUE)
dxy_cutoff_1 <- quantile(dxy2[[dxy_col]], 0.99, na.rm = TRUE)

cat("\n==============================\n")
cat("TOP 1% CUTOFFS\n")
cat("==============================\n")
cat("FST cutoff:", fst_cutoff_1, "\n")
cat("dXY cutoff:", dxy_cutoff_1, "\n")

fst_outliers_1 <- fst2 %>%
  filter(.data[[fst_col]] >= fst_cutoff_1) %>%
  arrange(desc(.data[[fst_col]])) %>%
  mutate(comparison = "Human_vs_Pig") %>%
  select(comparison, chromosome, window_pos_1, window_pos_2, midpoint, all_of(fst_col))

dxy_outliers_1 <- dxy2 %>%
  filter(.data[[dxy_col]] >= dxy_cutoff_1) %>%
  arrange(desc(.data[[dxy_col]])) %>%
  mutate(comparison = "Human_vs_Pig") %>%
  select(comparison, chromosome, window_pos_1, window_pos_2, midpoint, all_of(dxy_col))

#========================================================
# 11. EXACT OVERLAP OF TOP 1% OUTLIERS
#========================================================

shared_outliers_1 <- fst_outliers_1 %>%
  inner_join(
    dxy_outliers_1,
    by = c("comparison", "chromosome", "window_pos_1", "window_pos_2", "midpoint")
  ) %>%
  arrange(desc(avg_wc_fst), desc(avg_dxy))

cat("\n==============================\n")
cat("TOP 1% FST OUTLIERS\n")
cat("==============================\n")
print(fst_outliers_1, n = 20)

cat("\n==============================\n")
cat("TOP 1% dXY OUTLIERS\n")
cat("==============================\n")
print(dxy_outliers_1, n = 20)

cat("\n==============================\n")
cat("SHARED TOP 1% WINDOWS\n")
cat("==============================\n")
print(shared_outliers_1, n = 50)
cat("Number of shared top 1% windows:", nrow(shared_outliers_1), "\n")

write_csv(fst_outliers_1, file.path(outdir, "human_vs_pig_fst_top1pct_windows.csv"))
write_csv(dxy_outliers_1, file.path(outdir, "human_vs_pig_dxy_top1pct_windows.csv"))
write_csv(shared_outliers_1, file.path(outdir, "human_vs_pig_shared_top1pct_windows.csv"))

#========================================================
# 12. TOP 5% OUTLIERS
#========================================================

fst_cutoff_5 <- quantile(fst2[[fst_col]], 0.95, na.rm = TRUE)
dxy_cutoff_5 <- quantile(dxy2[[dxy_col]], 0.95, na.rm = TRUE)

cat("\n==============================\n")
cat("TOP 5% CUTOFFS\n")
cat("==============================\n")
cat("FST cutoff:", fst_cutoff_5, "\n")
cat("dXY cutoff:", dxy_cutoff_5, "\n")

fst_outliers_5 <- fst2 %>%
  filter(.data[[fst_col]] >= fst_cutoff_5) %>%
  arrange(desc(.data[[fst_col]])) %>%
  mutate(comparison = "Human_vs_Pig") %>%
  select(comparison, chromosome, window_pos_1, window_pos_2, midpoint, all_of(fst_col))

dxy_outliers_5 <- dxy2 %>%
  filter(.data[[dxy_col]] >= dxy_cutoff_5) %>%
  arrange(desc(.data[[dxy_col]])) %>%
  mutate(comparison = "Human_vs_Pig") %>%
  select(comparison, chromosome, window_pos_1, window_pos_2, midpoint, all_of(dxy_col))

#========================================================
# 13. EXACT OVERLAP OF TOP 5% OUTLIERS
#========================================================

shared_outliers_5 <- fst_outliers_5 %>%
  inner_join(
    dxy_outliers_5,
    by = c("comparison", "chromosome", "window_pos_1", "window_pos_2", "midpoint")
  ) %>%
  arrange(desc(avg_wc_fst), desc(avg_dxy))

cat("\n==============================\n")
cat("SHARED TOP 5% WINDOWS\n")
cat("==============================\n")
print(shared_outliers_5, n = 50)
cat("Number of shared top 5% windows:", nrow(shared_outliers_5), "\n")

write_csv(fst_outliers_5, file.path(outdir, "human_vs_pig_fst_top5pct_windows.csv"))
write_csv(dxy_outliers_5, file.path(outdir, "human_vs_pig_dxy_top5pct_windows.csv"))
write_csv(shared_outliers_5, file.path(outdir, "human_vs_pig_shared_top5pct_windows.csv"))

#========================================================
# 14. CHROMOSOME-LEVEL CLUSTERING: TOP 1%
#========================================================

fst_chr_counts_1 <- fst_outliers_1 %>%
  count(chromosome, sort = TRUE) %>%
  rename(fst_top1_count = n)

dxy_chr_counts_1 <- dxy_outliers_1 %>%
  count(chromosome, sort = TRUE) %>%
  rename(dxy_top1_count = n)

chr_overlap_1 <- full_join(fst_chr_counts_1, dxy_chr_counts_1, by = "chromosome") %>%
  mutate(
    fst_top1_count = replace_na(fst_top1_count, 0),
    dxy_top1_count = replace_na(dxy_top1_count, 0),
    shared_chromosome = fst_top1_count > 0 & dxy_top1_count > 0
  ) %>%
  arrange(desc(shared_chromosome), desc(fst_top1_count + dxy_top1_count))

cat("\n==============================\n")
cat("CHROMOSOME OVERLAP FOR TOP 1% OUTLIERS\n")
cat("==============================\n")
print(chr_overlap_1, n = Inf)

write_csv(fst_chr_counts_1, file.path(outdir, "human_vs_pig_fst_top1pct_counts_by_chromosome.csv"))
write_csv(dxy_chr_counts_1, file.path(outdir, "human_vs_pig_dxy_top1pct_counts_by_chromosome.csv"))
write_csv(chr_overlap_1, file.path(outdir, "human_vs_pig_top1pct_chromosome_overlap.csv"))

#========================================================
# 15. CHROMOSOME-LEVEL CLUSTERING: TOP 5%
#========================================================

fst_chr_counts_5 <- fst_outliers_5 %>%
  count(chromosome, sort = TRUE) %>%
  rename(fst_top5_count = n)

dxy_chr_counts_5 <- dxy_outliers_5 %>%
  count(chromosome, sort = TRUE) %>%
  rename(dxy_top5_count = n)

chr_overlap_5 <- full_join(fst_chr_counts_5, dxy_chr_counts_5, by = "chromosome") %>%
  mutate(
    fst_top5_count = replace_na(fst_top5_count, 0),
    dxy_top5_count = replace_na(dxy_top5_count, 0),
    shared_chromosome = fst_top5_count > 0 & dxy_top5_count > 0
  ) %>%
  arrange(desc(shared_chromosome), desc(fst_top5_count + dxy_top5_count))

cat("\n==============================\n")
cat("CHROMOSOME OVERLAP FOR TOP 5% OUTLIERS\n")
cat("==============================\n")
print(chr_overlap_5, n = Inf)

write_csv(fst_chr_counts_5, file.path(outdir, "human_vs_pig_fst_top5pct_counts_by_chromosome.csv"))
write_csv(dxy_chr_counts_5, file.path(outdir, "human_vs_pig_dxy_top5pct_counts_by_chromosome.csv"))
write_csv(chr_overlap_5, file.path(outdir, "human_vs_pig_top5pct_chromosome_overlap.csv"))

#========================================================
# 16. FUNCTION TO TEST NEARBY OVERLAPS
#========================================================

find_nearby_overlaps <- function(fst_df, dxy_df, max_dist = 10000) {
  inner_join(
    fst_df %>% rename(fst_midpoint = midpoint, fst_value = avg_wc_fst),
    dxy_df %>% rename(dxy_midpoint = midpoint, dxy_value = avg_dxy),
    by = c("comparison", "chromosome"),
    relationship = "many-to-many"
  ) %>%
    mutate(distance = abs(fst_midpoint - dxy_midpoint)) %>%
    filter(distance <= max_dist) %>%
    arrange(chromosome, distance, desc(fst_value), desc(dxy_value))
}

#========================================================
# 17. NEARBY OVERLAP: TOP 1%
#========================================================

nearby_overlap_1_10kb <- find_nearby_overlaps(fst_outliers_1, dxy_outliers_1, max_dist = 10000)
nearby_overlap_1_50kb <- find_nearby_overlaps(fst_outliers_1, dxy_outliers_1, max_dist = 50000)

cat("\n==============================\n")
cat("NEARBY OVERLAP FOR TOP 1% OUTLIERS (WITHIN 10 kb)\n")
cat("==============================\n")
print(nearby_overlap_1_10kb, n = 50)
cat("Number of nearby top 1% overlaps within 10 kb:", nrow(nearby_overlap_1_10kb), "\n")

cat("\n==============================\n")
cat("NEARBY OVERLAP FOR TOP 1% OUTLIERS (WITHIN 50 kb)\n")
cat("==============================\n")
print(nearby_overlap_1_50kb, n = 50)
cat("Number of nearby top 1% overlaps within 50 kb:", nrow(nearby_overlap_1_50kb), "\n")

write_csv(nearby_overlap_1_10kb, file.path(outdir, "human_vs_pig_nearby_overlap_top1pct_10kb.csv"))
write_csv(nearby_overlap_1_50kb, file.path(outdir, "human_vs_pig_nearby_overlap_top1pct_50kb.csv"))

#========================================================
# 18. NEARBY OVERLAP: TOP 5%
#========================================================

nearby_overlap_5_10kb <- find_nearby_overlaps(fst_outliers_5, dxy_outliers_5, max_dist = 10000)
nearby_overlap_5_50kb <- find_nearby_overlaps(fst_outliers_5, dxy_outliers_5, max_dist = 50000)

cat("\n==============================\n")
cat("NEARBY OVERLAP FOR TOP 5% OUTLIERS (WITHIN 10 kb)\n")
cat("==============================\n")
print(nearby_overlap_5_10kb, n = 50)
cat("Number of nearby top 5% overlaps within 10 kb:", nrow(nearby_overlap_5_10kb), "\n")

cat("\n==============================\n")
cat("NEARBY OVERLAP FOR TOP 5% OUTLIERS (WITHIN 50 kb)\n")
cat("==============================\n")
print(nearby_overlap_5_50kb, n = 50)
cat("Number of nearby top 5% overlaps within 50 kb:", nrow(nearby_overlap_5_50kb), "\n")

write_csv(nearby_overlap_5_10kb, file.path(outdir, "human_vs_pig_nearby_overlap_top5pct_10kb.csv"))
write_csv(nearby_overlap_5_50kb, file.path(outdir, "human_vs_pig_nearby_overlap_top5pct_50kb.csv"))

#========================================================
# 19. OPTIONAL: COUNT UNIQUE WINDOWS WITH NEARBY MATCHES
#========================================================

nearby_fst_top1_10kb <- nearby_overlap_1_10kb %>%
  distinct(comparison, chromosome, window_pos_1.x, window_pos_2.x, fst_midpoint, fst_value)

nearby_dxy_top1_10kb <- nearby_overlap_1_10kb %>%
  distinct(comparison, chromosome, window_pos_1.y, window_pos_2.y, dxy_midpoint, dxy_value)

cat("\n==============================\n")
cat("UNIQUE TOP 1% WINDOWS WITH A NEARBY MATCH WITHIN 10 kb\n")
cat("==============================\n")
cat("Unique FST windows with a nearby dXY peak:", nrow(nearby_fst_top1_10kb), "\n")
cat("Unique dXY windows with a nearby FST peak:", nrow(nearby_dxy_top1_10kb), "\n")

write_csv(nearby_fst_top1_10kb, file.path(outdir, "human_vs_pig_unique_fst_top1pct_with_nearby_dxy_10kb.csv"))
write_csv(nearby_dxy_top1_10kb, file.path(outdir, "human_vs_pig_unique_dxy_top1pct_with_nearby_fst_10kb.csv"))

#========================================================
# 20. RANK EXACT SHARED TOP 5% WINDOWS
#========================================================

shared_top5_ranked <- shared_outliers_5 %>%
  mutate(
    fst_rank = min_rank(desc(avg_wc_fst)),
    dxy_rank = min_rank(desc(avg_dxy)),
    combined_rank_score = fst_rank + dxy_rank
  ) %>%
  arrange(combined_rank_score, desc(avg_wc_fst), desc(avg_dxy))

cat("\n==============================\n")
cat("RANKED SHARED TOP 5% WINDOWS\n")
cat("==============================\n")
print(shared_top5_ranked, n = 50)

write_csv(shared_top5_ranked, file.path(outdir, "human_vs_pig_ranked_shared_top5pct_windows.csv"))

#========================================================
# 21. RANK NEARBY TOP 1% OVERLAPS WITHIN 10 kb
#========================================================

nearby_top1_ranked_10kb <- nearby_overlap_1_10kb %>%
  mutate(
    fst_rank = min_rank(desc(fst_value)),
    dxy_rank = min_rank(desc(dxy_value)),
    dist_rank = min_rank(distance),
    combined_rank_score = fst_rank + dxy_rank + dist_rank
  ) %>%
  arrange(combined_rank_score, distance, desc(fst_value), desc(dxy_value))

cat("\n==============================\n")
cat("RANKED NEARBY TOP 1% OVERLAPS (10 kb)\n")
cat("==============================\n")
print(nearby_top1_ranked_10kb, n = 50)

write_csv(nearby_top1_ranked_10kb, file.path(outdir, "human_vs_pig_ranked_nearby_top1pct_10kb.csv"))

#========================================================
# 22. FUNCTION TO MERGE NEARBY WINDOWS INTO REGIONS
#========================================================

merge_windows_to_regions <- function(df, gap_allowance = 10000) {
  df %>%
    arrange(chromosome, window_pos_1, window_pos_2) %>%
    group_by(chromosome) %>%
    mutate(
      prev_end = lag(window_pos_2),
      new_region = if_else(
        is.na(prev_end) | window_pos_1 > (prev_end + gap_allowance),
        1L, 0L
      ),
      region_id = cumsum(new_region)
    ) %>%
    group_by(chromosome, region_id) %>%
    summarise(
      region_start = min(window_pos_1),
      region_end = max(window_pos_2),
      n_windows = n(),
      mean_fst = mean(avg_wc_fst, na.rm = TRUE),
      max_fst = max(avg_wc_fst, na.rm = TRUE),
      mean_dxy = mean(avg_dxy, na.rm = TRUE),
      max_dxy = max(avg_dxy, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(region_width_bp = region_end - region_start + 1)
}

#========================================================
# 23. CANDIDATE REGIONS FROM SHARED TOP 5% WINDOWS
#========================================================

shared_top5_regions <- merge_windows_to_regions(shared_outliers_5, gap_allowance = 10000) %>%
  arrange(desc(max_fst), desc(max_dxy), desc(n_windows))

cat("\n==============================\n")
cat("CANDIDATE REGIONS FROM SHARED TOP 5% WINDOWS\n")
cat("==============================\n")
print(shared_top5_regions, n = 50)

write_csv(shared_top5_regions, file.path(outdir, "human_vs_pig_candidate_regions_shared_top5pct.csv"))

#========================================================
# 24. BUILD A SHORTLIST OF TOP CANDIDATE REGIONS
#========================================================

priority_candidates <- shared_top5_regions %>%
  mutate(
    priority_score = rank(-max_fst, ties.method = "min") +
      rank(-max_dxy, ties.method = "min") +
      rank(-n_windows, ties.method = "min")
  ) %>%
  arrange(priority_score, desc(max_fst), desc(max_dxy), desc(n_windows))

cat("\n==============================\n")
cat("PRIORITY CANDIDATE REGIONS\n")
cat("==============================\n")
print(priority_candidates, n = 20)

write_csv(priority_candidates, file.path(outdir, "human_vs_pig_priority_candidate_regions.csv"))

#========================================================
# 25. LOAD GENE ANNOTATION FILE
#========================================================

genes <- read_excel(file.path(base_dir, "Ascaris_genes.15714.annotation_expression.xlsx"))

cat("\n==============================\n")
cat("GENE FILE COLUMNS\n")
cat("==============================\n")
print(colnames(genes))

#========================================================
# 26. BUILD CHROMOSOME LOOKUP TABLE
#========================================================

chr_key <- tibble(
  gene_chr = c(
    "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
    "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
    "chrX1","chrX2","chrX3","chrX4","chrX5","chrMT"
  ),
  chromosome = c(
    "CM024168.1","CM024169.1","CM024170.1","CM024171.1","CM024172.1",
    "CM024173.1","CM024174.1","CM024175.1","CM024176.1","CM024177.1",
    "CM024178.1","CM024179.1","CM024180.1","CM024181.1","CM024182.1",
    "CM024183.1","CM024184.1","CM024185.1","CM024186.1",
    "CM024187.1","CM024188.1","CM024189.1","CM024190.1","CM024191.1",
    "CM024213.1"
  )
)

#========================================================
# 27. KEEP AND RENAME RELEVANT GENE COLUMNS
#========================================================

genes2 <- genes %>%
  select(
    chr,
    start,
    end,
    name,
    `Original gene description`,
    `Protein function group name`,
    `Curated annotation`,
    `AL's updated description`,
    `AL's comments`,
    `AL's functional class`,
    ASC10_names,
    tx61820_names,
    `ISO-seq_names`,
    expression,
    enriched_genes,
    average,
    max
  ) %>%
  rename(
    gene_chr   = chr,
    gene_start = start,
    gene_end   = end,
    gene_name  = name
  ) %>%
  left_join(chr_key, by = "gene_chr")

cat("\n==============================\n")
cat("CHROMOSOME NAMING CHECK\n")
cat("==============================\n")
cat("Example chromosomes in genes file:\n")
print(head(unique(genes2$gene_chr), 20))

cat("\nMapped chromosomes in genes file:\n")
print(head(unique(genes2$chromosome), 20))

cat("\nExample chromosomes in candidate regions:\n")
print(head(unique(priority_candidates$chromosome), 20))

cat("\nNumber of genes with missing chromosome mapping:\n")
print(sum(is.na(genes2$chromosome)))

unmapped_gene_chromosomes <- genes2 %>%
  filter(is.na(chromosome)) %>%
  distinct(gene_chr)

cat("\nUnmapped gene chromosome names (if any):\n")
print(unmapped_gene_chromosomes, n = Inf)

#========================================================
# 28. GENES OVERLAPPING PRIORITY CANDIDATE REGIONS
#========================================================

genes_overlapping_priority <- genes2 %>%
  filter(!is.na(chromosome)) %>%
  inner_join(priority_candidates, by = "chromosome", relationship = "many-to-many") %>%
  filter(
    gene_start <= region_end,
    gene_end >= region_start
  ) %>%
  mutate(
    overlap_type = "overlapping",
    overlap_start = pmax(gene_start, region_start),
    overlap_end   = pmin(gene_end, region_end),
    overlap_bp    = overlap_end - overlap_start + 1
  ) %>%
  arrange(priority_score, chromosome, region_start, gene_start)

cat("\n==============================\n")
cat("GENES OVERLAPPING PRIORITY CANDIDATE REGIONS\n")
cat("==============================\n")
print(genes_overlapping_priority, n = 50)

write_csv(
  genes_overlapping_priority,
  file.path(outdir, "human_vs_pig_genes_overlapping_priority_regions.csv")
)

#========================================================
# 29. GENES NEAR PRIORITY CANDIDATE REGIONS (WITHIN 10 kb)
#========================================================

genes_near_priority_10kb <- genes2 %>%
  filter(!is.na(chromosome)) %>%
  inner_join(priority_candidates, by = "chromosome", relationship = "many-to-many") %>%
  mutate(
    distance_to_region = case_when(
      gene_end < region_start ~ region_start - gene_end,
      gene_start > region_end ~ gene_start - region_end,
      TRUE ~ 0
    )
  ) %>%
  filter(distance_to_region <= 10000) %>%
  mutate(
    overlap_type = case_when(
      distance_to_region == 0 ~ "overlapping",
      TRUE ~ "nearby_10kb"
    )
  ) %>%
  arrange(priority_score, distance_to_region, chromosome, region_start, gene_start)

cat("\n==============================\n")
cat("GENES OVERLAPPING OR NEAR PRIORITY REGIONS (<=10 kb)\n")
cat("==============================\n")
print(genes_near_priority_10kb, n = 50)

write_csv(
  genes_near_priority_10kb,
  file.path(outdir, "human_vs_pig_genes_overlapping_or_near_priority_regions_10kb.csv")
)

#========================================================
# 30. CREATE ONE SUMMARY ROW PER GENE
#========================================================

priority_gene_summary <- genes_near_priority_10kb %>%
  group_by(
    chromosome,
    gene_chr,
    gene_start,
    gene_end,
    gene_name,
    `Original gene description`,
    `Protein function group name`,
    `Curated annotation`,
    `AL's updated description`,
    `AL's comments`,
    `AL's functional class`,
    ASC10_names,
    tx61820_names,
    `ISO-seq_names`,
    expression,
    enriched_genes,
    average,
    max
  ) %>%
  summarise(
    n_matching_regions = n(),
    best_priority_score = min(priority_score, na.rm = TRUE),
    best_max_fst = max(max_fst, na.rm = TRUE),
    best_max_dxy = max(max_dxy, na.rm = TRUE),
    total_windows_in_matched_regions = sum(n_windows, na.rm = TRUE),
    min_distance_to_region = min(distance_to_region, na.rm = TRUE),
    match_type = case_when(
      min(distance_to_region, na.rm = TRUE) == 0 ~ "overlapping",
      TRUE ~ "nearby_10kb"
    ),
    .groups = "drop"
  ) %>%
  arrange(best_priority_score, min_distance_to_region, desc(best_max_fst), desc(best_max_dxy))

cat("\n==============================\n")
cat("PRIORITISED CANDIDATE GENE SUMMARY\n")
cat("==============================\n")
print(priority_gene_summary, n = 58)

write_csv(
  priority_gene_summary,
  file.path(outdir, "human_vs_pig_priority_candidate_genes_summary.csv")
)

#========================================================
# 31. RETAIN THE FULL 58-GENE CANDIDATE SET
#========================================================
# priority_gene_summary contained 58 candidate genes.
# All 58 were retained; this was not a subset selected
# from a larger candidate-gene set.
#
# The legacy object and output filename containing "top58"
# are retained so that this documented script matches the
# analysis outputs already generated.

if (nrow(priority_gene_summary) != 58) {
  stop(
    "Expected 58 human-versus-pig candidate genes, but found ",
    nrow(priority_gene_summary),
    "."
  )
}

top_candidate_genes <- priority_gene_summary %>%
  arrange(
    best_priority_score,
    min_distance_to_region,
    desc(best_max_fst),
    desc(best_max_dxy)
  ) %>%
  slice_head(n = 58)

cat("\n==============================\n")
cat("FULL HUMAN-VERSUS-PIG CANDIDATE-GENE SET (n = 58)\n")
cat("==============================\n")
print(top_candidate_genes, n = 58)

write_csv(
  top_candidate_genes,
  file.path(outdir, "human_vs_pig_top58_candidate_genes.csv")
)

#========================================================
# 32. LOAD TRANSCRIPT ANNOTATION FILE
#========================================================

tx_annot_raw <- read_excel(file.path(base_dir, "Ascaris_transcripts_2017_version.xlsx"))

cat("\nTRANSCRIPT FILE COLUMNS:\n")
print(colnames(tx_annot_raw))

tx_annot <- tx_annot_raw %>%
  rename_with(~tolower(.)) %>%
  rename(
    transcript_name = transcript_name,
    tx_chromosome = chromosome,
    tx_start = start,
    tx_end = end,
    tx_annotation = annotation
  ) %>%
  mutate(
    transcript_name = trimws(transcript_name)
  )

#========================================================
# 33. EXPAND tx61820 LINKS
#========================================================

candidate_tx_links <- top_candidate_genes %>%
  select(gene_name, tx61820_names) %>%
  separate_rows(tx61820_names, sep = ",\\s*") %>%
  mutate(tx61820_names = trimws(tx61820_names)) %>%
  filter(!is.na(tx61820_names),
         tx61820_names != "",
         tx61820_names != "NA")

cat("\n==============================\n")
cat("EXPANDED TRANSCRIPT LINKS\n")
cat("==============================\n")
print(candidate_tx_links %>% count(gene_name), n = 20)

#========================================================
# 34. JOIN TO TRANSCRIPT ANNOTATION
#========================================================

candidate_tx_annot <- candidate_tx_links %>%
  left_join(
    tx_annot,
    by = c("tx61820_names" = "transcript_name"),
    keep = TRUE
  ) %>%
  mutate(
    tx_gene_id = sub("_t[0-9]+$", "", tx61820_names)
  )

cat("\n==============================\n")
cat("TRANSCRIPT MATCHING SUMMARY\n")
cat("==============================\n")
print(
  candidate_tx_annot %>%
    summarise(
      total_rows = n(),
      matched_rows = sum(!is.na(tx_annotation)),
      unmatched_rows = sum(is.na(tx_annotation)),
      matched_fraction = matched_rows / total_rows
    )
)

cat("\nExample joined rows:\n")
print(
  candidate_tx_annot %>%
    select(gene_name, tx61820_names, transcript_name, tx_annotation) %>%
    head(20)
)

write_csv(
  candidate_tx_annot,
  file.path(outdir, "human_vs_pig_candidate_gene_transcript_annotation_long.csv")
)

#========================================================
# 35. SUMMARISE TRANSCRIPT ANNOTATION PER GENE
#========================================================

transcript_annotation_summary <- candidate_tx_annot %>%
  group_by(gene_name) %>%
  summarise(
    n_tx61820_links = n_distinct(tx61820_names),
    n_matched_transcripts = n_distinct(transcript_name[!is.na(transcript_name)]),
    recovered_transcript_names = paste(
      unique(transcript_name[!is.na(transcript_name)]),
      collapse = "; "
    ),
    recovered_annotations = tx_annotation[!is.na(tx_annotation) & tx_annotation != ""] %>%
      unique() %>%
      .[!grepl("hypothetical|uncharacterized|unknown", ., ignore.case = TRUE)] %>%
      { if (length(.) == 0) "hypothetical protein" else paste(., collapse = "; ") },
    .groups = "drop"
  )

write_csv(
  transcript_annotation_summary,
  file.path(outdir, "human_vs_pig_candidate_gene_transcript_annotation_summary.csv")
)

#========================================================
# 36. MERGE BACK INTO THE FULL 58-GENE SET + ANNOTATION LOGIC
#========================================================

top_candidate_genes_enriched <- top_candidate_genes %>%
  left_join(transcript_annotation_summary, by = "gene_name") %>%
  mutate(
    preferred_annotation = case_when(
      !is.na(`Curated annotation`) &
        `Curated annotation` != "" &
        `Curated annotation` != "NA" ~ `Curated annotation`,
      !is.na(`AL's updated description`) &
        `AL's updated description` != "" &
        `AL's updated description` != "NA" ~ `AL's updated description`,
      !is.na(recovered_annotations) &
        recovered_annotations != "" &
        recovered_annotations != "NA" ~ recovered_annotations,
      !is.na(`Original gene description`) &
        `Original gene description` != "" &
        `Original gene description` != "NA" ~ `Original gene description`,
      TRUE ~ "Uncharacterised"
    )
  )

cat("\n==============================\n")
cat("ENRICHED FULL 58-GENE CANDIDATE SET\n")
cat("==============================\n")
print(top_candidate_genes_enriched, n = 58)

write_csv(
  top_candidate_genes_enriched,
  file.path(outdir, "human_vs_pig_top58_candidate_genes_enriched.csv")
)

#========================================================
# 37. FUNCTIONAL CLASS SUMMARY
#========================================================

top58_summary <- top_candidate_genes_enriched %>%
  mutate(
    annotation_class = case_when(
      is.na(preferred_annotation) ~ "Uncharacterised",
      grepl("hypothetical|uncharacterized|unknown", preferred_annotation, ignore.case = TRUE) ~ "Uncharacterised",
      grepl("secreted", `AL's functional class`, ignore.case = TRUE) ~ "Secreted",
      TRUE ~ "Annotated"
    )
  ) %>%
  count(annotation_class)

cat("\n==============================\n")
cat("ANNOTATION SUMMARY: FULL 58-GENE SET\n")
cat("==============================\n")
print(top58_summary)

cat("\nTotal genes counted: ", sum(top58_summary$n), "\n")

#========================================================
# 38. PLOT: CHROMOSOME DISTRIBUTION
#========================================================

p_chr_dist <- ggplot(top_candidate_genes_enriched %>% count(gene_chr),
                     aes(x = reorder(gene_chr, n), y = n)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Distribution of all 58 candidate genes across chromosomes",
    x = "Chromosome",
    y = "Number of genes"
  )

print(p_chr_dist)

ggsave(
  filename = file.path(outdir, "human_vs_pig_top58_chr_distribution.png"),
  plot = p_chr_dist,
  width = 8,
  height = 5
)

ggsave(
  filename = file.path(outdir, "human_vs_pig_top58_chr_distribution.svg"),
  plot = p_chr_dist,
  width = 8,
  height = 5
)

#========================================================
# 39. PLOT: FST vs dXY
#========================================================

p_fst_dxy <- ggplot(top_candidate_genes_enriched,
                    aes(x = best_max_fst, y = best_max_dxy)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Differentiation vs divergence (all 58 candidate genes)",
    x = "Max FST",
    y = "Max dXY"
  )

print(p_fst_dxy)

ggsave(
  filename = file.path(outdir, "human_vs_pig_top58_fst_vs_dxy.png"),
  plot = p_fst_dxy,
  width = 6,
  height = 5
)

ggsave(
  filename = file.path(outdir, "human_vs_pig_top58_fst_vs_dxy.svg"),
  plot = p_fst_dxy,
  width = 6,
  height = 5
)

#========================================================
# 40. PLOT: DISTANCE TO REGION (ALL UNIQUE CANDIDATE GENES)
#========================================================

distance_plot_df <- priority_gene_summary %>%
  mutate(
    distance_category = case_when(
      min_distance_to_region == 0 ~ "Overlapping",
      min_distance_to_region > 0 & min_distance_to_region <= 1000 ~ "0–1 kb",
      min_distance_to_region > 1000 & min_distance_to_region <= 5000 ~ "1–5 kb",
      min_distance_to_region > 5000 & min_distance_to_region <= 10000 ~ "5–10 kb",
      TRUE ~ NA_character_
    ),
    distance_category = factor(
      distance_category,
      levels = c("Overlapping", "0–1 kb", "1–5 kb", "5–10 kb")
    )
  ) %>%
  count(distance_category, .drop = FALSE)

p_distance <- ggplot(distance_plot_df, aes(x = distance_category, y = n)) +
  geom_col(fill = "grey30", width = 0.6) +
  geom_text(aes(label = n), vjust = -0.4, size = 4) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Distance of candidate genes to prioritised differentiated regions",
    x = "Distance to nearest region",
    y = "Number of genes"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

print(p_distance)

ggsave(
  filename = file.path(outdir, "human_vs_pig_candidate_gene_distance_to_region.png"),
  plot = p_distance,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  filename = file.path(outdir, "human_vs_pig_candidate_gene_distance_to_region.svg"),
  plot = p_distance,
  width = 7,
  height = 5
)

#========================================================
# 41. QC
#========================================================

table(genes_near_priority_10kb$overlap_type)

genes_near_priority_10kb %>%
  select(gene_name, chromosome, gene_start, gene_end,
         region_start, region_end, distance_to_region, overlap_type) %>%
  arrange(chromosome, region_start, gene_start) %>%
  print(n = 20)

genes_near_priority_10kb %>%
  mutate(check_overlap = gene_start <= region_end & gene_end >= region_start) %>%
  count(overlap_type, check_overlap)

genes2 %>%
  count(gene_chr, chromosome, sort = TRUE) %>%
  print(n = 30)

genes2 %>%
  filter(is.na(chromosome)) %>%
  count(gene_chr, sort = TRUE) %>%
  print(n = Inf)

candidate_tx_links %>%
  count(gene_name, sort = TRUE) %>%
  print(n = 20)

candidate_tx_annot %>%
  summarise(
    total_rows = n(),
    matched_rows = sum(!is.na(transcript_name)),
    unmatched_rows = sum(is.na(transcript_name)),
    matched_fraction = matched_rows / total_rows
  )

candidate_tx_annot %>%
  filter(gene_name %in% c("ascaris_gene_04441", "ascaris_gene_04443")) %>%
  select(gene_name, tx61820_names, tx_gene_id, transcript_name, tx_annotation) %>%
  print(n = 50)

top_candidate_genes_enriched %>%
  mutate(
    annotation_source = case_when(
      !is.na(`Curated annotation`) & `Curated annotation` != "" & `Curated annotation` != "NA" ~ "Curated annotation",
      !is.na(`AL's updated description`) & `AL's updated description` != "" & `AL's updated description` != "NA" ~ "AL updated description",
      !is.na(recovered_annotations) & recovered_annotations != "" & recovered_annotations != "NA" ~ "Recovered transcript annotation",
      TRUE ~ "Original gene description"
    )
  ) %>%
  count(annotation_source)

gene_to_check <- "ascaris_gene_04441"

priority_gene_summary %>%
  filter(gene_name == gene_to_check)

candidate_tx_links %>%
  filter(gene_name == gene_to_check)

candidate_tx_annot %>%
  filter(gene_name == gene_to_check) %>%
  select(gene_name, tx61820_names, tx_gene_id, transcript_name, tx_annotation)

top_candidate_genes_enriched %>%
  filter(gene_name == gene_to_check) %>%
  select(gene_name, preferred_annotation, recovered_annotations)

#========================================================
# 42. SUMMARY OF TRANSCRIPT-LEVEL ANNOTATION GAIN
#========================================================

tx_gain_summary <- top_candidate_genes_enriched %>%
  mutate(
    gained_tx_annotation = !is.na(recovered_annotations) &
      recovered_annotations != "" &
      !grepl("hypothetical|uncharacterized|unknown", recovered_annotations, ignore.case = TRUE)
  ) %>%
  count(gained_tx_annotation)

cat("\n==============================\n")
cat("TRANSCRIPT-LEVEL ANNOTATION GAIN\n")
cat("==============================\n")
print(tx_gain_summary)

cat("\nGenes with additional non-hypothetical transcript annotation: ",
    tx_gain_summary$n[tx_gain_summary$gained_tx_annotation], "\n")
cat("Genes without additional informative transcript annotation: ",
    tx_gain_summary$n[!tx_gain_summary$gained_tx_annotation], "\n")

#========================================================
# 43. ANNOTATION CATEGORIES FOR THE FULL 58-GENE SET
#========================================================

candidate_gene_categories <- top_candidate_genes_enriched %>%
  mutate(
    annotation_class = case_when(
      !is.na(`AL's functional class`) &
        grepl("secreted", `AL's functional class`, ignore.case = TRUE) ~ "Secreted",
      grepl("hypothetical|uncharacterized|unknown", preferred_annotation, ignore.case = TRUE) ~ "Uncharacterised",
      TRUE ~ "Annotated"
    )
  ) %>%
  count(annotation_class) %>%
  mutate(
    proportion = n / sum(n),
    label = paste0(annotation_class, " (", n, ")")
  )

cat("\n==============================\n")
cat("ANNOTATION CATEGORIES FOR THE FULL 58-GENE SET\n")
cat("==============================\n")
print(candidate_gene_categories)

#========================================================
# 44. LOAD GROUPED FUNCTIONAL METADATA
#========================================================

functional_annotation_grouped <- read_excel(
  file.path(base_dir, "functional_annotation_grouped_human_vs_pig_58.xlsx")
)

cat("\n==============================\n")
cat("FUNCTIONAL GROUP FILE COLUMNS\n")
cat("==============================\n")
print(colnames(functional_annotation_grouped))

functional_annotation_clean <- functional_annotation_grouped %>%
  rename_with(~tolower(.x)) %>%
  rename(
    gene_name = gene_name,
    Group = group
  ) %>%
  distinct(gene_name, .keep_all = TRUE)

cat("\n==============================\n")
cat("FUNCTIONAL GROUP COUNTS\n")
cat("==============================\n")
print(functional_annotation_clean %>% count(Group, sort = TRUE), n = Inf)

#========================================================
# 45. BUILD FINAL CLEAN SUMMARY TABLE FOR FIGURE/TABLE
#========================================================

candidate_gene_summary_table_for_figure <- top_candidate_genes_enriched %>%
  left_join(functional_annotation_clean, by = "gene_name") %>%
  transmute(
    `Gene Name` = gene_name,
    `Max fst` = round(best_max_fst, 4),
    `Max dxy` = round(best_max_dxy, 4),
    `Protein Annotation` = preferred_annotation,
    `Functional Class` = replace_na(`AL's functional class`, "Unknown"),
    Group
  ) %>%
  mutate(
    Group = replace_na(Group, "Unknown / poorly characterised"),
    `Protein Annotation` = case_when(
      `Gene Name` == "ascaris_gene_04439" ~ "mRNA-decapping enzyme 2",
      `Gene Name` == "ascaris_gene_06316" ~ "AATF protein",
      `Gene Name` == "ascaris_gene_06314" ~ "Pseudouridine-metabolizing protein",
      `Gene Name` == "ascaris_gene_03967" ~ "39S ribosomal protein",
      `Gene Name` == "ascaris_gene_04441" ~ "Uncharacterized protein C28H8.8",
      `Gene Name` == "ascaris_gene_06317" ~ "Domain-containing protein",
      `Gene Name` == "ascaris_gene_10841" ~ "Venom allergen 5",
      `Gene Name` == "ascaris_gene_03968" ~ "Ligand-gated ion channel 50",
      `Gene Name` == "ascaris_gene_03948" ~ "F-box WD repeat-containing protein 7",
      `Gene Name` == "ascaris_gene_03947" ~ "Ubiquinol-cytochrome-c reductase complex assembly factor 2",
      `Gene Name` == "ascaris_gene_06212" ~ "Serpentine receptor class X45",
      `Gene Name` == "ascaris_gene_09448" ~ "UDP-glucuronosyltransferase ugt-48",
      `Gene Name` == "ascaris_gene_09449" ~ "Acyl-CoA-binding domain-containing protein 5",
      `Gene Name` == "ascaris_gene_09447" ~ "Transmembrane protease serine 6",
      `Gene Name` == "ascaris_gene_01953" ~ "Bifunctional nitrilase/nitrile hydratase NIT4A",
      `Gene Name` == "ascaris_gene_01954" ~ "Chymotrypsin/elastase isoinhibitor 1",
      str_detect(str_to_lower(`Protein Annotation`), "hypothetical") ~ "Hypothetical protein",
      TRUE ~ `Protein Annotation`
    ),
    `Functional Class` = case_when(
      `Gene Name` %in% c("ascaris_gene_04439", "ascaris_gene_06316") ~ "Transcription Machinery",
      `Gene Name` == "ascaris_gene_06314" ~ "RNA metabolism / regulation",
      `Gene Name` == "ascaris_gene_03967" ~ "Protein Synthesis",
      `Gene Name` == "ascaris_gene_04441" ~ "Secreted / Unknown",
      `Gene Name` == "ascaris_gene_06317" ~ "Secreted",
      `Gene Name` == "ascaris_gene_10841" ~ "Secreted / Unknown",
      `Gene Name` %in% c("ascaris_gene_03947", "ascaris_gene_00147",
                         "ascaris_gene_05833", "ascaris_gene_09448",
                         "ascaris_gene_09449", "ascaris_gene_01953") ~ "Energy / intermediary metabolism",
      `Gene Name` %in% c("ascaris_gene_06212", "ascaris_gene_01628",
                         "ascaris_gene_05834") ~ "Signal Transduction / membrane",
      `Gene Name` == "ascaris_gene_00888" ~ "Protein Export / membrane trafficking",
      `Gene Name` == "ascaris_gene_01627" ~ "Extracellular matrix / membrane",
      `Gene Name` %in% c("ascaris_gene_08124", "ascaris_gene_08125") ~ "Surface-associated receptor",
      `Gene Name` == "ascaris_gene_01954" ~ "Protease inhibitor / secreted",
      `Gene Name` == "ascaris_gene_09447" ~ "Secreted / Protease",
      str_detect(str_to_lower(`Protein Annotation`), "hypothetical") ~ "Unknown",
      TRUE ~ `Functional Class`
    )
  ) %>%
  arrange(Group, desc(`Max fst`), desc(`Max dxy`))

cat("\n==============================\n")
cat("FINAL CANDIDATE GENE SUMMARY TABLE\n")
cat("==============================\n")
print(candidate_gene_summary_table_for_figure, n = 58)

write_csv(
  candidate_gene_summary_table_for_figure,
  file.path(base_dir, "candidate_gene_summary_table_for_figure_58.csv")
)

#========================================================
# 45B. BUILD GENOMIC PLOTTING TABLE (WITH COMBINED SIGNAL)
#========================================================

candidate_gene_plot_df <- top_candidate_genes_enriched %>%
  left_join(functional_annotation_clean, by = "gene_name") %>%
  mutate(
    Group = replace_na(Group, "Unknown / poorly characterised"),
    gene_midpoint = (gene_start + gene_end) / 2,
    combined_signal = best_max_fst * best_max_dxy
  ) %>%
  select(
    gene_chr,
    chromosome,
    gene_name,
    gene_start,
    gene_end,
    gene_midpoint,
    Group,
    best_max_fst,
    best_max_dxy,
    combined_signal,
    best_priority_score,
    min_distance_to_region,
    match_type
  )

# Use actual chromosome sizes
chr_lengths <- tibble(
  gene_chr = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
               "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
               "chr18","chr19"),
  chromosome = c("CM024168.1","CM024169.1","CM024170.1","CM024171.1","CM024172.1",
                 "CM024173.1","CM024174.1","CM024175.1","CM024176.1","CM024177.1",
                 "CM024178.1","CM024179.1","CM024180.1","CM024181.1","CM024182.1",
                 "CM024183.1","CM024184.1","CM024185.1","CM024186.1"),
  chr_length = c(18835767,14205397,14128337,13918991,13427824,13248486,13137961,
                 12190940,11532163,11047307,10268763,9833456,8681520,8487591,
                 8253480,8141563,7348195,7084669,7082418)
)

candidate_gene_plot_df <- candidate_gene_plot_df %>%
  left_join(chr_lengths %>% select(gene_chr, chr_length), by = "gene_chr") %>%
  mutate(
    relative_pos = gene_midpoint / chr_length
  )

chr_levels_present <- c(
  "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
  "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19"
)

candidate_gene_plot_df <- candidate_gene_plot_df %>%
  filter(gene_chr %in% chr_levels_present) %>%
  mutate(
    gene_chr = factor(gene_chr, levels = chr_levels_present)
  ) %>%
  arrange(gene_chr, relative_pos)

write_csv(
  candidate_gene_plot_df,
  file.path(base_dir, "candidate_gene_plot_df_58.csv")
)

custom_cols_chrplot <- c(
  "Unknown / poorly characterised" = "#fd91a0",
  "Host–parasite interaction / secreted / interface" = "#fdc170",
  "Metabolism, mitochondrial function, and redox homeostasis" = "#ffec95",
  "Membrane transport and surface-associated proteins" = "#a2daff",
  "Gene expression and genome regulation" = "#c7fa97",
  "Protein synthesis and turnover" = "#dfabff",
  "Signalling and sensory / membrane receptors" = "#8fd3c7",
  "Structural / cell division" = "#ad63ba"
)

p_chr_location <- ggplot(
  candidate_gene_plot_df,
  aes(
    x = relative_pos,
    y = gene_chr,
    fill = Group,
    size = combined_signal
  )
) +
  geom_point(
    shape = 21,
    colour = "black",
    stroke = 0.35,
    alpha = 1,
    position = position_jitter(width = 0, height = 0.12)
  ) +
  scale_fill_manual(values = custom_cols_chrplot, drop = FALSE) +
  scale_size_continuous(range = c(1, 3.2)) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0)
  ) +
  labs(
    title = "Chromosomal location of all 58 candidate genes",
    x = "Relative position along chromosome",
    y = "Chromosome",
    fill = "Functional category",
    size = "Combined FST–dXY signal"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p_chr_location)

ggsave(
  filename = file.path(outdir, "top58_chromosomal_location_combined_signal.png"),
  plot = p_chr_location,
  width = 11,
  height = 6.5,
  dpi = 600
)

ggsave(
  filename = file.path(outdir, "top58_chromosomal_location_combined_signal.svg"),
  plot = p_chr_location,
  width = 11,
  height = 6.5
)

#========================================================
# 46. PIE CHART: FUNCTIONAL CATEGORIES
#========================================================

group_counts <- candidate_gene_summary_table_for_figure %>%
  count(Group) %>%
  mutate(
    proportion = n / sum(n),
    percent = round(proportion * 100, 1),
    Group_label = if_else(
      Group == "Unknown / poorly characterised",
      "Poorly characterised",
      Group
    ),
    slice_label = paste0(percent, "%")
  ) %>%
  arrange(desc(n)) %>%
  mutate(Group_label = factor(Group_label, levels = Group_label))

custom_cols <- c(
  "Poorly characterised" = "#fd91a0",
  "Host–parasite interaction / secreted / interface" = "#fdc170",
  "Metabolism, mitochondrial function, and redox homeostasis" = "#ffec95",
  "Membrane transport and surface-associated proteins" = "#a2daff",
  "Gene expression and genome regulation" = "#c7fa97",
  "Protein synthesis and turnover" = "#dfabff",
  "Signalling and sensory / membrane receptors" = "#8fd3c7",
  "Structural / cell division" = "#ad63ba"
)

p_pie <- ggplot(group_counts, aes(x = "", y = n, fill = Group_label)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y", start = 0) +
  geom_text(
    aes(label = slice_label),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_fill_manual(values = custom_cols, drop = FALSE) +
  theme_void() +
  labs(
    title = "Functional annotation of all 58 candidate genes",
    fill = "Category"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9)
  )

print(p_pie)

ggsave(
  filename = file.path(outdir, "top58_functional_categories_piechart.png"),
  plot = p_pie,
  width = 9,
  height = 6
)

ggsave(
  filename = file.path(outdir, "top58_functional_categories_piechart.svg"),
  plot = p_pie,
  width = 9,
  height = 6
)

#========================================================
# 47. PIE CHART: OVERLAPPING VS NEARBY GENES
#========================================================

candidate_location_categories <- top_candidate_genes_enriched %>%
  count(match_type) %>%
  mutate(
    label = paste0(match_type, " (", n, ")")
  )

cat("\n==============================\n")
cat("LOCATION CATEGORIES FOR ALL 58 CANDIDATE GENES\n")
cat("==============================\n")
print(candidate_location_categories)

p_loc <- ggplot(candidate_location_categories, aes(x = "", y = n, fill = match_type)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  theme_void() +
  labs(title = "Location of all 58 candidate genes relative to candidate regions")

print(p_loc)

ggsave(
  filename = file.path(outdir, "top58_overlap_vs_nearby_piechart.png"),
  plot = p_loc,
  width = 7,
  height = 6
)

ggsave(
  filename = file.path(outdir, "top58_overlap_vs_nearby_piechart.svg"),
  plot = p_loc,
  width = 7,
  height = 6
)

#========================================================
# 48. BAR CHART: ANNOTATION CATEGORIES
#========================================================

p_annot_bar <- ggplot(candidate_gene_categories, aes(x = annotation_class, y = n)) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "Annotation categories of all 58 candidate genes",
    x = "Category",
    y = "Number of genes"
  )

print(p_annot_bar)

ggsave(
  filename = file.path(outdir, "top58_annotation_categories_barchart.png"),
  plot = p_annot_bar,
  width = 7,
  height = 5
)

ggsave(
  filename = file.path(outdir, "top58_annotation_categories_barchart.svg"),
  plot = p_annot_bar,
  width = 7,
  height = 5
)

#========================================================
# 49. LIST ANNOTATED GENES
#========================================================

annotated_genes <- top_candidate_genes_enriched %>%
  filter(
    !grepl("hypothetical|uncharacterized|unknown", preferred_annotation, ignore.case = TRUE) &
      !grepl("secreted", `AL's functional class`, ignore.case = TRUE)
  ) %>%
  select(
    gene_name,
    preferred_annotation,
    `AL's functional class`,
    best_max_fst,
    best_max_dxy
  ) %>%
  arrange(desc(best_max_fst))

cat("\n==============================\n")
cat("ANNOTATED GENES IN THE FULL 58-GENE SET\n")
cat("==============================\n")
print(annotated_genes, n = 58)

write_csv(
  annotated_genes,
  file.path(outdir, "top58_annotated_genes.csv")
)

#========================================================
# 50. LIST SECRETED / INTERFACE GENES
#========================================================

secreted_genes <- top_candidate_genes_enriched %>%
  filter(
    !is.na(`AL's functional class`) &
      grepl("secreted|extracellular|protease", `AL's functional class`, ignore.case = TRUE)
  ) %>%
  select(
    gene_name,
    preferred_annotation,
    `AL's functional class`,
    best_max_fst,
    best_max_dxy
  ) %>%
  arrange(desc(best_max_fst))

cat("\n==============================\n")
cat("SECRETED / INTERFACE GENES IN THE FULL 58-GENE SET\n")
cat("==============================\n")
print(secreted_genes, n = 58)

write_csv(
  secreted_genes,
  file.path(outdir, "top58_secreted_interface_genes.csv")
)

#========================================================
# 51. LIST UNCHARACTERISED GENES
#========================================================

uncharacterised_genes <- top_candidate_genes_enriched %>%
  filter(
    is.na(preferred_annotation) |
      grepl("hypothetical|uncharacterized|unknown", preferred_annotation, ignore.case = TRUE)
  ) %>%
  select(
    gene_name,
    preferred_annotation,
    `Original gene description`,
    `AL's functional class`,
    best_max_fst,
    best_max_dxy
  ) %>%
  arrange(desc(best_max_fst))

cat("\n==============================\n")
cat("UNCHARACTERISED GENES IN THE FULL 58-GENE SET\n")
cat("==============================\n")
print(uncharacterised_genes, n = 58)

write_csv(
  uncharacterised_genes,
  file.path(outdir, "top58_uncharacterised_genes.csv")
```


## Human vs pig gene-level π and Tajima’s D extraction

```r

#========================================================
# GENE-LEVEL PI AND TAJIMA'S D SUMMARY
# HUMAN VS PIG
#========================================================

library(tidyverse)

#========================================================
# 0. DEFINE DIRECTORIES
#========================================================

base_dir <- "S:/Research/Lauren_PopGen/pi TajimasD/human_pig"
div_dir  <- file.path(base_dir, "candidate_gene_outputs")
outdir   <- div_dir

dir.create(outdir, showWarnings = FALSE)

cat("\n==============================\n")
cat("DIRECTORY CHECK\n")
cat("==============================\n")
cat("Base dir:", base_dir, "\n")
cat("Diversity/output dir:", div_dir, "\n")

cat("\nFiles in base directory:\n")
print(list.files(base_dir))

cat("\nFiles in candidate_gene_outputs:\n")
print(list.files(div_dir))

#========================================================
# 1. LOAD INPUT FILES
#========================================================

# Full human-versus-pig candidate-gene set (58 genes); legacy filename retained
top_candidate_genes_enriched <- read_csv(
  file.path(div_dir, "human_vs_pig_top58_candidate_genes_enriched.csv"),
  show_col_types = FALSE
)

# Diversity files are in base_dir, not div_dir
pi_human <- read_tsv(
  file.path(base_dir, "human_pi.txt"),
  show_col_types = FALSE
)

pi_pig <- read_tsv(
  file.path(base_dir, "pig_pi.txt"),
  show_col_types = FALSE
)

tajd_human <- read_tsv(
  file.path(base_dir, "human_tajima_d.txt"),
  show_col_types = FALSE
)

tajd_pig <- read_tsv(
  file.path(base_dir, "pig_tajima_d.txt"),
  show_col_types = FALSE
)

#========================================================
# 2. QUICK SANITY CHECKS
#========================================================

cat("\n==============================\n")
cat("TOP CANDIDATE GENES COLUMNS\n")
cat("==============================\n")
print(colnames(top_candidate_genes_enriched))

cat("\n==============================\n")
cat("PI HUMAN COLUMNS\n")
cat("==============================\n")
print(colnames(pi_human))

cat("\n==============================\n")
cat("PI PIG COLUMNS\n")
cat("==============================\n")
print(colnames(pi_pig))

cat("\n==============================\n")
cat("TAJIMA'S D HUMAN COLUMNS\n")
cat("==============================\n")
print(colnames(tajd_human))

cat("\n==============================\n")
cat("TAJIMA'S D PIG COLUMNS\n")
cat("==============================\n")
print(colnames(tajd_pig))

#========================================================
# 3. STANDARDISE COLUMN NAMES
#========================================================

pi_human2 <- pi_human %>%
  rename(
    chromosome = chromosome,
    window_pos_1 = window_pos_1,
    window_pos_2 = window_pos_2,
    pi = avg_pi
  ) %>%
  mutate(population = "Human")

pi_pig2 <- pi_pig %>%
  rename(
    chromosome = chromosome,
    window_pos_1 = window_pos_1,
    window_pos_2 = window_pos_2,
    pi = avg_pi
  ) %>%
  mutate(population = "Pig")

tajd_human2 <- tajd_human %>%
  rename(
    chromosome = chromosome,
    window_pos_1 = window_pos_1,
    window_pos_2 = window_pos_2,
    tajimas_d = tajima_d
  ) %>%
  mutate(population = "Human")

tajd_pig2 <- tajd_pig %>%
  rename(
    chromosome = chromosome,
    window_pos_1 = window_pos_1,
    window_pos_2 = window_pos_2,
    tajimas_d = tajima_d
  ) %>%
  mutate(population = "Pig")

pi_all <- bind_rows(pi_human2, pi_pig2)
tajd_all <- bind_rows(tajd_human2, tajd_pig2)

cat("\n==============================\n")
cat("CHECK STANDARDISED PI DATA\n")
cat("==============================\n")
print(pi_all %>% select(chromosome, window_pos_1, window_pos_2, pi, population) %>% head())

cat("\n==============================\n")
cat("CHECK STANDARDISED TAJIMA'S D DATA\n")
cat("==============================\n")
print(tajd_all %>% select(chromosome, window_pos_1, window_pos_2, tajimas_d, population) %>% head())

#========================================================
# 4. KEEP ONLY GENE COORDINATES NEEDED FOR OVERLAP
#========================================================

candidate_genes_for_overlap <- top_candidate_genes_enriched %>%
  select(
    gene_name,
    gene_chr,
    chromosome,
    gene_start,
    gene_end
  ) %>%
  distinct()

cat("\n==============================\n")
cat("CANDIDATE GENES FOR OVERLAP\n")
cat("==============================\n")
print(candidate_genes_for_overlap, n = 10)

#========================================================
# 5. OVERLAP PI WINDOWS WITH CANDIDATE GENES
#========================================================

gene_pi_overlap <- candidate_genes_for_overlap %>%
  inner_join(pi_all, by = "chromosome") %>%
  filter(
    gene_start <= window_pos_2,
    gene_end >= window_pos_1
  )

cat("\n==============================\n")
cat("GENE-PI OVERLAPS\n")
cat("==============================\n")
print(gene_pi_overlap, n = 20)

#========================================================
# 6. OVERLAP TAJIMA'S D WINDOWS WITH CANDIDATE GENES
#========================================================

gene_tajd_overlap <- candidate_genes_for_overlap %>%
  inner_join(tajd_all, by = "chromosome") %>%
  filter(
    gene_start <= window_pos_2,
    gene_end >= window_pos_1
  )

cat("\n==============================\n")
cat("GENE-TAJIMA'S D OVERLAPS\n")
cat("==============================\n")
print(gene_tajd_overlap, n = 20)

#========================================================
# 7. SUMMARISE PI PER GENE
#========================================================

gene_pi_summary <- gene_pi_overlap %>%
  group_by(gene_name, population) %>%
  summarise(
    n_pi_windows = n(),
    mean_pi = mean(pi, na.rm = TRUE),
    median_pi = median(pi, na.rm = TRUE),
    min_pi = min(pi, na.rm = TRUE),
    max_pi = max(pi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = population,
    values_from = c(n_pi_windows, mean_pi, median_pi, min_pi, max_pi),
    names_sep = "_"
  )

cat("\n==============================\n")
cat("GENE-LEVEL PI SUMMARY\n")
cat("==============================\n")
print(gene_pi_summary, n = 20)

#========================================================
# 8. SUMMARISE TAJIMA'S D PER GENE
#========================================================

gene_tajd_summary <- gene_tajd_overlap %>%
  group_by(gene_name, population) %>%
  summarise(
    n_tajd_windows = n(),
    mean_tajd = mean(tajimas_d, na.rm = TRUE),
    median_tajd = median(tajimas_d, na.rm = TRUE),
    min_tajd = min(tajimas_d, na.rm = TRUE),
    max_tajd = max(tajimas_d, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = population,
    values_from = c(n_tajd_windows, mean_tajd, median_tajd, min_tajd, max_tajd),
    names_sep = "_"
  )

cat("\n==============================\n")
cat("GENE-LEVEL TAJIMA'S D SUMMARY\n")
cat("==============================\n")
print(gene_tajd_summary, n = 20)

#========================================================
# 9. MERGE DIVERSITY SUMMARIES BACK INTO THE FULL 58-GENE TABLE
#========================================================

top_candidate_genes_with_diversity <- top_candidate_genes_enriched %>%
  left_join(gene_pi_summary, by = "gene_name") %>%
  left_join(gene_tajd_summary, by = "gene_name") %>%
  mutate(
    delta_mean_pi = mean_pi_Pig - mean_pi_Human,
    delta_mean_tajd = mean_tajd_Pig - mean_tajd_Human
  )

cat("\n==============================\n")
cat("TOP CANDIDATE GENES WITH DIVERSITY SUMMARY\n")
cat("==============================\n")
print(top_candidate_genes_with_diversity, n = 20)

write_csv(
  top_candidate_genes_with_diversity,
  file.path(outdir, "human_vs_pig_top58_candidate_genes_enriched_with_diversity.csv")
)

write_csv(
  gene_pi_summary,
  file.path(outdir, "human_vs_pig_top58_gene_pi_summary.csv")
)

write_csv(
  gene_tajd_summary,
  file.path(outdir, "human_vs_pig_top58_gene_tajd_summary.csv")
)

write_csv(
  gene_pi_overlap,
  file.path(outdir, "human_vs_pig_top58_gene_pi_window_overlaps.csv")
)

write_csv(
  gene_tajd_overlap,
  file.path(outdir, "human_vs_pig_top58_gene_tajd_window_overlaps.csv")
)

top58_clean_diversity <- top_candidate_genes_with_diversity %>%
  select(
    gene_name,
    preferred_annotation,
    best_max_fst,
    best_max_dxy,
    mean_pi_Human,
    mean_pi_Pig,
    delta_mean_pi,
    mean_tajd_Human,
    mean_tajd_Pig,
    delta_mean_tajd
  ) %>%
  mutate(
    across(
      c(best_max_fst, best_max_dxy,
        mean_pi_Human, mean_pi_Pig, delta_mean_pi,
        mean_tajd_Human, mean_tajd_Pig, delta_mean_tajd),
      ~ round(.x, 4)
    )
  ) %>%
  arrange(desc(best_max_fst), desc(best_max_dxy))

print(top58_clean_diversity, n = 58)

write_csv(
  top58_clean_diversity,
  file.path(outdir, "human_vs_pig_top58_candidate_genes_clean_with_diversity.csv")
)
```


## Africa vs non-Africa gene-level π and Tajima’s D extraction


```r

#========================================================
# GENE-LEVEL PI AND TAJIMA'S D SUMMARY
# AFRICA VS NON-AFRICA
# FULL CANDIDATE SET
#========================================================

library(tidyverse)

#========================================================
# 0. DEFINE DIRECTORIES
#========================================================

base_dir <- "S:/Research/Lauren_PopGen/pi TajimasD"
div_dir  <- file.path(base_dir, "candidate_gene_outputs")
outdir   <- div_dir

dir.create(outdir, showWarnings = FALSE)

cat("\n==============================\n")
cat("DIRECTORY CHECK\n")
cat("==============================\n")
cat("Base dir:", base_dir, "\n")
cat("Diversity dir:", div_dir, "\n")

cat("\nFiles in base directory:\n")
print(list.files(base_dir))

cat("\nFiles in diversity/output directory:\n")
print(list.files(div_dir))

#========================================================
# 1. LOAD INPUT FILES
#========================================================

# FULL Africa vs non-Africa candidate genes
africa_candidate_genes <- read_csv(
  file.path(div_dir, "africa_vs_nonafrica_priority_candidate_genes_summary.csv"),
  show_col_types = FALSE
)

# Diversity files
pi_africa <- read_tsv(
  file.path(div_dir, "africa_pi.txt"),
  show_col_types = FALSE
)

pi_nonafrica <- read_tsv(
  file.path(div_dir, "nonafrica_pi.txt"),
  show_col_types = FALSE
)

tajd_africa <- read_tsv(
  file.path(div_dir, "africa_tajima_d.txt"),
  show_col_types = FALSE
)

tajd_nonafrica <- read_tsv(
  file.path(div_dir, "nonafrica_tajima_d.txt"),
  show_col_types = FALSE
)

#========================================================
# 2. QUICK SANITY CHECKS
#========================================================

cat("\n==============================\n")
cat("AFRICA FULL CANDIDATE GENE COLUMNS\n")
cat("==============================\n")
print(colnames(africa_candidate_genes))

cat("\n==============================\n")
cat("PI AFRICA COLUMNS\n")
cat("==============================\n")
print(colnames(pi_africa))

cat("\n==============================\n")
cat("PI NON-AFRICA COLUMNS\n")
cat("==============================\n")
print(colnames(pi_nonafrica))

cat("\n==============================\n")
cat("TAJIMA'S D AFRICA COLUMNS\n")
cat("==============================\n")
print(colnames(tajd_africa))

cat("\n==============================\n")
cat("TAJIMA'S D NON-AFRICA COLUMNS\n")
cat("==============================\n")
print(colnames(tajd_nonafrica))

#========================================================
# 3. STANDARDISE COLUMN NAMES
#========================================================

pi_africa2 <- pi_africa %>%
  rename(
    chromosome = chromosome,
    window_pos_1 = window_pos_1,
    window_pos_2 = window_pos_2,
    pi = avg_pi
  ) %>%
  mutate(population = "Africa")

pi_nonafrica2 <- pi_nonafrica %>%
  rename(
    chromosome = chromosome,
    window_pos_1 = window_pos_1,
    window_pos_2 = window_pos_2,
    pi = avg_pi
  ) %>%
  mutate(population = "nonAfrica")

tajd_africa2 <- tajd_africa %>%
  rename(
    chromosome = chromosome,
    window_pos_1 = window_pos_1,
    window_pos_2 = window_pos_2,
    tajimas_d = tajima_d
  ) %>%
  mutate(population = "Africa")

tajd_nonafrica2 <- tajd_nonafrica %>%
  rename(
    chromosome = chromosome,
    window_pos_1 = window_pos_1,
    window_pos_2 = window_pos_2,
    tajimas_d = tajima_d
  ) %>%
  mutate(population = "nonAfrica")

pi_all <- bind_rows(pi_africa2, pi_nonafrica2)
tajd_all <- bind_rows(tajd_africa2, tajd_nonafrica2)

cat("\n==============================\n")
cat("CHECK STANDARDISED PI DATA\n")
cat("==============================\n")
print(pi_all %>% select(chromosome, window_pos_1, window_pos_2, pi, population) %>% head())

cat("\n==============================\n")
cat("CHECK STANDARDISED TAJIMA'S D DATA\n")
cat("==============================\n")
print(tajd_all %>% select(chromosome, window_pos_1, window_pos_2, tajimas_d, population) %>% head())

#========================================================
# 4. KEEP ONLY GENE COORDINATES NEEDED FOR OVERLAP
#========================================================

candidate_genes_for_overlap <- africa_candidate_genes %>%
  select(
    gene_name,
    gene_chr,
    chromosome,
    gene_start,
    gene_end
  ) %>%
  distinct()

cat("\n==============================\n")
cat("FULL AFRICA CANDIDATE GENES FOR OVERLAP\n")
cat("==============================\n")
print(candidate_genes_for_overlap, n = 20)

cat("\nNumber of unique Africa candidate genes:\n")
print(nrow(candidate_genes_for_overlap))

#========================================================
# 5. OVERLAP PI WINDOWS WITH CANDIDATE GENES
#========================================================

gene_pi_overlap <- candidate_genes_for_overlap %>%
  inner_join(pi_all, by = "chromosome") %>%
  filter(
    gene_start <= window_pos_2,
    gene_end >= window_pos_1
  )

cat("\n==============================\n")
cat("GENE-PI OVERLAPS\n")
cat("==============================\n")
print(gene_pi_overlap, n = 20)

#========================================================
# 6. OVERLAP TAJIMA'S D WINDOWS WITH CANDIDATE GENES
#========================================================

gene_tajd_overlap <- candidate_genes_for_overlap %>%
  inner_join(tajd_all, by = "chromosome") %>%
  filter(
    gene_start <= window_pos_2,
    gene_end >= window_pos_1
  )

cat("\n==============================\n")
cat("GENE-TAJIMA'S D OVERLAPS\n")
cat("==============================\n")
print(gene_tajd_overlap, n = 20)

#========================================================
# 7. SUMMARISE PI PER GENE
#========================================================

gene_pi_summary <- gene_pi_overlap %>%
  group_by(gene_name, population) %>%
  summarise(
    n_pi_windows = n(),
    mean_pi = mean(pi, na.rm = TRUE),
    median_pi = median(pi, na.rm = TRUE),
    min_pi = min(pi, na.rm = TRUE),
    max_pi = max(pi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = population,
    values_from = c(n_pi_windows, mean_pi, median_pi, min_pi, max_pi),
    names_sep = "_"
  )

cat("\n==============================\n")
cat("GENE-LEVEL PI SUMMARY\n")
cat("==============================\n")
print(gene_pi_summary, n = 20)

#========================================================
# 8. SUMMARISE TAJIMA'S D PER GENE
#========================================================

gene_tajd_summary <- gene_tajd_overlap %>%
  group_by(gene_name, population) %>%
  summarise(
    n_tajd_windows = n(),
    mean_tajd = mean(tajimas_d, na.rm = TRUE),
    median_tajd = median(tajimas_d, na.rm = TRUE),
    min_tajd = min(tajimas_d, na.rm = TRUE),
    max_tajd = max(tajimas_d, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = population,
    values_from = c(n_tajd_windows, mean_tajd, median_tajd, min_tajd, max_tajd),
    names_sep = "_"
  )

cat("\n==============================\n")
cat("GENE-LEVEL TAJIMA'S D SUMMARY\n")
cat("==============================\n")
print(gene_tajd_summary, n = 20)

#========================================================
# 9. MERGE DIVERSITY SUMMARIES BACK INTO FULL AFRICA TABLE
#========================================================

africa_candidate_genes_with_diversity <- africa_candidate_genes %>%
  left_join(gene_pi_summary, by = "gene_name") %>%
  left_join(gene_tajd_summary, by = "gene_name") %>%
  mutate(
    delta_mean_pi = mean_pi_nonAfrica - mean_pi_Africa,
    delta_mean_tajd = mean_tajd_nonAfrica - mean_tajd_Africa
  )

cat("\n==============================\n")
cat("FULL AFRICA CANDIDATE GENES WITH DIVERSITY SUMMARY\n")
cat("==============================\n")
print(africa_candidate_genes_with_diversity, n = 20)

write_csv(
  africa_candidate_genes_with_diversity,
  file.path(outdir, "africa_vs_nonafrica_full_candidate_genes_with_diversity.csv")
)

write_csv(
  gene_pi_summary,
  file.path(outdir, "africa_vs_nonafrica_full_gene_pi_summary.csv")
)

write_csv(
  gene_tajd_summary,
  file.path(outdir, "africa_vs_nonafrica_full_gene_tajd_summary.csv")
)

write_csv(
  gene_pi_overlap,
  file.path(outdir, "africa_vs_nonafrica_full_gene_pi_window_overlaps.csv")
)

write_csv(
  gene_tajd_overlap,
  file.path(outdir, "africa_vs_nonafrica_full_gene_tajd_window_overlaps.csv")
)

#========================================================
# 10. CLEAN SUMMARY TABLE
#========================================================

africa_clean_diversity <- africa_candidate_genes_with_diversity %>%
  select(
    gene_name,
    `Original gene description`,
    `Curated annotation`,
    `AL's updated description`,
    best_max_fst,
    best_max_dxy,
    mean_pi_Africa,
    mean_pi_nonAfrica,
    delta_mean_pi,
    mean_tajd_Africa,
    mean_tajd_nonAfrica,
    delta_mean_tajd
  ) %>%
  mutate(
    preferred_annotation = case_when(
      !is.na(`Curated annotation`) & `Curated annotation` != "" & `Curated annotation` != "NA" ~ `Curated annotation`,
      !is.na(`AL's updated description`) & `AL's updated description` != "" & `AL's updated description` != "NA" ~ `AL's updated description`,
      !is.na(`Original gene description`) & `Original gene description` != "" & `Original gene description` != "NA" ~ `Original gene description`,
      TRUE ~ "Uncharacterised"
    )
  ) %>%
  select(
    gene_name,
    preferred_annotation,
    best_max_fst,
    best_max_dxy,
    mean_pi_Africa,
    mean_pi_nonAfrica,
    delta_mean_pi,
    mean_tajd_Africa,
    mean_tajd_nonAfrica,
    delta_mean_tajd
  ) %>%
  mutate(
    across(
      c(best_max_fst, best_max_dxy,
        mean_pi_Africa, mean_pi_nonAfrica, delta_mean_pi,
        mean_tajd_Africa, mean_tajd_nonAfrica, delta_mean_tajd),
      ~ round(.x, 4)
    )
  ) %>%
  arrange(desc(best_max_fst), desc(best_max_dxy))

print(africa_clean_diversity, n = 100)

write_csv(
  africa_clean_diversity,
  file.path(outdir, "africa_vs_nonafrica_full_candidate_genes_clean_with_diversity.csv")
)
```


## Candidate-gene FST vs dXY shared/unique comparison figure


```r

library(tidyverse)
library(patchwork)
library(gridExtra)
library(grid)

#========================================================
# 10. PLOT: CANDIDATE-GENE FST vs dXY FOR BOTH COMPARISONS
#========================================================

# IMPORTANT:
# Africa vs non-Africa files are in the parent pi TajimasD folder.
# Human vs Pig files are in the human_pig folder.

project_dir <- "S:/Research/Lauren_PopGen/pi TajimasD"
human_pig_base_dir <- file.path(project_dir, "human_pig")
outdir <- file.path(human_pig_base_dir, "candidate_gene_outputs")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#------------------------------
# 10A. LOAD FULL CANDIDATE GENE TABLES
#------------------------------

africa_candidates <- read_csv(
  file.path(
    project_dir,
    "candidate_gene_outputs",
    "africa_vs_nonafrica_priority_candidate_genes_summary.csv"
  ),
  show_col_types = FALSE
) %>%
  mutate(comparison = "Africa vs non-Africa")

human_pig_candidates <- read_csv(
  file.path(
    human_pig_base_dir,
    "candidate_gene_outputs",
    "human_vs_pig_top58_candidate_genes_enriched.csv"
  ),
  show_col_types = FALSE
) %>%
  mutate(comparison = "Human vs Pig")

cat("\n==============================\n")
cat("FILE CHECK\n")
cat("==============================\n")
cat("Africa rows:", nrow(africa_candidates), "\n")
cat("Human/Pig rows:", nrow(human_pig_candidates), "\n")

#------------------------------
# 10B. COMBINE AND CLASSIFY SHARED/UNIQUE GENES
#------------------------------

africa_gene_names <- africa_candidates %>%
  distinct(gene_name) %>%
  pull(gene_name)

human_pig_gene_names <- human_pig_candidates %>%
  distinct(gene_name) %>%
  pull(gene_name)

candidate_fst_dxy <- bind_rows(
  africa_candidates,
  human_pig_candidates
) %>%
  filter(
    is.finite(best_max_fst),
    is.finite(best_max_dxy),
    !is.na(best_max_fst),
    !is.na(best_max_dxy)
  ) %>%
  mutate(
    comparison = factor(
      comparison,
      levels = c("Africa vs non-Africa", "Human vs Pig")
    ),
    candidate_status = case_when(
      gene_name %in% africa_gene_names &
        gene_name %in% human_pig_gene_names ~ "Shared",
      comparison == "Africa vs non-Africa" ~ "Africa vs non-Africa only",
      comparison == "Human vs Pig" ~ "Human vs Pig only",
      TRUE ~ "Unclassified"
    ),
    candidate_status = factor(
      candidate_status,
      levels = c(
        "Africa vs non-Africa only",
        "Human vs Pig only",
        "Shared"
      )
    )
  )

cat("\n==============================\n")
cat("CANDIDATE COUNTS\n")
cat("==============================\n")
print(candidate_fst_dxy %>% count(comparison))
print(candidate_fst_dxy %>% count(candidate_status))

# Check stacked points: multiple genes with same FST/dXY
stacked_points_check <- candidate_fst_dxy %>%
  count(comparison, best_max_fst, best_max_dxy, candidate_status, sort = TRUE) %>%
  filter(n > 1)

cat("\n==============================\n")
cat("STACKED POINTS CHECK\n")
cat("==============================\n")
print(stacked_points_check, n = 50)

write_csv(
  stacked_points_check,
  file.path(outdir, "candidate_genes_stacked_fst_dxy_points_check.csv")
)

#------------------------------
# 10C. DESCRIPTIVE SPEARMAN RELATIONSHIP BETWEEN FST AND dXY
#------------------------------
# Each row in candidate_fst_dxy represents a candidate gene.
# Several genes can be associated with the same differentiated
# region and can therefore inherit identical max FST and dXY
# values. Those stacked points are valid for displaying the
# candidate-gene set, but they are not fully independent
# observations for an inferential significance test.
#
# Spearman's rho is therefore retained as a descriptive summary
# of the monotonic pattern shown in the gene-level scatter plot.
# P-values are not calculated or interpreted here. A rho close
# to zero indicates little evidence of a consistent monotonic
# relationship; it should not be described as proof that no
# relationship exists.

fst_dxy_correlations <- candidate_fst_dxy %>%
  group_by(comparison) %>%
  summarise(
    n_genes = n(),

    n_unique_fst_dxy_pairs = n_distinct(
      paste(best_max_fst, best_max_dxy, sep = "_")
    ),

    spearman_rho = cor(
      best_max_fst,
      best_max_dxy,
      method = "spearman",
      use = "complete.obs"
    ),

    .groups = "drop"
  )

cat("\n==============================\n")
cat("DESCRIPTIVE GENE-LEVEL FST-dXY SPEARMAN CORRELATIONS\n")
cat("==============================\n")
print(fst_dxy_correlations)

write_csv(
  fst_dxy_correlations,
  file.path(outdir, "candidate_genes_fst_dxy_correlations.csv")
)

#------------------------------
# 10D. PLOT: SHARED VS COMPARISON-SPECIFIC CANDIDATE GENES
#      all points as circles
#------------------------------

custom_cols_fst_dxy <- c(
  "Africa vs non-Africa only" = "#ee4266",
  "Human vs Pig only" = "#ffd23f",
  "Shared" = "#540d6e"
)

p_fst_dxy_candidates <- ggplot(
  candidate_fst_dxy,
  aes(
    x = best_max_fst,
    y = best_max_dxy,
    fill = candidate_status
  )
) +
  geom_jitter(
    shape = 21,             
    width = 0.004,
    height = 0.0004,
    size = 3,
    alpha = 0.8,
    colour = "grey60",
    stroke = 0.25
  ) +
  scale_fill_manual(values = custom_cols_fst_dxy, drop = FALSE) +
  theme_minimal(base_size = 12) +
  labs(
    title = expression(paste(F[ST], " vs ", d[XY], " in candidate genes")),
    x = expression(paste("Max ", F[ST])),
    y = expression(paste("Max ", d[XY])),
    fill = "Candidate status"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p_fst_dxy_candidates)

ggsave(
  filename = file.path(outdir, "candidate_genes_fst_vs_dxy_shared_unique.png"),
  plot = p_fst_dxy_candidates,
  width = 9,
  height = 7,
  dpi = 600
)

ggsave(
  filename = file.path(outdir, "candidate_genes_fst_vs_dxy_shared_unique.svg"),
  plot = p_fst_dxy_candidates,
  width = 7,
  height = 5
)

#------------------------------
# 10E. COMBINE SCATTER PLOT WITH CORRELATION TABLE
#------------------------------

# Make a clean descriptive correlation table for plotting.
# No inferential p-values are shown because multiple genes can
# share the same candidate-region FST and dXY signal.
cor_table_plot <- fst_dxy_correlations %>%
  mutate(
    comparison = as.character(comparison),
    spearman_rho = round(spearman_rho, 3)
  ) %>%
  rename(
    `Comparison` = comparison,
    `n genes` = n_genes,
    `Unique FST-dXY pairs` = n_unique_fst_dxy_pairs,
    `Spearman rho` = spearman_rho
  )

# Convert table to a grob
cor_table_grob <- tableGrob(
  cor_table_plot,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 10,
    core = list(
      fg_params = list(hjust = 0.5, x = 0.5)
    ),
    colhead = list(
      fg_params = list(fontface = "bold", hjust = 0.5, x = 0.5)
    )
  )
)

# Add a title above the table
cor_table_title <- textGrob(
  expression(paste("Descriptive Spearman relationship between ", F[ST], " and ", d[XY])),
  gp = gpar(fontsize = 11, fontface = "bold"),
  x = 0,
  hjust = 0
)

cor_table_with_title <- arrangeGrob(
  cor_table_title,
  cor_table_grob,
  ncol = 1,
  heights = c(0.2, 1)
)

# Combine scatter plot and table
combined_fst_dxy_figure <- p_fst_dxy_candidates /
  wrap_elements(cor_table_with_title) +
  plot_layout(heights = c(3.2, 1))

print(combined_fst_dxy_figure)

# Save combined figure
ggsave(
  filename = file.path(outdir, "candidate_genes_fst_vs_dxy_with_correlation_table.png"),
  plot = combined_fst_dxy_figure,
  width = 8,
  height = 6.5,
  dpi = 600
)

ggsave(
  filename = file.path(outdir, "candidate_genes_fst_vs_dxy_with_correlation_table.svg"),
  plot = combined_fst_dxy_figure,
  width = 8,
  height = 6.5
)
```


## Genome-wide comparison figure

```r
library(tidyverse)
library(patchwork)
library(svglite)

setwd("S:/Research/Lauren_PopGen/pi TajimasD")

#========================================================
# GENOME-WIDE COMPARISON FIGURE
# Shared pi + Africa vs non-Africa + Human vs Pig
#========================================================

#---------------------------
# read files
#---------------------------
pi <- read_tsv("all_146_pi.txt", show_col_types = FALSE)

fst_africa <- read_tsv("africa_vs_nonafrica_fst.txt", show_col_types = FALSE)
dxy_africa <- read_tsv("africa_vs_nonafrica_dxy.txt", show_col_types = FALSE)

fst_hp <- read_tsv("human_pig/human_vs_pig_fst.txt", show_col_types = FALSE)
dxy_hp <- read_tsv("human_pig/human_vs_pig_dxy.txt", show_col_types = FALSE)

#---------------------------
# chromosome lookup and true assembly lengths
#---------------------------
# True chromosome lengths are used only to place chromosomes
# end-to-end on the cumulative x-axis. The plotted y-values remain
# the observed 10 kb window estimates, coloured by chromosome.

chr_lookup <- tibble(
  chromosome = c(
    "CM024168.1", "CM024169.1", "CM024170.1", "CM024171.1",
    "CM024172.1", "CM024173.1", "CM024174.1", "CM024175.1",
    "CM024176.1", "CM024177.1", "CM024178.1", "CM024179.1",
    "CM024180.1", "CM024181.1", "CM024182.1", "CM024183.1",
    "CM024184.1", "CM024185.1", "CM024186.1"
  ),
  chr_label = paste0("Chr", 1:19),
  chr_length_bp = c(
    18835767, 14205397, 14128337, 13918991, 13427824,
    13248486, 13137961, 12190940, 11532163, 11047307,
    10268763, 9833456, 8681520, 8487591, 8253480,
    8141563, 7348195, 7084669, 7082418
  )
)

#========================================================
# 1. PREPARE SHARED PI PANEL
#========================================================

pi_plot <- pi %>%
  left_join(chr_lookup, by = "chromosome") %>%
  mutate(midpoint = (window_pos_1 + window_pos_2) / 2) %>%
  filter(!is.na(avg_pi)) %>%
  distinct(chromosome, window_pos_1, window_pos_2, .keep_all = TRUE)

chrom_df <- chr_lookup %>%
  mutate(
    tot = cumsum(chr_length_bp) - chr_length_bp
  ) %>%
  select(
    chromosome,
    tot
  )

pi_plot2 <- pi_plot %>%
  left_join(chrom_df, by = "chromosome") %>%
  mutate(bp_cum = midpoint + tot)

axis_df <- chr_lookup %>%
  left_join(chrom_df, by = "chromosome") %>%
  transmute(
    chr_label,
    center = tot + chr_length_bp / 2
  )

mean_pi <- mean(pi_plot2$avg_pi, na.rm = TRUE)

#========================================================
# 2. HELPER FUNCTION FOR FST / DXY PANELS
#========================================================

prepare_genomewide_plot <- function(df, value_col) {
  df %>%
    left_join(chr_lookup, by = "chromosome") %>%
    mutate(midpoint = (window_pos_1 + window_pos_2) / 2) %>%
    filter(!is.na(.data[[value_col]])) %>%
    distinct(chromosome, window_pos_1, window_pos_2, .keep_all = TRUE) %>%
    left_join(chrom_df, by = "chromosome") %>%
    mutate(bp_cum = midpoint + tot)
}

#========================================================
# 3. PREPARE AFRICA VS NON-AFRICA
#========================================================

fst_africa2 <- prepare_genomewide_plot(fst_africa, "avg_wc_fst")
dxy_africa2 <- prepare_genomewide_plot(dxy_africa, "avg_dxy")

mean_fst_africa <- mean(fst_africa2$avg_wc_fst, na.rm = TRUE)
mean_dxy_africa <- mean(dxy_africa2$avg_dxy, na.rm = TRUE)

#========================================================
# 4. PREPARE HUMAN VS PIG
#========================================================

fst_hp2 <- fst_hp %>%
  filter((pop1 == "Human" & pop2 == "Pig") | (pop1 == "Pig" & pop2 == "Human")) %>%
  prepare_genomewide_plot("avg_wc_fst")

dxy_hp2 <- dxy_hp %>%
  filter((pop1 == "Human" & pop2 == "Pig") | (pop1 == "Pig" & pop2 == "Human")) %>%
  prepare_genomewide_plot("avg_dxy")

mean_fst_hp <- mean(fst_hp2$avg_wc_fst, na.rm = TRUE)
mean_dxy_hp <- mean(dxy_hp2$avg_dxy, na.rm = TRUE)

#========================================================
# 5. BUILD PANELS
#========================================================

# A. pi
p_pi <- ggplot(pi_plot2, aes(x = bp_cum / 1e6, y = avg_pi, colour = chromosome)) +
  geom_point(size = 0.25, alpha = 0.45) +
  geom_hline(
    yintercept = mean_pi,
    linetype = "dashed",
    colour = "black",
    linewidth = 0.8
  ) +
  scale_colour_hue(c = 80) +
  scale_x_continuous(
    breaks = axis_df$center / 1e6,
    labels = axis_df$chr_label,
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = expression(pi),
    title = "A. Genome-wide nucleotide diversity"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold", size = 11)
  )

# B. Africa vs non-Africa FST
p_fst_africa <- ggplot(fst_africa2, aes(x = bp_cum / 1e6, y = avg_wc_fst, colour = chromosome)) +
  geom_point(size = 0.25, alpha = 0.45) +
  geom_hline(
    yintercept = mean_fst_africa,
    linetype = "dashed",
    colour = "grey20",
    linewidth = 0.8
  ) +
  scale_colour_hue(c = 80) +
  scale_x_continuous(
    breaks = axis_df$center / 1e6,
    labels = axis_df$chr_label,
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(-0.3, 1.0)) +
  labs(
    x = NULL,
    y = expression(F[ST]),
    title = "B. Genome-wide genetic differentiation (Africa vs non-Africa)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold", size = 11)
  )

# C. Africa vs non-Africa dXY
p_dxy_africa <- ggplot(dxy_africa2, aes(x = bp_cum / 1e6, y = avg_dxy, colour = chromosome)) +
  geom_point(size = 0.25, alpha = 0.45) +
  geom_hline(
    yintercept = mean_dxy_africa,
    linetype = "dashed",
    colour = "grey20",
    linewidth = 0.8
  ) +
  scale_colour_hue(c = 80) +
  scale_x_continuous(
    breaks = axis_df$center / 1e6,
    labels = axis_df$chr_label,
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "dXY",
    title = "C. Genome-wide absolute divergence (Africa vs non-Africa)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold", size = 11)
  )

# D. Human vs Pig FST
p_fst_hp <- ggplot(fst_hp2, aes(x = bp_cum / 1e6, y = avg_wc_fst, colour = chromosome)) +
  geom_point(size = 0.25, alpha = 0.45) +
  geom_hline(
    yintercept = mean_fst_hp,
    linetype = "dashed",
    colour = "grey20",
    linewidth = 0.8
  ) +
  scale_colour_hue(c = 80) +
  scale_x_continuous(
    breaks = axis_df$center / 1e6,
    labels = axis_df$chr_label,
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(-0.3, 1.0)) +
  labs(
    x = NULL,
    y = expression(F[ST]),
    title = "D. Genome-wide genetic differentiation (Human vs Pig)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold", size = 11)
  )

# E. Human vs Pig dXY
p_dxy_hp <- ggplot(dxy_hp2, aes(x = bp_cum / 1e6, y = avg_dxy, colour = chromosome)) +
  geom_point(size = 0.25, alpha = 0.45) +
  geom_hline(
    yintercept = mean_dxy_hp,
    linetype = "dashed",
    colour = "grey20",
    linewidth = 0.8
  ) +
  scale_colour_hue(c = 80) +
  scale_x_continuous(
    breaks = axis_df$center / 1e6,
    labels = axis_df$chr_label,
    expand = c(0, 0)
  ) +
  labs(
    x = "Chromosome",
    y = "dXY",
    title = "E. Genome-wide absolute divergence (Human vs Pig)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 11)
  )

#========================================================
# 6. COMBINE FIGURE
#========================================================

comparison_fig <- p_pi / p_fst_africa / p_dxy_africa / p_fst_hp / p_dxy_hp +
  plot_layout(heights = c(1, 1, 1, 1, 1)) +
  plot_annotation(
    title = "Genome-wide diversity and differentiation across comparison schemes",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

print(comparison_fig)

#========================================================
# 7. SAVE OUTPUTS
#========================================================

ggsave(
  filename = "human_pig/Genome_wide_pi_FST_dXY_comparison_5panel.png",
  plot = comparison_fig,
  width = 12,
  height = 17,
  dpi = 600
)

ggsave(
  filename = "human_pig/Genome_wide_pi_FST_dXY_comparison_5panel.pdf",
  plot = comparison_fig,
  width = 12,
  height = 17
)

ggsave(
  filename = "human_pig/Genome_wide_pi_FST_dXY_comparison_5panel.svg",
  plot = comparison_fig,
  width = 12,
  height = 17,
  device = svglite
)

#========================================================
# 8. QC
#========================================================

cat("\n==============================\n")
cat("ROW COUNTS\n")
cat("==============================\n")
cat("pi rows:", nrow(pi_plot2), "\n")
cat("Africa FST rows:", nrow(fst_africa2), "\n")
cat("Africa dXY rows:", nrow(dxy_africa2), "\n")
cat("HumanPig FST rows:", nrow(fst_hp2), "\n")
cat("HumanPig dXY rows:", nrow(dxy_hp2), "\n")

cat("\n==============================\n")
cat("MEAN VALUES\n")
cat("==============================\n")
cat("Mean pi:", mean_pi, "\n")
cat("Mean Africa FST:", mean_fst_africa, "\n")
cat("Mean Africa dXY:", mean_dxy_africa, "\n")
cat("Mean HumanPig FST:", mean_fst_hp, "\n")
cat("Mean HumanPig dXY:", mean_dxy_hp, "\n")
```


## Supplementary genome-wide diversity/differentiation figure

```r
#========================================================
# SUPPLEMENTARY FIGURE:
# GENOME-WIDE DIVERSITY AND DIFFERENTIATION SUMMARIES
# INCLUDING HUMAN VS PIG FST AND dXY
#========================================================

library(tidyverse)
library(patchwork)

#========================================================
# 0. DEFINE DIRECTORIES
#========================================================

project_dir <- "S:/Research/Lauren_PopGen/pi TajimasD"
text_dir    <- file.path(project_dir, "text files")
outdir      <- file.path(project_dir, "candidate_gene_outputs")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("\n==============================\n")
cat("DIRECTORY CHECK\n")
cat("==============================\n")
cat("Project dir:", project_dir, "\n")
cat("Text input dir:", text_dir, "\n")
cat("Output dir:", outdir, "\n")

cat("\nFiles in text files directory:\n")
print(list.files(text_dir))

#========================================================
# 1. LOAD INPUT FILES FROM TEXT FILES FOLDER
#========================================================

# Nucleotide diversity
pi_africa <- read_tsv(file.path(text_dir, "africa_pi.txt"), show_col_types = FALSE)
pi_asia   <- read_tsv(file.path(text_dir, "asia_pi.txt"), show_col_types = FALSE)
pi_europe <- read_tsv(file.path(text_dir, "europe_pi.txt"), show_col_types = FALSE)

# Tajima's D
tajd_africa <- read_tsv(file.path(text_dir, "africa_tajima_d.txt"), show_col_types = FALSE)
tajd_asia   <- read_tsv(file.path(text_dir, "asia_tajima_d.txt"), show_col_types = FALSE)
tajd_europe <- read_tsv(file.path(text_dir, "europe_tajima_d.txt"), show_col_types = FALSE)

# Africa vs non-Africa FST and dXY
fst_af <- read_tsv(file.path(text_dir, "africa_vs_nonafrica_fst.txt"), show_col_types = FALSE)
dxy_af <- read_tsv(file.path(text_dir, "africa_vs_nonafrica_dxy.txt"), show_col_types = FALSE)

# Human vs Pig FST and dXY
# These input files contain only the Human-versus-Pig comparison,
# so no additional pop1/pop2 filtering is required here.
fst_hp <- read_tsv(file.path(text_dir, "human_vs_pig_fst.txt"), show_col_types = FALSE)
dxy_hp <- read_tsv(file.path(text_dir, "human_vs_pig_dxy.txt"), show_col_types = FALSE)

#========================================================
# 2. QUICK COLUMN CHECKS
#========================================================

cat("\n==============================\n")
cat("COLUMN CHECKS\n")
cat("==============================\n")

cat("\nAfrica pi columns:\n")
print(colnames(pi_africa))

cat("\nAfrica Tajima's D columns:\n")
print(colnames(tajd_africa))

cat("\nAfrica vs non-Africa FST columns:\n")
print(colnames(fst_af))

cat("\nAfrica vs non-Africa dXY columns:\n")
print(colnames(dxy_af))

cat("\nHuman vs Pig FST columns:\n")
print(colnames(fst_hp))

cat("\nHuman vs Pig dXY columns:\n")
print(colnames(dxy_hp))

#========================================================
# 3. STANDARDISE AND FILTER DATA
#========================================================

pi_all <- bind_rows(
  pi_africa %>% mutate(group = "Africa"),
  pi_asia   %>% mutate(group = "Asia"),
  pi_europe %>% mutate(group = "Europe")
) %>%
  filter(
    is.finite(avg_pi),
    !is.na(avg_pi)
  )

tajd_all <- bind_rows(
  tajd_africa %>% mutate(group = "Africa"),
  tajd_asia   %>% mutate(group = "Asia"),
  tajd_europe %>% mutate(group = "Europe")
) %>%
  filter(
    is.finite(tajima_d),
    !is.na(tajima_d)
  )

fst_af2 <- fst_af %>%
  filter(
    is.finite(avg_wc_fst),
    !is.na(avg_wc_fst)
  )

dxy_af2 <- dxy_af %>%
  filter(
    is.finite(avg_dxy),
    !is.na(avg_dxy)
  )

fst_hp2 <- fst_hp %>%
  filter(
    is.finite(avg_wc_fst),
    !is.na(avg_wc_fst)
  )

dxy_hp2 <- dxy_hp %>%
  filter(
    is.finite(avg_dxy),
    !is.na(avg_dxy)
  )

cat("\n==============================\n")
cat("ROW COUNTS AFTER FILTERING\n")
cat("==============================\n")
cat("pi windows:", nrow(pi_all), "\n")
cat("Tajima's D windows:", nrow(tajd_all), "\n")
cat("Africa vs non-Africa FST windows:", nrow(fst_af2), "\n")
cat("Africa vs non-Africa dXY windows:", nrow(dxy_af2), "\n")
cat("Human vs Pig FST windows:", nrow(fst_hp2), "\n")
cat("Human vs Pig dXY windows:", nrow(dxy_hp2), "\n")

#========================================================
# 4. CHROMOSOME LABELS FOR PI BOXPLOT
#========================================================

chr_key <- tibble(
  chromosome = c(
    "CM024168.1","CM024169.1","CM024170.1","CM024171.1","CM024172.1",
    "CM024173.1","CM024174.1","CM024175.1","CM024176.1","CM024177.1",
    "CM024178.1","CM024179.1","CM024180.1","CM024181.1","CM024182.1",
    "CM024183.1","CM024184.1","CM024185.1","CM024186.1"
  ),
  chr_label = paste0("Chr", 1:19)
)

pi_all_chr <- pi_all %>%
  left_join(chr_key, by = "chromosome") %>%
  filter(!is.na(chr_label)) %>%
  mutate(
    chr_label = factor(chr_label, levels = paste0("Chr", 1:19))
  )

#========================================================
# 5. MAKE INDIVIDUAL PANELS
#========================================================

p_a <- ggplot(pi_all, aes(x = avg_pi)) +
  geom_histogram(
    bins = 70,
    fill = "#4DA3D9",
    colour = "black",
    linewidth = 0.2
  ) +
  theme_classic(base_size = 11) +
  labs(
    title = "A. Nucleotide diversity distribution",
    x = expression(pi),
    y = "Number of windows"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 11)
  )

p_b <- ggplot(tajd_all, aes(x = tajima_d)) +
  geom_histogram(
    bins = 70,
    fill = "#F28E1C",
    colour = "black",
    linewidth = 0.2
  ) +
  theme_classic(base_size = 11) +
  labs(
    title = "B. Tajima's D distribution",
    x = "Tajima's D",
    y = "Number of windows"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 11)
  )

p_c <- ggplot(fst_af2, aes(x = avg_wc_fst)) +
  geom_histogram(
    bins = 70,
    fill = "#8E4497",
    colour = "black",
    linewidth = 0.2
  ) +
  theme_classic(base_size = 11) +
  labs(
    title = "C. FST distribution (Africa vs non-Africa)",
    x = expression(F[ST]),
    y = "Number of windows"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 11)
  )

p_d <- ggplot(dxy_af2, aes(x = avg_dxy)) +
  geom_histogram(
    bins = 70,
    fill = "#43D17B",
    colour = "black",
    linewidth = 0.2
  ) +
  theme_classic(base_size = 11) +
  labs(
    title = "D. dXY distribution (Africa vs non-Africa)",
    x = expression(d[XY]),
    y = "Number of windows"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 11)
  )

p_f <- ggplot(fst_hp2, aes(x = avg_wc_fst)) +
  geom_histogram(
    bins = 70,
    fill = "#C77CFF",
    colour = "black",
    linewidth = 0.2
  ) +
  theme_classic(base_size = 11) +
  labs(
    title = "E. FST distribution (Human vs Pig)",
    x = expression(F[ST]),
    y = "Number of windows"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 11)
  )

p_g <- ggplot(dxy_hp2, aes(x = avg_dxy)) +
  geom_histogram(
    bins = 70,
    fill = "#00BFC4",
    colour = "black",
    linewidth = 0.2
  ) +
  theme_classic(base_size = 11) +
  labs(
    title = "F. dXY distribution (Human vs Pig)",
    x = expression(d[XY]),
    y = "Number of windows"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 11)
  )

p_e <- ggplot(pi_all_chr, aes(x = chr_label, y = avg_pi, fill = chr_label)) +
  geom_boxplot(
    outlier.size = 0.4,
    linewidth = 0.3
  ) +
  theme_classic(base_size = 11) +
  labs(
    title = "G. Chromosome-level nucleotide diversity",
    x = "Chromosome",
    y = expression(pi)
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

#========================================================
# 6. COMBINE PANELS
#========================================================

supp_fig <- (
  (p_a | p_b) /
    (p_c | p_d) /
    (p_f | p_g) /
    p_e
) +
  plot_layout(
    heights = c(1, 1, 1, 1.15)
  ) +
  plot_annotation(
    title = "Supplementary genome-wide diversity and differentiation summaries",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0
      )
    )
  )

print(supp_fig)

#========================================================
# 7. SAVE OUTPUTS
#========================================================

ggsave(
  filename = file.path(
    outdir,
    "supplementary_genome_wide_diversity_differentiation_with_human_pig.png"
  ),
  plot = supp_fig,
  width = 11,
  height = 13,
  dpi = 600
)

ggsave(
  filename = file.path(
    outdir,
    "supplementary_genome_wide_diversity_differentiation_with_human_pig.svg"
  ),
  plot = supp_fig,
  width = 11,
  height = 13
)

cat("\n==============================\n")
cat("OUTPUT FILES WRITTEN\n")
cat("==============================\n")
cat(
  file.path(
    outdir,
    "supplementary_genome_wide_diversity_differentiation_with_human_pig.png"
  ),
  "\n"
)
cat(
  file.path(
    outdir,
    "supplementary_genome_wide_diversity_differentiation_with_human_pig.svg"
  ),
  "\n"
```

---

## Outputs and QC checks

Top-1%/5% windows, nearby regions, merged intervals and candidate genes. Consult recorded logs and intermediate tables to distinguish expected settings from actual execution success.

## Provenance and limitations

- Some R blocks rely on pre-defined working directory/outdir; historical candidate coordinates are provisional pending PIXY rerun.
- Sample identifiers and host/country metadata were manually curated in `sample_QC_list.csv`; it is not generated by an algorithm in these files.
- Commands, versions, seeds and paths are retained where present; annotations are not evidence that all jobs have been freshly rerun.
- For the discovery and correction of the PIXY all-sites input issue, see [`99_PIXY_allsites_correction.md`](99_PIXY_allsites_correction.md).
