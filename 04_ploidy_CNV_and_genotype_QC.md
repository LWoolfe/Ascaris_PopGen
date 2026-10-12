# 04 — Ploidy/CNV and genotype-composition QC

## Purpose

Document the supplementary **chromosome-level relative copy-number (CN)** and **SNP genotype-composition** analyses performed for the *Ascaris* nuclear WGS dataset. The CN branch uses mapped reads rather than variant calls; the genotype-composition branch uses the final filtered **146-sample SNP VCF**. These are distinct downstream QC/description branches, not interchangeable ploidy estimates.

## Inputs

| Branch | Input | Origin |
|---|---|---|
| Relative CN: per-sample | `03_Mapping/3a_markdup/*.markdup.bam` | Mapping + `MarkDuplicates`, documented in `03_mapping_and_mapping_QC.md` |
| Relative CN: 25-kb windows | `01_Reference/25kb_mappability/Ascaris.autosomes19.25kb.unique90.bed` | Autosomal windows with ≥90% uniquely mappable bases |
| CN chromosome summary | `*.25kb.relative_CN.tsv`, `*.ploidy_summary.tsv` | Per-sample mosdepth job below |
| CN visualisation | `all_samples.depth5x.chromosome_CN.tsv`, `autosome_order.txt` | Chromosome summary job below |
| Genotype composition | `05b_146_data/Ascaris.autosomes.SNPs.146.maxmiss0.90.vcf.gz` | Final VCF documented in `05_variant_calling_and_SNP_QC.md` |
| Plot annotations | `sample_QC_list.csv` | **Manually curated** metadata; source file must be retained alongside code |

**Workflows were previously run.** All embedded code is for provenance/documentation, not a request to rerun analyses. SLURM job directives, sample counts and absolute paths are kept as used. Scripts below are separate jobs/analyses and should not be submitted as one combined shell file.

## Workflow

```text
156 marked-duplicate BAMs ── mosdepth (MAPQ ≥20; 25-kb ≥90% uniquely mappable windows)
             ↓
      depth per window ── median autosomal window depth per sample
             ↓
      relative_CN = 2 × window depth / sample median
             ↓
      chromosome median relative CN (≥5× main; ≥10× sensitivity)
             ↓
      heatmap by host and country

146-sample filtered SNP VCF ── bcftools stats / PSC genotype counts
             ↓
      per-sample genotype counts and percentages
             ↓
      QC of 146 rows and count formulae
             ↓
      country/host genotype composition figures
```

## Annotated scripts

### 04A. Per-sample relative CN: `ploidy25kb` array job

**Purpose.** For each of the 156 duplicate-marked BAM files, calculate mean depth in the autosomal 25-kb windows, then divide each window depth by the sample's median autosomal window depth and scale to diploid CN = 2.

**Key parameters.** `--array=1-156%20`, `mosdepth --mapq 20 --no-per-base`, 25-kb BED with ≥90% unique mappability. The script checks for the expected BAM count and warns (rather than stopping) if mosdepth outputs an unexpected number of windows.

**Status: USER-SUPPLIED ORIGINAL SCRIPT.** Formatting normalised; commands and thresholds preserved.

```bash
#!/bin/bash
#SBATCH --job-name=ploidy25kb
#SBATCH --output=/parallel_scratch/lw01165/03_Mapping/3d_25kb_ploidy/logs/ploidy25kb_%A_%a.out
#SBATCH --error=/parallel_scratch/lw01165/03_Mapping/3d_25kb_ploidy/logs/ploidy25kb_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-156%20

set -eo pipefail

echo "============================================================"
echo "25-kb autosomal depth / relative copy-number analysis"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "============================================================"

# ------------------------------------------------------------
# Environment
# ------------------------------------------------------------

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env

echo
echo "mosdepth:"
mosdepth --version

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

BAMDIR="/parallel_scratch/lw01165/03_Mapping/3a_markdup"
WINDOWS="/parallel_scratch/lw01165/01_Reference/25kb_mappability/Ascaris.autosomes19.25kb.unique90.bed"
OUTDIR="/parallel_scratch/lw01165/03_Mapping/3d_25kb_ploidy"
LOGDIR="${OUTDIR}/logs"

mkdir -p "${OUTDIR}" "${LOGDIR}"

# ------------------------------------------------------------
# Build sorted list of final markdup BAMs
# ------------------------------------------------------------

mapfile -t BAMS < <(
    find "${BAMDIR}" \
        -maxdepth 1 \
        -type f \
        -name "*.markdup.bam" \
        | sort
)

N_BAMS=${#BAMS[@]}

echo
echo "BAM files found: ${N_BAMS}"

if [[ "${N_BAMS}" -ne 156 ]]; then
    echo "ERROR: Expected 156 markdup BAMs but found ${N_BAMS}"
    exit 1
fi

# SLURM array is 1-based; bash arrays are 0-based
BAM="${BAMS[$((SLURM_ARRAY_TASK_ID - 1))]}"
SAMPLE=$(basename "${BAM}" .markdup.bam)

echo "Sample: ${SAMPLE}"
echo "BAM: ${BAM}"

# ------------------------------------------------------------
# Check BAM index
# ------------------------------------------------------------

if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
    echo "ERROR: BAM index not found:"
    echo "${BAM}"
    exit 1
fi

# ------------------------------------------------------------
# Check window BED
# ------------------------------------------------------------

if [[ ! -s "${WINDOWS}" ]]; then
    echo "ERROR: Window BED missing or empty:"
    echo "${WINDOWS}"
    exit 1
fi

EXPECTED_WINDOWS=$(wc -l < "${WINDOWS}")

echo
echo "Autosomal >=90% uniquely mappable windows: ${EXPECTED_WINDOWS}"

# ------------------------------------------------------------
# Run mosdepth
#
# --by             calculate mean depth within each BED region
# --no-per-base    don't produce enormous per-base output
# --mapq 20        retain reasonably well-mapped reads
#
# mosdepth automatically excludes reads marked as:
# unmapped, secondary, QC-fail and duplicate.
# ------------------------------------------------------------

PREFIX="${OUTDIR}/${SAMPLE}"

echo
echo "Calculating mean depth..."

mosdepth \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --by "${WINDOWS}" \
    --mapq 20 \
    --no-per-base \
    "${PREFIX}" \
    "${BAM}"

REGIONS="${PREFIX}.regions.bed.gz"

if [[ ! -s "${REGIONS}" ]]; then
    echo "ERROR: mosdepth regions output not produced."
    exit 1
fi

# ------------------------------------------------------------
# Save human-readable mean-depth table
# ------------------------------------------------------------

DEPTH_OUT="${OUTDIR}/${SAMPLE}.25kb.mean_depth.tsv"

{
    echo -e "chrom\tstart\tend\tmean_depth"
    zcat "${REGIONS}"
} > "${DEPTH_OUT}"

N_WINDOWS=$(zcat "${REGIONS}" | wc -l)

echo
echo "Windows returned by mosdepth: ${N_WINDOWS}"

if [[ "${N_WINDOWS}" -ne "${EXPECTED_WINDOWS}" ]]; then
    echo "WARNING: Expected ${EXPECTED_WINDOWS} windows but got ${N_WINDOWS}"
fi

# ------------------------------------------------------------
# Calculate median autosomal window depth
# ------------------------------------------------------------

DEPTH_TMP="${OUTDIR}/${SAMPLE}.depth_values.tmp"

zcat "${REGIONS}" \
    | awk '{print $4}' \
    | sort -n \
    > "${DEPTH_TMP}"

N=$(wc -l < "${DEPTH_TMP}")

if (( N % 2 == 1 )); then

    MID=$(( (N + 1) / 2 ))

    MEDIAN=$(sed -n "${MID}p" "${DEPTH_TMP}")

else

    MID1=$(( N / 2 ))
    MID2=$(( MID1 + 1 ))

    V1=$(sed -n "${MID1}p" "${DEPTH_TMP}")
    V2=$(sed -n "${MID2}p" "${DEPTH_TMP}")

    MEDIAN=$(awk -v a="${V1}" -v b="${V2}" \
        'BEGIN {printf "%.8f", (a+b)/2}')

fi

rm -f "${DEPTH_TMP}"

echo
echo "Median autosomal depth: ${MEDIAN}"

# Prevent division by zero
awk -v x="${MEDIAN}" 'BEGIN {
    if (x <= 0) exit 1
}' || {
    echo "ERROR: Median depth is zero or invalid."
    exit 1
}

# ------------------------------------------------------------
# Convert relative depth to estimated relative copy number
#
# Diploid autosomal baseline:
# relative_CN = 2 * window_depth / median_autosomal_depth
# ------------------------------------------------------------

CN_OUT="${OUTDIR}/${SAMPLE}.25kb.relative_CN.tsv"

{
    echo -e "chrom\tstart\tend\tmean_depth\tmedian_autosomal_depth\trelative_CN\tsample"

    zcat "${REGIONS}" \
        | awk -v median="${MEDIAN}" -v sample="${SAMPLE}" '
          BEGIN {OFS="\t"}
          {
              cn = 2 * ($4 / median)
              print $1,$2,$3,$4,median,cn,sample
          }'

} > "${CN_OUT}"

# ------------------------------------------------------------
# Sample-level summary
# ------------------------------------------------------------

SUMMARY="${OUTDIR}/${SAMPLE}.ploidy_summary.tsv"

{
    echo -e "sample\tn_windows\tmedian_autosomal_depth"
    echo -e "${SAMPLE}\t${N_WINDOWS}\t${MEDIAN}"
} > "${SUMMARY}"

# ------------------------------------------------------------
# Preview
# ------------------------------------------------------------

echo
echo "First five relative-CN windows:"
head -6 "${CN_OUT}"

echo
echo "Outputs:"
echo "${DEPTH_OUT}"
echo "${CN_OUT}"
echo "${SUMMARY}"

echo
echo "Finished: $(date)"
echo "============================================================"

```

**Outputs per sample:** `SAMPLE.25kb.mean_depth.tsv`, `SAMPLE.25kb.relative_CN.tsv`, `SAMPLE.ploidy_summary.tsv` (plus mosdepth native outputs).

### 04B. Across-sample chromosome CN summary: `chrCNsummary`

**Purpose.** Validate all 156 per-sample files, derive ≥5× and ≥10× sample sets, calculate median relative CN for each of 19 autosomes per sample and summarise exploratory CN outliers. The script expects 7,111 windows per sample.

**Status: USER-SUPPLIED ORIGINAL SCRIPT.** Do not confuse the exploratory thresholds `<1.5` / `>2.5` with a formal aneuploidy diagnosis.

**Caution:** the original job clears its `chromosome_CN_analysis/` output directory using `rm -rf`. Review before re-execution; retained verbatim to preserve history.

```bash
#!/bin/bash
#SBATCH --job-name=chrCNsummary
#SBATCH --output=/parallel_scratch/lw01165/03_Mapping/3d_25kb_ploidy/logs/chrCNsummary_%j.out
#SBATCH --error=/parallel_scratch/lw01165/03_Mapping/3d_25kb_ploidy/logs/chrCNsummary_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -eo pipefail

echo "============================================================"
echo "Ascaris chromosome-level relative copy-number analysis"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "============================================================"


# ============================================================
# 1. PATHS AND SETTINGS
# ============================================================

WORKDIR="/parallel_scratch/lw01165/03_Mapping/3d_25kb_ploidy"

WINDOWS="/parallel_scratch/lw01165/01_Reference/25kb_mappability/Ascaris.autosomes19.25kb.unique90.bed"

OUTDIR="${WORKDIR}/chromosome_CN_analysis"

EXPECTED_SAMPLES=156
EXPECTED_WINDOWS=7111
EXPECTED_CHROMOSOMES=19

# Main analysis threshold
PRIMARY_DEPTH=5

# More stringent sensitivity analysis
SENSITIVITY_DEPTH=10

# Exploratory chromosome-level CN flagging thresholds
LOW_CN=1.5
HIGH_CN=2.5


# ============================================================
# 2. CHECK REQUIRED INPUTS
# ============================================================

if [[ ! -d "${WORKDIR}" ]]; then
    echo "ERROR: Working directory not found:"
    echo "${WORKDIR}"
    exit 1
fi

if [[ ! -s "${WINDOWS}" ]]; then
    echo "ERROR: Window BED missing or empty:"
    echo "${WINDOWS}"
    exit 1
fi

if ! command -v gawk >/dev/null 2>&1; then
    echo "ERROR: gawk not found."
    exit 1
fi


# ============================================================
# 3. CREATE CLEAN OUTPUT DIRECTORY
# ============================================================

rm -rf "${OUTDIR}"
mkdir -p "${OUTDIR}"

cd "${WORKDIR}"


echo
echo "Analysis parameters:"
echo "  Expected samples:          ${EXPECTED_SAMPLES}"
echo "  Autosomal windows/sample:  ${EXPECTED_WINDOWS}"
echo "  Autosomes:                 ${EXPECTED_CHROMOSOMES}"
echo "  Window size:               25 kb"
echo "  Unique mappability:        >=90%"
echo "  Primary depth threshold:   >=${PRIMARY_DEPTH}x"
echo "  Sensitivity threshold:     >=${SENSITIVITY_DEPTH}x"
echo "  Exploratory CN flag:       <${LOW_CN} or >${HIGH_CN}"
echo


# ============================================================
# 4. AUTOSOME ORDER
#
# Preserve chromosome order as it appears in the BED file.
# ============================================================

AUTOSOME_ORDER="${OUTDIR}/autosome_order.txt"

awk '!seen[$1]++ {print $1}' "${WINDOWS}" \
    > "${AUTOSOME_ORDER}"

N_AUTOSOMES=$(wc -l < "${AUTOSOME_ORDER}")

echo "Autosomes detected in BED: ${N_AUTOSOMES}"

if [[ "${N_AUTOSOMES}" -ne "${EXPECTED_CHROMOSOMES}" ]]; then
    echo "ERROR: Expected ${EXPECTED_CHROMOSOMES} autosomes."
    exit 1
fi


# ============================================================
# 5. VERIFY INDIVIDUAL MOSDEPTH OUTPUTS
# ============================================================

N_SUMMARY=$(find "${WORKDIR}" \
    -maxdepth 1 \
    -type f \
    -name "*.ploidy_summary.tsv" \
    ! -name "all_samples*" \
    | wc -l)

N_CN=$(find "${WORKDIR}" \
    -maxdepth 1 \
    -type f \
    -name "*.25kb.relative_CN.tsv" \
    | wc -l)

echo
echo "Input files:"
echo "  Individual sample summaries: ${N_SUMMARY}"
echo "  Relative-CN files:           ${N_CN}"

if [[ "${N_SUMMARY}" -ne "${EXPECTED_SAMPLES}" ]]; then
    echo "ERROR: Expected ${EXPECTED_SAMPLES} individual summary files."
    exit 1
fi

if [[ "${N_CN}" -ne "${EXPECTED_SAMPLES}" ]]; then
    echo "ERROR: Expected ${EXPECTED_SAMPLES} relative-CN files."
    exit 1
fi


# ============================================================
# 6. REBUILD COMBINED SAMPLE DEPTH SUMMARY
#
# This avoids relying on manually created combined files.
# ============================================================

COMBINED_SUMMARY="${OUTDIR}/all_samples.ploidy_summary.tsv"

echo -e \
"sample\tn_windows\tmedian_autosomal_depth" \
> "${COMBINED_SUMMARY}"


while IFS= read -r FILE
do

    awk 'NR==2 {print}' "${FILE}" \
        >> "${COMBINED_SUMMARY}"

done < <(
    find "${WORKDIR}" \
        -maxdepth 1 \
        -type f \
        -name "*.ploidy_summary.tsv" \
        ! -name "all_samples*" \
        | sort
)


N_COMBINED=$(( $(wc -l < "${COMBINED_SUMMARY}") - 1 ))

echo
echo "Combined sample summary: ${N_COMBINED} samples"

if [[ "${N_COMBINED}" -ne "${EXPECTED_SAMPLES}" ]]; then
    echo "ERROR: Combined summary does not contain ${EXPECTED_SAMPLES} samples."
    exit 1
fi


# ============================================================
# 7. COVERAGE DISTRIBUTION
# ============================================================

COVERAGE_SUMMARY="${OUTDIR}/coverage_summary.txt"

awk '
NR > 1 {

    n++
    sum += $3

    if ($3 < 1)  lt1++
    if ($3 < 2)  lt2++
    if ($3 < 3)  lt3++
    if ($3 < 5)  lt5++
    if ($3 < 10) lt10++
}

END {

    print "Samples:", n
    printf "Mean median autosomal depth: %.4f\n", sum/n
    print "<1x:",  lt1+0
    print "<2x:",  lt2+0
    print "<3x:",  lt3+0
    print "<5x:",  lt5+0
    print "<10x:", lt10+0
}
' "${COMBINED_SUMMARY}" \
> "${COVERAGE_SUMMARY}"


# ============================================================
# 8. CREATE >=5x AND >=10x SAMPLE LISTS
# ============================================================

awk -v d="${PRIMARY_DEPTH}" '
NR > 1 && $3 >= d {
    print $1
}
' "${COMBINED_SUMMARY}" \
> "${OUTDIR}/samples.depth5x.txt"


awk -v d="${SENSITIVITY_DEPTH}" '
NR > 1 && $3 >= d {
    print $1
}
' "${COMBINED_SUMMARY}" \
> "${OUTDIR}/samples.depth10x.txt"


N5=$(wc -l < "${OUTDIR}/samples.depth5x.txt")
N10=$(wc -l < "${OUTDIR}/samples.depth10x.txt")


# ============================================================
# 9. RECORD EXCLUDED SAMPLES
# ============================================================

{
    echo -e "sample\tmedian_autosomal_depth"

    awk -v d="${PRIMARY_DEPTH}" '
    NR > 1 && $3 < d {
        print $1 "\t" $3
    }
    ' "${COMBINED_SUMMARY}"

} > "${OUTDIR}/samples.excluded_below5x.tsv"


{
    echo -e "sample\tmedian_autosomal_depth"

    awk -v d="${SENSITIVITY_DEPTH}" '
    NR > 1 && $3 < d {
        print $1 "\t" $3
    }
    ' "${COMBINED_SUMMARY}"

} > "${OUTDIR}/samples.excluded_below10x.tsv"


echo
echo "Samples retained:"
echo "  >=5x:  ${N5}"
echo "  >=10x: ${N10}"


# ============================================================
# 10. FUNCTION TO CALCULATE MEDIAN CN PER CHROMOSOME
#
# Each sample already contains:
#
# chrom
# start
# end
# mean_depth
# median_autosomal_depth
# relative_CN
# sample
#
# For each chromosome we take the MEDIAN relative_CN across
# all retained 25-kb windows.
# ============================================================

make_chr_summary () {

    SAMPLE_LIST="$1"
    OUTPUT="$2"

    echo -e \
"sample\tchrom\tn_windows\tmedian_relative_CN\tmedian_autosomal_depth" \
    > "${OUTPUT}"

    while read -r SAMPLE
    do

        FILE="${WORKDIR}/${SAMPLE}.25kb.relative_CN.tsv"
        SUMMARY_FILE="${WORKDIR}/${SAMPLE}.ploidy_summary.tsv"

        if [[ ! -s "${FILE}" ]]; then
            echo "ERROR: Missing CN file for ${SAMPLE}" >&2
            exit 1
        fi

        if [[ ! -s "${SUMMARY_FILE}" ]]; then
            echo "ERROR: Missing summary file for ${SAMPLE}" >&2
            exit 1
        fi

        N_FILE_WINDOWS=$(( $(wc -l < "${FILE}") - 1 ))

        if [[ "${N_FILE_WINDOWS}" -ne "${EXPECTED_WINDOWS}" ]]; then
            echo "ERROR: ${SAMPLE} contains ${N_FILE_WINDOWS} windows; expected ${EXPECTED_WINDOWS}" >&2
            exit 1
        fi

        SAMPLE_DEPTH=$(awk 'NR==2 {print $3}' "${SUMMARY_FILE}")

        gawk \
            -v sample="${SAMPLE}" \
            -v depth="${SAMPLE_DEPTH}" '

        BEGIN {
            OFS="\t"
        }

        NR == 1 {
            next
        }

        {
            chr=$1

            # Record chromosome order based on first appearance
            if (!(chr in seen_chr)) {
                seen_chr[chr]=1
                chr_order[++nchr]=chr
            }

            value[chr, ++n[chr]]=$6
        }

        END {

            for (c=1; c<=nchr; c++) {

                chr=chr_order[c]

                delete x

                for (i=1; i<=n[chr]; i++) {
                    x[i]=value[chr,i]
                }

                m=asort(x, x, "@val_num_asc")

                if (m % 2 == 1) {
                    med=x[(m+1)/2]
                }
                else {
                    med=(x[m/2] + x[m/2+1])/2
                }

                print sample,chr,m,med,depth
            }
        }
        ' "${FILE}" >> "${OUTPUT}"

    done < "${SAMPLE_LIST}"
}


# ============================================================
# 11. GENERATE CHROMOSOME-LEVEL CN TABLES
# ============================================================

echo
echo "Generating >=5x chromosome-level CN table..."

make_chr_summary \
    "${OUTDIR}/samples.depth5x.txt" \
    "${OUTDIR}/all_samples.depth5x.chromosome_CN.tsv"


echo "Generating >=10x chromosome-level CN table..."

make_chr_summary \
    "${OUTDIR}/samples.depth10x.txt" \
    "${OUTDIR}/all_samples.depth10x.chromosome_CN.tsv"


# ============================================================
# 12. VALIDATE CHROMOSOME TABLES
# ============================================================

for DEPTH in 5 10
do

    FILE="${OUTDIR}/all_samples.depth${DEPTH}x.chromosome_CN.tsv"

    N_SAMPLES=$(awk 'NR>1 {print $1}' "${FILE}" \
        | sort -u \
        | wc -l)

    N_CHR=$(awk 'NR>1 {print $2}' "${FILE}" \
        | sort -u \
        | wc -l)

    N_ROWS=$(( $(wc -l < "${FILE}") - 1 ))

    EXPECTED_FOR_FILTER=$(( N_SAMPLES * EXPECTED_CHROMOSOMES ))

    echo
    echo ">=${DEPTH}x table QC:"
    echo "  Samples:            ${N_SAMPLES}"
    echo "  Chromosomes:        ${N_CHR}"
    echo "  Data rows:          ${N_ROWS}"
    echo "  Expected data rows: ${EXPECTED_FOR_FILTER}"

    if [[ "${N_CHR}" -ne "${EXPECTED_CHROMOSOMES}" ]]; then
        echo "ERROR: Expected ${EXPECTED_CHROMOSOMES} autosomes."
        exit 1
    fi

    if [[ "${N_ROWS}" -ne "${EXPECTED_FOR_FILTER}" ]]; then
        echo "ERROR: Unexpected number of chromosome-level rows."
        exit 1
    fi

done


# ============================================================
# 13. IDENTIFY EXPLORATORY CHROMOSOME CN OUTLIERS
#
# These are FLAGGING thresholds only, not formal biological
# definitions of aneuploidy.
# ============================================================

for DEPTH in 5 10
do

    INPUT="${OUTDIR}/all_samples.depth${DEPTH}x.chromosome_CN.tsv"

    OUTPUT="${OUTDIR}/depth${DEPTH}x.chromosome_CN_outliers.tsv"

    awk \
        -v low="${LOW_CN}" \
        -v high="${HIGH_CN}" '

    BEGIN {
        OFS="\t"
    }

    NR == 1 {
        print
        next
    }

    $4 < low || $4 > high {
        print
    }

    ' "${INPUT}" > "${OUTPUT}"

done


# ============================================================
# 14. RANK LOWEST/HIGHEST CHROMOSOME CN VALUES
#
# IMPORTANT:
# We use sed rather than head after sort so pipefail does not
# terminate the job because of SIGPIPE.
# ============================================================

for DEPTH in 5 10
do

    INPUT="${OUTDIR}/all_samples.depth${DEPTH}x.chromosome_CN.tsv"


    # Lowest 20

    TEMP_LOW="${OUTDIR}/depth${DEPTH}x.lowest.sorted.tmp"

    awk 'NR>1 {
        print $4 "\t" $1 "\t" $2 "\t" $5
    }' "${INPUT}" \
    | sort -k1,1n \
    > "${TEMP_LOW}"

    {
        echo -e \
"median_relative_CN\tsample\tchrom\tmedian_autosomal_depth"

        sed -n '1,20p' "${TEMP_LOW}"

    } > "${OUTDIR}/depth${DEPTH}x.lowest_20_chromosome_CN.tsv"

    rm -f "${TEMP_LOW}"


    # Highest 20

    TEMP_HIGH="${OUTDIR}/depth${DEPTH}x.highest.sorted.tmp"

    awk 'NR>1 {
        print $4 "\t" $1 "\t" $2 "\t" $5
    }' "${INPUT}" \
    | sort -k1,1nr \
    > "${TEMP_HIGH}"

    {
        echo -e \
"median_relative_CN\tsample\tchrom\tmedian_autosomal_depth"

        sed -n '1,20p' "${TEMP_HIGH}"

    } > "${OUTDIR}/depth${DEPTH}x.highest_20_chromosome_CN.tsv"

    rm -f "${TEMP_HIGH}"

done


# ============================================================
# 15. DATASET-WIDE CN RANGE SUMMARY
# ============================================================

RANGE_OUT="${OUTDIR}/chromosome_CN_ranges.tsv"

echo -e \
"depth_filter\tn_samples\tn_chromosome_values\tmin_CN\tmean_CN\tmax_CN" \
> "${RANGE_OUT}"


for DEPTH in 5 10
do

    INPUT="${OUTDIR}/all_samples.depth${DEPTH}x.chromosome_CN.tsv"

    awk -v label=">=${DEPTH}x" '

    NR > 1 {

        n++
        sum += $4

        if (n == 1 || $4 < min)
            min=$4

        if (n == 1 || $4 > max)
            max=$4

        sample[$1]=1
    }

    END {

        ns=0

        for (s in sample)
            ns++

        printf "%s\t%d\t%d\t%.6f\t%.6f\t%.6f\n",
               label,ns,n,min,sum/n,max
    }

    ' "${INPUT}" >> "${RANGE_OUT}"

done


# ============================================================
# 16. PER-CHROMOSOME SUMMARY ACROSS SAMPLES
#
# Useful for identifying chromosomes that systematically have
# higher/lower relative depth across the dataset.
# ============================================================

for DEPTH in 5 10
do

    INPUT="${OUTDIR}/all_samples.depth${DEPTH}x.chromosome_CN.tsv"

    OUTPUT="${OUTDIR}/depth${DEPTH}x.per_chromosome_summary.tsv"

    gawk \
        -v low="${LOW_CN}" \
        -v high="${HIGH_CN}" '

    BEGIN {
        OFS="\t"

        print "chrom",
              "n_samples",
              "mean_relative_CN",
              "median_relative_CN",
              "min_relative_CN",
              "max_relative_CN",
              "n_flagged"
    }

    NR == 1 {
        next
    }

    {
        chr=$2
        cn=$4

        if (!(chr in seen)) {
            seen[chr]=1
            order[++nchr]=chr
        }

        value[chr, ++n[chr]]=cn
        sum[chr]+=cn

        if (n[chr] == 1 || cn < min[chr])
            min[chr]=cn

        if (n[chr] == 1 || cn > max[chr])
            max[chr]=cn

        if (cn < low || cn > high)
            flagged[chr]++
    }

    END {

        for (c=1; c<=nchr; c++) {

            chr=order[c]

            delete x

            for (i=1; i<=n[chr]; i++)
                x[i]=value[chr,i]

            m=asort(x, x, "@val_num_asc")

            if (m % 2 == 1)
                med=x[(m+1)/2]
            else
                med=(x[m/2] + x[m/2+1])/2

            print chr,
                  n[chr],
                  sum[chr]/n[chr],
                  med,
                  min[chr],
                  max[chr],
                  flagged[chr]+0
        }
    }

    ' "${INPUT}" > "${OUTPUT}"

done


# ============================================================
# 17. COUNT FLAGGED VALUES AND SAMPLES
# ============================================================

N_OUT5=$(( $(wc -l < \
    "${OUTDIR}/depth5x.chromosome_CN_outliers.tsv") - 1 ))

N_OUT10=$(( $(wc -l < \
    "${OUTDIR}/depth10x.chromosome_CN_outliers.tsv") - 1 ))


N_OUT_SAMPLES5=$(awk 'NR>1 {print $1}' \
    "${OUTDIR}/depth5x.chromosome_CN_outliers.tsv" \
    | sort -u \
    | wc -l)

N_OUT_SAMPLES10=$(awk 'NR>1 {print $1}' \
    "${OUTDIR}/depth10x.chromosome_CN_outliers.tsv" \
    | sort -u \
    | wc -l)


# ============================================================
# 18. WRITE ANALYSIS PARAMETERS FILE
# ============================================================

PARAMETERS="${OUTDIR}/analysis_parameters.tsv"

{
    echo -e "parameter\tvalue"
    echo -e "window_size_bp\t25000"
    echo -e "minimum_unique_mappability\t0.90"
    echo -e "expected_autosomes\t${EXPECTED_CHROMOSOMES}"
    echo -e "expected_windows_per_sample\t${EXPECTED_WINDOWS}"
    echo -e "total_samples\t${EXPECTED_SAMPLES}"
    echo -e "primary_minimum_depth\t${PRIMARY_DEPTH}"
    echo -e "sensitivity_minimum_depth\t${SENSITIVITY_DEPTH}"
    echo -e "exploratory_low_CN_flag\t${LOW_CN}"
    echo -e "exploratory_high_CN_flag\t${HIGH_CN}"
    echo -e "relative_CN_formula\t2_x_window_depth_divided_by_sample_median_autosomal_depth"
    echo -e "chromosome_summary_statistic\tmedian_relative_CN"

} > "${PARAMETERS}"


# ============================================================
# 19. FINAL REPORT
# ============================================================

echo
echo "============================================================"
echo "FINAL ANALYSIS SUMMARY"
echo "============================================================"

echo
cat "${COVERAGE_SUMMARY}"

echo
echo "------------------------------------------------------------"
echo "Primary analysis (>=5x)"
echo "------------------------------------------------------------"
echo "Samples retained: ${N5}/${EXPECTED_SAMPLES}"
echo "Chromosome values outside ${LOW_CN}-${HIGH_CN}: ${N_OUT5}"
echo "Samples containing >=1 flagged chromosome: ${N_OUT_SAMPLES5}"

echo
echo "------------------------------------------------------------"
echo "Sensitivity analysis (>=10x)"
echo "------------------------------------------------------------"
echo "Samples retained: ${N10}/${EXPECTED_SAMPLES}"
echo "Chromosome values outside ${LOW_CN}-${HIGH_CN}: ${N_OUT10}"
echo "Samples containing >=1 flagged chromosome: ${N_OUT_SAMPLES10}"

echo
echo "------------------------------------------------------------"
echo "Chromosome CN ranges"
echo "------------------------------------------------------------"
cat "${RANGE_OUT}"

echo
echo "------------------------------------------------------------"
echo "Samples excluded from primary >=5x analysis"
echo "------------------------------------------------------------"
cat "${OUTDIR}/samples.excluded_below5x.tsv"

echo
echo "------------------------------------------------------------"
echo "Highest chromosome-level CN values (>=5x)"
echo "------------------------------------------------------------"
sed -n '1,6p' \
    "${OUTDIR}/depth5x.highest_20_chromosome_CN.tsv"

echo
echo "------------------------------------------------------------"
echo "Lowest chromosome-level CN values (>=5x)"
echo "------------------------------------------------------------"
sed -n '1,6p' \
    "${OUTDIR}/depth5x.lowest_20_chromosome_CN.tsv"

echo
echo "------------------------------------------------------------"
echo "Outputs written to:"
echo "${OUTDIR}"
echo "------------------------------------------------------------"

echo
echo "Finished: $(date)"
echo "============================================================"
```

**Outputs:** `all_samples.depth5x.chromosome_CN.tsv`, `all_samples.depth10x.chromosome_CN.tsv`, per-chromosome summaries, depth-based sample lists, QC reports, and exploratory CN outliers.

### 04C. Per-sample SNP genotype counts: `SNP_genotype_counts`

**Purpose.** Generate per-sample genotype counts from the final 146-sample filtered SNP VCF using `bcftools stats -s -`. Parse the `PSC` sample-level records and calculate the proportion of called SNP loci that are heterozygous or homozygous alternative.

**Status: USER-SUPPLIED ORIGINAL SCRIPT (RECOVERED).** Previously recorded only as a partial command; the full original SLURM job has now been supplied and is retained below. Formatting has been normalised without changing the calculation or thresholds.

**Inputs:** `05b_146_data/Ascaris.autosomes.SNPs.146.maxmiss0.90.vcf.gz` (146 samples).

**Outputs:** `Ascaris.146.bcftools_stats.txt`, `Ascaris.146.per_sample_genotype_counts.tsv` and `Ascaris.146.genotype_summary.txt`.

**QC note:** The original script identifies `PSC` columns by their positions (not by parsing the `bcftools stats` descriptive header). When reproducing, confirm these positions match the installed version's `PSC` header. The independent `QC_genotypes` validation follows in 04D.

```bash
#!/bin/bash
#SBATCH --job-name=SNP_genotype_counts
#SBATCH --output=SNP_genotype_counts_%j.out
#SBATCH --error=SNP_genotype_counts_%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

set -eo pipefail

# ============================================================
# Per-sample heterozygous and homozygous-alt SNP counts
# ============================================================

VCF="/parallel_scratch/lw01165/05_analysis/05b_146_data/Ascaris.autosomes.SNPs.146.maxmiss0.90.vcf.gz"
OUTDIR="/parallel_scratch/lw01165/05_analysis/05b_146_data/genotype_counts"

mkdir -p "${OUTDIR}"

# Activate environment containing bcftools
source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env

echo "============================================================"
echo "Per-sample SNP genotype counts"
echo "Started: $(date)"
echo "VCF: ${VCF}"
echo "============================================================"

bcftools --version | head -1

# ============================================================
# 1. Run bcftools stats
# ============================================================

bcftools stats \
    -s - \
    "${VCF}" \
    > "${OUTDIR}/Ascaris.146.bcftools_stats.txt"

# ============================================================
# 2. Extract per-sample genotype counts
#
# PSC columns from bcftools stats:
#
# $3 = sample
# $4 = nRefHom
# $5 = nNonRefHom
# $6 = nHets
# $14 = nMissing
#
# Since the input contains SNPs only:
#
# called loci = nRefHom + nNonRefHom + nHets
# heterozygous percentage = 100 * nHets / called loci
# homozygous-alt percentage = 100 * nNonRefHom / called loci
# ============================================================

awk 'BEGIN {
    OFS="\t";

    print "sample",
          "n_ref_hom",
          "n_hom_alt",
          "n_het",
          "n_missing",
          "n_called",
          "het_percent_called",
          "hom_alt_percent_called"
}

$1=="PSC" {

    sample=$3;
    refhom=$4;
    homalt=$5;
    het=$6;
    missing=$14;

    called=refhom+homalt+het;

    if (called > 0) {
        het_pct=(het/called)*100;
        homalt_pct=(homalt/called)*100;
    }
    else {
        het_pct="NA";
        homalt_pct="NA";
    }

    print sample,
          refhom,
          homalt,
          het,
          missing,
          called,
          het_pct,
          homalt_pct
}' "${OUTDIR}/Ascaris.146.bcftools_stats.txt" \
> "${OUTDIR}/Ascaris.146.per_sample_genotype_counts.tsv"

# ============================================================
# 3. Dataset summary
# ============================================================

awk 'BEGIN {
    FS=OFS="\t"
}

NR==1 {
    next
}

{
    samples++;
    het_sum += $4;
    homalt_sum += $3;

    if (samples==1 || $4 < het_min) {
        het_min=$4
    }
    if (samples==1 || $4 > het_max) {
        het_max=$4
    }
    if (samples==1 || $3 < homalt_min) {
        homalt_min=$3
    }
    if (samples==1 || $3 > homalt_max) {
        homalt_max=$3
    }
    if (samples==1 || $7 < hetpct_min) {
        hetpct_min=$7
    }
    if (samples==1 || $7 > hetpct_max) {
        hetpct_max=$7
    }
    hetpct_sum += $7;
}

END {
    print "Samples:", samples;
    print "Heterozygous SNP count range:", het_min, "-", het_max;
    print "Mean heterozygous SNP count:", het_sum/samples;
    print "Homozygous-alt SNP count range:", homalt_min, "-", homalt_max;
    print "Mean homozygous-alt SNP count:", homalt_sum/samples;
    print "Heterozygous genotype percentage range:", hetpct_min, "-", hetpct_max;
    print "Mean heterozygous genotype percentage:", hetpct_sum/samples;
}' "${OUTDIR}/Ascaris.146.per_sample_genotype_counts.tsv" \
> "${OUTDIR}/Ascaris.146.genotype_summary.txt"

# ============================================================
# 4. Print summary
# ============================================================

echo
echo "============================================================"
echo "SUMMARY"
echo "============================================================"

cat "${OUTDIR}/Ascaris.146.genotype_summary.txt"

echo
echo "First 10 samples:"
column -t \
    "${OUTDIR}/Ascaris.146.per_sample_genotype_counts.tsv" \
    | head -11 || true

echo
echo "Finished: $(date)"
echo "Output directory: ${OUTDIR}"
echo "============================================================"
```

### 04D. Per-sample genotype-count QA: `QC_genotypes`

**Purpose.** Audit existing `bcftools` PSC records and genotype counts; verify 146 samples, unique sample IDs, and count/percentage formulae. It does **not** rerun `bcftools stats`.

**Status: USER-SUPPLIED ORIGINAL SCRIPT.**

```bash
#!/bin/bash
#SBATCH --job-name=QC_genotypes
#SBATCH --output=QC_genotypes_%j.out
#SBATCH --error=QC_genotypes_%j.err
#SBATCH --time=00:20:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1

set -eo pipefail

# ============================================================
# QC and summarise per-sample SNP genotype counts
#
# This job DOES NOT rerun bcftools stats.
#
# It checks the existing outputs from the previous job and
# regenerates the dataset-level genotype summary.
# ============================================================


# ============================================================
# 1. INPUTS
# ============================================================

OUTDIR="/parallel_scratch/lw01165/05_analysis/05b_146_data/genotype_counts"

STATS="${OUTDIR}/Ascaris.146.bcftools_stats.txt"

COUNTS="${OUTDIR}/Ascaris.146.per_sample_genotype_counts.tsv"

SUMMARY="${OUTDIR}/Ascaris.146.genotype_summary.txt"

QC_REPORT="${OUTDIR}/Ascaris.146.genotype_QC.txt"

EXPECTED_SAMPLES=146


echo "============================================================"
echo "Ascaris SNP genotype-count QC"
echo "Started: $(date)"
echo "============================================================"
echo


# ============================================================
# 2. CHECK REQUIRED FILES EXIST
# ============================================================

echo "Checking input files..."

if [[ ! -s "${STATS}" ]]; then
    echo "ERROR: Missing or empty file:"
    echo "${STATS}"
    exit 1
fi

if [[ ! -s "${COUNTS}" ]]; then
    echo "ERROR: Missing or empty file:"
    echo "${COUNTS}"
    exit 1
fi

echo "PASS: bcftools stats file exists and is non-empty."
echo "PASS: per-sample genotype-count file exists and is non-empty."
echo


# ============================================================
# 3. CHECK NUMBER OF PSC RECORDS IN BCFTOOLS STATS
# ============================================================

N_PSC=$(grep -c $'^PSC\t' "${STATS}" || true)

echo "Per-sample PSC records in bcftools stats: ${N_PSC}"

if [[ "${N_PSC}" -ne "${EXPECTED_SAMPLES}" ]]; then
    echo "ERROR: Expected ${EXPECTED_SAMPLES} PSC records but found ${N_PSC}."
    exit 1
fi

echo "PASS: bcftools stats contains ${EXPECTED_SAMPLES} samples."
echo


# ============================================================
# 4. CHECK PER-SAMPLE TSV HEADER
# ============================================================

EXPECTED_HEADER=$'sample\tn_ref_hom\tn_hom_alt\tn_het\tn_missing\tn_called\thet_percent_called\thom_alt_percent_called'

OBSERVED_HEADER=$(head -n 1 "${COUNTS}")

if [[ "${OBSERVED_HEADER}" != "${EXPECTED_HEADER}" ]]; then

    echo "ERROR: Unexpected header in:"
    echo "${COUNTS}"
    echo
    echo "Observed:"
    echo "${OBSERVED_HEADER}"
    echo
    echo "Expected:"
    echo "${EXPECTED_HEADER}"

    exit 1
fi

echo "PASS: per-sample TSV header is correct."
echo


# ============================================================
# 5. CHECK NUMBER OF SAMPLE ROWS
# ============================================================

N_SAMPLES=$(awk 'NR > 1 {n++} END {print n+0}' "${COUNTS}")

echo "Samples in per-sample TSV: ${N_SAMPLES}"

if [[ "${N_SAMPLES}" -ne "${EXPECTED_SAMPLES}" ]]; then
    echo "ERROR: Expected ${EXPECTED_SAMPLES} samples but found ${N_SAMPLES}."
    exit 1
fi

echo "PASS: expected sample count recovered."
echo


# ============================================================
# 6. CHECK FOR DUPLICATE SAMPLE IDS
# ============================================================

DUPLICATES=$(awk -F'\t' '
NR > 1 {
    count[$1]++
}
END {
    for (sample in count) {
        if (count[sample] > 1) {
            print sample "\t" count[sample]
        }
    }
}' "${COUNTS}")


if [[ -n "${DUPLICATES}" ]]; then

    echo "ERROR: Duplicate sample IDs detected:"
    echo "${DUPLICATES}"
    exit 1

fi

echo "PASS: no duplicate sample IDs."
echo


# ============================================================
# 7. CHECK TABLE STRUCTURE AND CALCULATIONS
#
# Columns:
#
# 1 sample
# 2 n_ref_hom
# 3 n_hom_alt
# 4 n_het
# 5 n_missing
# 6 n_called
# 7 het_percent_called
# 8 hom_alt_percent_called
#
# Checks:
#
# n_called = n_ref_hom + n_hom_alt + n_het
#
# het % =
#     100 * n_het / n_called
#
# hom-alt % =
#     100 * n_hom_alt / n_called
# ============================================================

echo "Checking genotype calculations..."

awk -F'\t' '
BEGIN {
    tolerance = 0.001
    bad = 0
}

NR == 1 {
    next
}

{
    # Check number of columns
    if (NF != 8) {
        print "ERROR: " $1 " has " NF " columns; expected 8." > "/dev/stderr"
        bad = 1
    }

    refhom  = $2 + 0
    homalt  = $3 + 0
    het     = $4 + 0
    missing = $5 + 0
    called  = $6 + 0

    # Genotype counts must not be negative
    if (refhom < 0 || homalt < 0 || het < 0 || missing < 0 || called < 0) {
        print "ERROR: Negative genotype count for " $1 > "/dev/stderr"
        bad = 1
    }

    # Called count must equal genotype categories
    expected_called = refhom + homalt + het

    if (called != expected_called) {
        print "ERROR: n_called mismatch for " $1 \
              ": reported=" called \
              ", calculated=" expected_called > "/dev/stderr"
        bad = 1
    }

    if (called <= 0) {
        print "ERROR: Zero successfully called loci for " $1 > "/dev/stderr"
        bad = 1
        next
    }

    # Recalculate percentages
    expected_het_pct = 100 * het / called
    expected_homalt_pct = 100 * homalt / called

    diff1 = $7 - expected_het_pct
    if (diff1 < 0) diff1 = -diff1

    diff2 = $8 - expected_homalt_pct
    if (diff2 < 0) diff2 = -diff2

    if (diff1 > tolerance) {
        print "ERROR: Heterozygous percentage mismatch for " $1 > "/dev/stderr"
        bad = 1
    }

    if (diff2 > tolerance) {
        print "ERROR: Homozygous-alt percentage mismatch for " $1 > "/dev/stderr"
        bad = 1
    }
}

END {
    exit bad
}' "${COUNTS}"

echo "PASS: n_called = n_ref_hom + n_hom_alt + n_het for all samples."
echo "PASS: heterozygous percentages recalculate correctly."
echo "PASS: homozygous-alt percentages recalculate correctly."
echo


# ============================================================
# 8. CHECK TOTAL NUMBER OF SNP LOCI IS CONSISTENT
#
# n_called + n_missing should describe the same SNP dataset
# for every sample.
# ============================================================

echo "Checking total SNP loci represented per sample..."

N_TOTAL_VALUES=$(awk -F'\t' '
NR > 1 {
    total = $5 + $6
    seen[total] = 1
}
END {
    n = 0
    for (x in seen) {
        n++
    }
    print n
}' "${COUNTS}")


if [[ "${N_TOTAL_VALUES}" -ne 1 ]]; then

    echo "WARNING: n_called + n_missing is not identical across all samples."
    echo "Inspect the per-sample table before interpretation."

else

    TOTAL_SNPS=$(awk -F'\t' '
    NR == 2 {
        print $5 + $6
        exit
    }' "${COUNTS}")

    echo "PASS: all samples represent the same ${TOTAL_SNPS} SNP loci."

fi

echo


# ============================================================
# 9. GENERATE DATASET-LEVEL SUMMARY
#
# Portable awk syntax used here to avoid the previous
# multiline-if problem.
# ============================================================

echo "Generating genotype summary..."

awk -F'\t' '
BEGIN {
    OFS = "\t"
}

NR == 1 {
    next
}

NR == 2 {

    samples = 1

    het_sum = $4
    homalt_sum = $3
    refhom_sum = $2

    hetpct_sum = $7
    homaltpct_sum = $8

    het_min = het_max = $4
    homalt_min = homalt_max = $3
    refhom_min = refhom_max = $2

    hetpct_min = hetpct_max = $7
    homaltpct_min = homaltpct_max = $8

    called_min = called_max = $6

    next
}

NR > 2 {

    samples++

    het_sum += $4
    homalt_sum += $3
    refhom_sum += $2

    hetpct_sum += $7
    homaltpct_sum += $8

    if ($4 < het_min) het_min = $4
    if ($4 > het_max) het_max = $4

    if ($3 < homalt_min) homalt_min = $3
    if ($3 > homalt_max) homalt_max = $3

    if ($2 < refhom_min) refhom_min = $2
    if ($2 > refhom_max) refhom_max = $2

    if ($7 < hetpct_min) hetpct_min = $7
    if ($7 > hetpct_max) hetpct_max = $7

    if ($8 < homaltpct_min) homaltpct_min = $8
    if ($8 > homaltpct_max) homaltpct_max = $8

    if ($6 < called_min) called_min = $6
    if ($6 > called_max) called_max = $6
}

END {

    print "Samples:", samples
    print ""

    print "Successfully called SNP loci range:", called_min, "-", called_max
    print ""

    print "Heterozygous SNP count range:", het_min, "-", het_max
    print "Mean heterozygous SNP count:", het_sum / samples
    print ""

    print "Homozygous-alt SNP count range:", homalt_min, "-", homalt_max
    print "Mean homozygous-alt SNP count:", homalt_sum / samples
    print ""

    print "Homozygous-reference SNP-locus count range:", refhom_min, "-", refhom_max
    print "Mean homozygous-reference SNP-locus count:", refhom_sum / samples
    print ""

    print "Heterozygous genotype percentage range:", hetpct_min, "-", hetpct_max
    print "Mean heterozygous genotype percentage:", hetpct_sum / samples
    print ""

    print "Homozygous-alt genotype percentage range:", homaltpct_min, "-", homaltpct_max
    print "Mean homozygous-alt genotype percentage:", homaltpct_sum / samples
}' "${COUNTS}" \
> "${SUMMARY}"


# ============================================================
# 10. WRITE QC REPORT
# ============================================================

{
    echo "Ascaris SNP genotype-count QC"
    echo "Date: $(date)"
    echo
    echo "Expected samples: ${EXPECTED_SAMPLES}"
    echo "PSC records: ${N_PSC}"
    echo "Per-sample TSV rows: ${N_SAMPLES}"
    echo
    echo "Input files:"
    echo "${STATS}"
    echo "${COUNTS}"
    echo
    echo "QC status: PASS"
    echo
    echo "Checks completed:"
    echo "- bcftools stats output exists"
    echo "- per-sample TSV exists"
    echo "- 146 PSC records detected"
    echo "- 146 per-sample records detected"
    echo "- no duplicate sample IDs"
    echo "- expected eight-column structure"
    echo "- n_called calculations validated"
    echo "- heterozygous percentages validated"
    echo "- homozygous-alt percentages validated"
} > "${QC_REPORT}"


# ============================================================
# 11. PRINT RESULTS TO SLURM OUTPUT
# ============================================================

echo "============================================================"
echo "QC COMPLETE"
echo "============================================================"
echo

cat "${QC_REPORT}"

echo
echo "============================================================"
echo "GENOTYPE SUMMARY"
echo "============================================================"
echo

cat "${SUMMARY}"

echo
echo "============================================================"
echo "FIRST 10 SAMPLES"
echo "============================================================"
echo

head -n 11 "${COUNTS}"

echo
echo "============================================================"
echo "FINISHED"
echo "============================================================"
echo "Finished: $(date)"
echo
echo "QC report:"
echo "${QC_REPORT}"
echo
echo "Summary:"
echo "${SUMMARY}"
echo
echo "Per-sample table:"
echo "${COUNTS}"
```

### 04E. CN heatmap and supplementary genotype-composition plots (R)

**Purpose.** Produce the heatmap from the chromosome-level CN table and the supplementary per-sample genotype-composition plots by country and host. The source R file contains both figure workflows in sequence, with detailed input validation and local output handling.

**Status: USER-SUPPLIED ORIGINAL SCRIPT.** The input metadata file was manually curated. File paths reflect the original RStudio environment; figure outputs were written locally when the `S:` drive was not writable.

```r
# ============================================================
# Ascaris chromosome-level relative copy-number heatmap
# ============================================================
#
# Input:
#   all_samples.depth5x.chromosome_CN.tsv
#   autosome_order.txt
#   sample_QC_list.csv
#
# Analysis:
#   - >=5x sample dataset
#   - 19 autosomes
#   - chromosome-level median relative copy number
#   - expected diploid relative CN = 2
#
# Heatmap annotations:
#   - Country
#   - Host
#
# NOTE:
#   Input files are read from the S: drive.
#   Figures are written locally because the current R session
#   can read from S: but does not have write permission to it.
# ============================================================


# ============================================================
# 1. LOAD PACKAGES
# ============================================================

library(tidyverse)
library(pheatmap)


# ============================================================
# 2. DEFINE INPUT DIRECTORY
# ============================================================

INPUT_DIR <- "S:/Research/Lauren_PopGen/Ploidy_CNV"

CN_FILE <- file.path(
  INPUT_DIR,
  "all_samples.depth5x.chromosome_CN.tsv"
)

CHR_FILE <- file.path(
  INPUT_DIR,
  "autosome_order.txt"
)

META_FILE <- file.path(
  INPUT_DIR,
  "sample_QC_list.csv"
)


# ============================================================
# 3. CHECK INPUT FILES EXIST
# ============================================================

input_files <- c(
  CN_FILE,
  CHR_FILE,
  META_FILE
)

if (!all(file.exists(input_files))) {
  
  missing_files <- input_files[
    !file.exists(input_files)
  ]
  
  stop(
    paste(
      "Missing input file(s):",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}

cat("All input files found.\n\n")


# ============================================================
# 4. READ INPUT FILES
# ============================================================

cn <- read.delim(
  CN_FILE,
  header = TRUE,
  stringsAsFactors = FALSE
)

chr_order <- readLines(
  CHR_FILE
)

meta <- read.csv(
  META_FILE,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# 5. INITIAL INPUT QC
# ============================================================

cat("============================================================\n")
cat("INPUT QC\n")
cat("============================================================\n\n")

cat("CN table dimensions:\n")
print(dim(cn))

cat("\nNumber of autosomes in chromosome-order file:\n")
print(length(chr_order))

cat("\nMetadata dimensions:\n")
print(dim(meta))

cat("\nMetadata columns:\n")
print(names(meta))

cat("\nFirst CN rows:\n")
print(head(cn))

cat("\nFirst metadata rows:\n")
print(head(meta))


# ============================================================
# 6. CLEAN METADATA
#
# Removes:
#   - leading/trailing whitespace
#   - repeated internal whitespace
#   - inconsistent capitalisation
#   - blank country/region/host values
#
# This fixes duplicated categories caused by hidden whitespace.
# ============================================================

meta_clean <- meta %>%
  mutate(
    
    # Sample ID
    sample_id = str_squish(sample_id),
    
    # Country
    country = str_squish(country),
    
    country = ifelse(
      is.na(country) | country == "",
      "Unknown",
      str_to_title(country)
    ),
    
    # Region
    # Retained in metadata but not used in final heatmap
    region = str_squish(region),
    
    region = ifelse(
      is.na(region) | region == "",
      "Unknown",
      str_to_title(region)
    ),
    
    # Host
    host = str_squish(host),
    
    host = ifelse(
      is.na(host) | host == "",
      "Unknown",
      str_to_title(host)
    )
  )


# ============================================================
# 7. CHECK CLEANED METADATA
# ============================================================

cat("\n============================================================\n")
cat("CLEANED METADATA\n")
cat("============================================================\n\n")

cat("Host categories:\n")
print(
  sort(
    unique(meta_clean$host)
  )
)

cat("\nRegion categories:\n")
print(
  sort(
    unique(meta_clean$region)
  )
)

cat("\nCountry categories:\n")
print(
  sort(
    unique(meta_clean$country)
  )
)

cat("\nHost counts:\n")
print(
  table(
    meta_clean$host,
    useNA = "ifany"
  )
)

cat("\nRegion counts:\n")
print(
  table(
    meta_clean$region,
    useNA = "ifany"
  )
)

cat("\nCountry counts:\n")
print(
  table(
    meta_clean$country,
    useNA = "ifany"
  )
)


# ============================================================
# 8. CREATE SAMPLE x CHROMOSOME CN MATRIX
# ============================================================

cn_matrix <- cn %>%
  select(
    sample,
    chrom,
    median_relative_CN
  ) %>%
  pivot_wider(
    names_from = chrom,
    values_from = median_relative_CN
  ) %>%
  column_to_rownames(
    "sample"
  )


# ============================================================
# 9. PUT CHROMOSOMES IN GENOMIC ORDER
# ============================================================

missing_chr <- setdiff(
  chr_order,
  colnames(cn_matrix)
)

if (length(missing_chr) > 0) {
  
  stop(
    paste(
      "Chromosome(s) from autosome_order.txt are missing:",
      paste(
        missing_chr,
        collapse = ", "
      )
    )
  )
}

cn_matrix <- cn_matrix[
  ,
  chr_order,
  drop = FALSE
]


# ============================================================
# 10. CHECK CN MATRIX
# ============================================================

cat("\n============================================================\n")
cat("CN MATRIX QC\n")
cat("============================================================\n\n")

cat("Matrix dimensions:\n")
print(dim(cn_matrix))

cat("\nCN range:\n")
print(
  range(
    cn_matrix,
    na.rm = TRUE
  )
)

cat("\nNumber of missing CN values:\n")
print(
  sum(
    is.na(cn_matrix)
  )
)


# ============================================================
# 11. CHECK SAMPLE IDs AGAINST METADATA
# ============================================================

missing_metadata <- setdiff(
  rownames(cn_matrix),
  meta_clean$sample_id
)

not_in_cn <- setdiff(
  meta_clean$sample_id,
  rownames(cn_matrix)
)

cat("\nSamples in CN matrix but absent from metadata:\n")
print(missing_metadata)

cat("\nMetadata samples absent from CN matrix:\n")
print(not_in_cn)

if (length(missing_metadata) > 0) {
  
  stop(
    "At least one CN sample has no corresponding metadata row."
  )
}


# ============================================================
# 12. CREATE HEATMAP ANNOTATION
#
# Only country and host are used.
#
# Region is omitted because country already provides the
# geographic information needed for the figure.
# ============================================================

annotation <- meta_clean %>%
  filter(
    sample_id %in% rownames(cn_matrix)
  ) %>%
  select(
    sample_id,
    country,
    host
  ) %>%
  column_to_rownames(
    "sample_id"
  )


# ============================================================
# 13. CHECK ANNOTATION
# ============================================================

cat("\n============================================================\n")
cat("ANNOTATION QC\n")
cat("============================================================\n\n")

cat("Annotation dimensions:\n")
print(dim(annotation))

cat("\nMissing annotation rows:\n")
print(
  setdiff(
    rownames(cn_matrix),
    rownames(annotation)
  )
)

cat("\nRemaining blank metadata fields:\n")
print(
  colSums(
    annotation == ""
  )
)

cat("\nERR068355 metadata:\n")

if ("ERR068355" %in% rownames(annotation)) {
  
  print(
    annotation[
      "ERR068355",
      ,
      drop = FALSE
    ]
  )
}


# ============================================================
# 14. ORDER SAMPLES
#
# Samples are ordered by:
#
#   Country
#     -> Host
#       -> Sample ID
#
# Unknown country is placed at the end.
#
# Samples are NOT clustered according to CN because the figure
# is intended to allow comparison with host/geography.
# ============================================================

sample_order <- meta_clean %>%
  filter(
    sample_id %in% rownames(cn_matrix)
  ) %>%
  arrange(
    country == "Unknown",
    country,
    host,
    sample_id
  ) %>%
  pull(
    sample_id
  )


# ============================================================
# 15. REORDER MATRIX AND ANNOTATION
# ============================================================

cn_matrix_ordered <- cn_matrix[
  sample_order,
  ,
  drop = FALSE
]

annotation_ordered <- annotation[
  sample_order,
  ,
  drop = FALSE
]


rows_match <- identical(
  rownames(cn_matrix_ordered),
  rownames(annotation_ordered)
)

cat("\nMatrix and annotation row order identical:\n")
print(rows_match)

if (!rows_match) {
  
  stop(
    "CN matrix and metadata annotation are not aligned."
  )
}


# ============================================================
# 16. RENAME AUTOSOMES
#
# Convert:
#
#   CM024168.1 ... CM024186.1
#
# to:
#
#   Chr1 ... Chr19
# ============================================================

colnames(cn_matrix_ordered) <- paste0(
  "Chr",
  seq_len(
    ncol(cn_matrix_ordered)
  )
)


# ============================================================
# 17. RENAME ANNOTATION HEADINGS
# ============================================================

colnames(annotation_ordered) <- c(
  "Country",
  "Host"
)


# ============================================================
# 18. CHECK FINAL SAMPLE ORDER
#
# No large white gaps are added between countries.
# Country grouping is instead shown using the annotation bar.
# ============================================================

country_summary <- annotation_ordered %>%
  as.data.frame() %>%
  count(
    Country,
    Host,
    name = "n"
  )

cat("\n============================================================\n")
cat("FINAL SAMPLE ORDER\n")
cat("============================================================\n\n")

print(country_summary)


# ============================================================
# 19. DEFINE COPY-NUMBER HEATMAP SCALE
#
# Relative CN = 2 represents expected diploid dosage.
#
# Blue  = lower relative copy number
# White = approximately diploid
# Red   = higher relative copy number
#
# Fixed scale:
#
# 1.5 ----------- 2.0 ----------- 2.5
#
# ============================================================

breaks <- seq(
  1.5,
  2.5,
  length.out = 101
)

heat_cols <- colorRampPalette(
  c(
    "#2166AC",
    "#F7F7F7",
    "#B2182B"
  )
)(100)


# ============================================================
# 20. DEFINE ANNOTATION COLOURS
#
# HOST:
#   Human   = black
#   Pig     = dark grey
#   Unknown = light grey
#
# COUNTRY:
#   brighter categorical colours for clear separation
# ============================================================


# ------------------------------------------------------------
# Host colours
# ------------------------------------------------------------

host_colours <- c(
  
  "Human" = "#000000",
  
  "Pig" = "#666666",
  
  "Unknown" = "#D3D3D3"
  
)


# ------------------------------------------------------------
# Country colours
# ------------------------------------------------------------

country_colours <- c(
  
  "Ethiopia" = "#E64B35",
  
  "Germany" = "#4DBBD5",
  
  "Hungary" = "#EFC000",
  
  "India" = "#DC7EC0",
  
  "Ireland" = "#00A087",
  
  "Kenya" = "#3C5488",
  
  "Philippines" = "#00BFC4",
  
  "Thailand" = "#F39B7F",
  
  "Uganda" = "#8491B4",
  
  "United Kingdom" = "#7E6148",
  
  "Unknown" = "#BDBDBD"
  
)


# ------------------------------------------------------------
# Combine annotation colours
# ------------------------------------------------------------

annotation_colours <- list(
  
  Country = country_colours,
  
  Host = host_colours
  
)


# ============================================================
# 21. CHECK ALL ANNOTATION CATEGORIES HAVE COLOURS
# ============================================================

missing_country_colours <- setdiff(
  unique(annotation_ordered$Country),
  names(country_colours)
)

missing_host_colours <- setdiff(
  unique(annotation_ordered$Host),
  names(host_colours)
)


if (length(missing_country_colours) > 0) {
  
  stop(
    paste(
      "Missing colour for country:",
      paste(
        missing_country_colours,
        collapse = ", "
      )
    )
  )
}


if (length(missing_host_colours) > 0) {
  
  stop(
    paste(
      "Missing colour for host:",
      paste(
        missing_host_colours,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 22. DEFINE LOCAL OUTPUT DIRECTORY
#
# The current R session can read from S: but cannot write
# there, so figures are saved locally.
# ============================================================

USER_DIR <- Sys.getenv(
  "USERPROFILE"
)

if (
  is.na(USER_DIR) ||
  USER_DIR == "" ||
  !dir.exists(USER_DIR)
) {
  
  stop(
    "Could not identify the Windows user-profile directory."
  )
}


OUTDIR <- file.path(
  USER_DIR,
  "CN_heatmap_outputs"
)


if (!dir.exists(OUTDIR)) {
  
  created <- dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = TRUE
  )
  
  if (
    !created &&
    !dir.exists(OUTDIR)
  ) {
    
    stop(
      paste(
        "Could not create output directory:",
        OUTDIR
      )
    )
  }
}


# ============================================================
# 23. TEST LOCAL WRITE PERMISSION
# ============================================================

write_status <- file.access(
  OUTDIR,
  mode = 2
)

cat("\n============================================================\n")
cat("OUTPUT DIRECTORY\n")
cat("============================================================\n\n")

cat(
  "Output directory:\n",
  OUTDIR,
  "\n\n"
)

cat(
  "Write-access status (0 = writable):\n"
)

print(
  write_status
)

if (write_status != 0) {
  
  stop(
    paste(
      "R cannot write to:",
      OUTDIR
    )
  )
}


# ============================================================
# 24. CREATE HEATMAP OBJECT
# ============================================================

heatmap_plot <- pheatmap(
  
  cn_matrix_ordered,
  
  # Keep samples in biological order
  cluster_rows = FALSE,
  
  # Keep chromosomes in genomic order
  cluster_cols = FALSE,
  
  # Values are already normalised relative CN
  scale = "none",
  
  # Relative CN colour scale
  color = heat_cols,
  breaks = breaks,
  
  # Country + Host annotations
  annotation_row = annotation_ordered,
  
  # Custom annotation colours
  annotation_colors = annotation_colours,
  
  # Individual sample IDs are not displayed
  show_rownames = FALSE,
  
  # Remove borders between cells
  border_color = NA,
  
  # Chromosome labels
  fontsize_col = 9,
  
  angle_col = 45,
  
  # Figure title
  main = "Autosomal relative copy number",
  
  # Create object without immediately drawing it
  silent = TRUE
)


# ============================================================
# 25. DISPLAY HEATMAP IF RUNNING INTERACTIVELY
# ============================================================

if (interactive()) {
  
  grid::grid.newpage()
  
  grid::grid.draw(
    heatmap_plot$gtable
  )
  
} else {
  
  cat(
    "\nR is running non-interactively, so no on-screen plot ",
    "window will be displayed.\n",
    "The PDF and PNG will still be written to disk.\n"
  )
}


# ============================================================
# 26. DEFINE OUTPUT FILES
# ============================================================

PDF_FILE <- file.path(
  OUTDIR,
  "Ascaris_autosomal_relative_CN_heatmap.pdf"
)

PNG_FILE <- file.path(
  OUTDIR,
  "Ascaris_autosomal_relative_CN_heatmap.png"
)


# ============================================================
# 27. SAVE PDF
# ============================================================

cat("\nWriting PDF:\n")
cat(PDF_FILE, "\n")

pheatmap(
  
  cn_matrix_ordered,
  
  cluster_rows = FALSE,
  
  cluster_cols = FALSE,
  
  scale = "none",
  
  color = heat_cols,
  
  breaks = breaks,
  
  annotation_row = annotation_ordered,
  
  annotation_colors = annotation_colours,
  
  show_rownames = FALSE,
  
  border_color = NA,
  
  fontsize_col = 9,
  
  angle_col = 45,
  
  main = "Autosomal relative copy number",
  
  filename = PDF_FILE,
  
  width = 9,
  
  height = 10
)


# ============================================================
# 28. SAVE PNG
# ============================================================

cat("\nWriting PNG:\n")
cat(PNG_FILE, "\n")

pheatmap(
  
  cn_matrix_ordered,
  
  cluster_rows = FALSE,
  
  cluster_cols = FALSE,
  
  scale = "none",
  
  color = heat_cols,
  
  breaks = breaks,
  
  annotation_row = annotation_ordered,
  
  annotation_colors = annotation_colours,
  
  show_rownames = FALSE,
  
  border_color = NA,
  
  fontsize_col = 9,
  
  angle_col = 45,
  
  main = "Autosomal relative copy number",
  
  filename = PNG_FILE,
  
  width = 9,
  
  height = 10
)


# ============================================================
# 29. CONFIRM FIGURES WERE CREATED
# ============================================================

PDF_EXISTS <- file.exists(
  PDF_FILE
)

PNG_EXISTS <- file.exists(
  PNG_FILE
)


cat("\n============================================================\n")
cat("FIGURE OUTPUT CHECK\n")
cat("============================================================\n\n")


cat(
  "PDF created:",
  PDF_EXISTS,
  "\n"
)

cat(
  "PNG created:",
  PNG_EXISTS,
  "\n"
)


cat(
  "\nPDF:\n",
  PDF_FILE,
  "\n"
)

cat(
  "\nPNG:\n",
  PNG_FILE,
  "\n"
)


if (!PDF_EXISTS) {
  
  warning(
    "PDF was not successfully created."
  )
}

if (!PNG_EXISTS) {
  
  warning(
    "PNG was not successfully created."
  )
}


# ============================================================
# 30. FINAL ANALYSIS QC
# ============================================================

cat("\n============================================================\n")
cat("FINAL HEATMAP SUMMARY\n")
cat("============================================================\n\n")


cat(
  "Samples plotted:",
  nrow(cn_matrix_ordered),
  "\n"
)

cat(
  "Autosomes plotted:",
  ncol(cn_matrix_ordered),
  "\n"
)

cat(
  "Minimum chromosome-level relative CN:",
  min(
    cn_matrix_ordered,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Maximum chromosome-level relative CN:",
  max(
    cn_matrix_ordered,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Median chromosome-level relative CN:",
  median(
    as.matrix(cn_matrix_ordered),
    na.rm = TRUE
  ),
  "\n"
)


cat("\nCountry counts in plotted dataset:\n")

print(
  table(
    annotation_ordered$Country
  )
)


cat("\nHost counts in plotted dataset:\n")

print(
  table(
    annotation_ordered$Host
  )
)


cat("\nFiles present in output directory:\n")

print(
  list.files(
    OUTDIR,
    full.names = TRUE
  )
)


cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")


cat(
  "\nFigures are saved locally in:\n",
  OUTDIR,
  "\n"
)


cat(
  "\nThey can now be copied manually to:\n",
  INPUT_DIR,
  "\n"
)


# ============================================================
# Ascaris per-sample genotype composition plots
# ============================================================

library(tidyverse)

# Check that patchwork is installed
if (!requireNamespace("patchwork", quietly = TRUE)) {
  
  stop(
    paste(
      "The 'patchwork' package is required.",
      "Install it using:",
      "install.packages('patchwork')"
    )
  )
}


# ============================================================
# 1. DEFINE INPUT DIRECTORY AND FILES
# ============================================================

INPUT_DIR <- "S:/Research/Lauren_PopGen/Ploidy_CNV"

GENO_FILE <- file.path(
  INPUT_DIR,
  "Ascaris.146.per_sample_genotype_counts.tsv"
)

META_FILE <- file.path(
  INPUT_DIR,
  "sample_QC_list.csv"
)


# ============================================================
# 2. CHECK INPUT FILES EXIST
# ============================================================

input_files <- c(
  GENO_FILE,
  META_FILE
)

if (!all(file.exists(input_files))) {
  
  missing_files <- input_files[
    !file.exists(input_files)
  ]
  
  stop(
    paste(
      "Missing input file(s):",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}


cat("============================================================\n")
cat("GENOTYPE COMPOSITION ANALYSIS\n")
cat("============================================================\n\n")

cat("All input files found.\n\n")


# ============================================================
# 3. READ DATA
# ============================================================

geno <- read.delim(
  GENO_FILE,
  header = TRUE,
  stringsAsFactors = FALSE
)

meta <- read.csv(
  META_FILE,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# 4. INITIAL INPUT QC
# ============================================================

cat("============================================================\n")
cat("INPUT QC\n")
cat("============================================================\n\n")

cat("Genotype table dimensions:\n")
print(dim(geno))

cat("\nMetadata dimensions:\n")
print(dim(meta))

cat("\nGenotype columns:\n")
print(names(geno))

cat("\nMetadata columns:\n")
print(names(meta))

cat("\nFirst genotype rows:\n")
print(head(geno))

cat("\nFirst metadata rows:\n")
print(head(meta))


# ============================================================
# 5. CHECK REQUIRED COLUMNS
# ============================================================

required_geno_columns <- c(
  "sample",
  "n_ref_hom",
  "n_hom_alt",
  "n_het",
  "n_missing",
  "n_called",
  "het_percent_called",
  "hom_alt_percent_called"
)

required_meta_columns <- c(
  "sample_id",
  "country",
  "region",
  "host"
)


missing_geno_columns <- setdiff(
  required_geno_columns,
  names(geno)
)

missing_meta_columns <- setdiff(
  required_meta_columns,
  names(meta)
)


if (length(missing_geno_columns) > 0) {
  
  stop(
    paste(
      "Missing genotype column(s):",
      paste(
        missing_geno_columns,
        collapse = ", "
      )
    )
  )
}


if (length(missing_meta_columns) > 0) {
  
  stop(
    paste(
      "Missing metadata column(s):",
      paste(
        missing_meta_columns,
        collapse = ", "
      )
    )
  )
}


cat("\nPASS: required columns found.\n")


# ============================================================
# 6. CLEAN METADATA
# ============================================================

meta_clean <- meta %>%
  mutate(
    
    # Sample ID
    sample_id = str_squish(sample_id),
    
    # Country
    country = str_squish(country),
    
    country = ifelse(
      is.na(country) | country == "",
      "Unknown",
      str_to_title(country)
    ),
    
    # Region
    region = str_squish(region),
    
    region = ifelse(
      is.na(region) | region == "",
      "Unknown",
      str_to_title(region)
    ),
    
    # Host
    host = str_squish(host),
    
    host = ifelse(
      is.na(host) | host == "",
      "Unknown",
      str_to_title(host)
    )
  )


# ============================================================
# 7. CLEAN GENOTYPE SAMPLE IDS
# ============================================================

geno <- geno %>%
  mutate(
    sample = str_squish(sample)
  )


# ============================================================
# 8. CHECK SAMPLE IDS
# ============================================================

cat("\n============================================================\n")
cat("SAMPLE-ID QC\n")
cat("============================================================\n\n")


duplicate_geno <- geno %>%
  count(sample) %>%
  filter(n > 1)

duplicate_meta <- meta_clean %>%
  count(sample_id) %>%
  filter(n > 1)


if (nrow(duplicate_geno) > 0) {
  
  cat("Duplicate genotype sample IDs:\n")
  print(duplicate_geno)
  
  stop(
    "Duplicate sample IDs detected in genotype table."
  )
}


if (nrow(duplicate_meta) > 0) {
  
  cat("Duplicate metadata sample IDs:\n")
  print(duplicate_meta)
  
  stop(
    "Duplicate sample IDs detected in metadata."
  )
}


cat("PASS: no duplicate sample IDs.\n")


# ============================================================
# 9. CHECK GENOTYPE SAMPLES AGAINST METADATA
# ============================================================

missing_metadata <- setdiff(
  geno$sample,
  meta_clean$sample_id
)

cat("\nGenotype samples absent from metadata:\n")
print(missing_metadata)


if (length(missing_metadata) > 0) {
  
  stop(
    "At least one genotype sample has no matching metadata row."
  )
}


cat("\nPASS: all genotype samples have metadata.\n")


# ============================================================
# 10. JOIN GENOTYPE DATA TO METADATA
# ============================================================

plot_data <- geno %>%
  left_join(
    meta_clean,
    by = c(
      "sample" = "sample_id"
    )
  )


# ============================================================
# 11. QC JOIN
# ============================================================

cat("\n============================================================\n")
cat("JOIN QC\n")
cat("============================================================\n\n")

cat(
  "Samples in genotype table:",
  nrow(geno),
  "\n"
)

cat(
  "Samples after metadata join:",
  nrow(plot_data),
  "\n"
)


if (nrow(plot_data) != nrow(geno)) {
  
  stop(
    "Number of rows changed after metadata join."
  )
}


cat("\nCountries:\n")

print(
  table(
    plot_data$country,
    useNA = "ifany"
  )
)


cat("\nHosts:\n")

print(
  table(
    plot_data$host,
    useNA = "ifany"
  )
)


cat("\nHeterozygous genotype percentage range:\n")

print(
  range(
    plot_data$het_percent_called,
    na.rm = TRUE
  )
)


cat("\nHomozygous-alt genotype percentage range:\n")

print(
  range(
    plot_data$hom_alt_percent_called,
    na.rm = TRUE
  )
)


# ============================================================
# 12. SET COUNTRY ORDER
# ============================================================

country_order <- c(
  "Ethiopia",
  "Kenya",
  "Uganda",
  "India",
  "Philippines",
  "Thailand",
  "Germany",
  "Hungary",
  "Ireland",
  "United Kingdom",
  "Unknown"
)


# Check for countries not already included in predefined order

extra_countries <- setdiff(
  unique(plot_data$country),
  country_order
)


if (length(extra_countries) > 0) {
  
  cat(
    "\nAdditional countries found and appended to plotting order:\n"
  )
  
  print(extra_countries)
  
  country_order <- c(
    country_order[
      country_order != "Unknown"
    ],
    sort(extra_countries),
    "Unknown"
  )
}


plot_data <- plot_data %>%
  mutate(
    
    country = factor(
      country,
      levels = country_order
    ),
    
    host = factor(
      host,
      levels = c(
        "Human",
        "Pig",
        "Unknown"
      )
    )
  )


# ============================================================
# 13. DEFINE HOST COLOURS AND SHAPES
# ============================================================

host_colours <- c(
  
  "Human" = "#000000",
  
  "Pig" = "#777777",
  
  "Unknown" = "#CFCFCF"
  
)


host_shapes <- c(
  
  "Human" = 16,
  
  "Pig" = 17,
  
  "Unknown" = 15
  
)


# ============================================================
# 14. CREATE COUNTRY PLOT
#
# All country categories retain a boxplot.
#
# For categories represented by one sample, the boxplot
# naturally collapses to a horizontal line. Individual points
# remain visible so sample size is clear.
# ============================================================

p_country <- ggplot(
  
  plot_data,
  
  aes(
    x = country,
    y = het_percent_called
  )
  
) +
  
  geom_boxplot(
    
    width = 0.65,
    
    outlier.shape = NA,
    
    fill = "white",
    
    colour = "black",
    
    linewidth = 0.5
    
  ) +
  
  geom_jitter(
    
    aes(
      colour = host,
      shape = host
    ),
    
    width = 0.16,
    
    height = 0,
    
    size = 2,
    
    alpha = 0.75
    
  ) +
  
  scale_colour_manual(
    
    values = host_colours,
    
    name = "Host",
    
    drop = TRUE
    
  ) +
  
  scale_shape_manual(
    
    values = host_shapes,
    
    name = "Host",
    
    drop = TRUE
    
  ) +
  
  labs(
    
    x = "Country",
    
    y = "Heterozygous genotypes at called polymorphic loci (%)"
    
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    legend.position = "right"
    
  )


# ============================================================
# 15. CREATE HOST PLOT
# ============================================================

p_host <- ggplot(
  
  plot_data,
  
  aes(
    x = host,
    y = het_percent_called
  )
  
) +
  
  geom_boxplot(
    
    width = 0.5,
    
    outlier.shape = NA,
    
    fill = "white",
    
    colour = "black",
    
    linewidth = 0.5
    
  ) +
  
  geom_jitter(
    
    aes(
      colour = host
    ),
    
    width = 0.12,
    
    height = 0,
    
    size = 2,
    
    alpha = 0.70
    
  ) +
  
  scale_colour_manual(
    
    values = host_colours,
    
    guide = "none",
    
    drop = TRUE
    
  ) +
  
  labs(
    
    x = "Host",
    
    y = "Heterozygous genotypes at called polymorphic loci (%)"
    
  ) +
  
  theme_classic(
    base_size = 12
  )


# ============================================================
# 16. DISPLAY INDIVIDUAL PLOTS
#
# print() is used explicitly so plots display when the script
# is sourced in RStudio.
# ============================================================

cat("\nDisplaying country plot...\n")

print(
  p_country
)


cat("\nDisplaying host plot...\n")

print(
  p_host
)


# ============================================================
# 17. SUMMARY STATISTICS BY COUNTRY
# ============================================================

country_summary <- plot_data %>%
  
  group_by(
    country
  ) %>%
  
  summarise(
    
    n = n(),
    
    mean_het_percent = mean(
      het_percent_called,
      na.rm = TRUE
    ),
    
    median_het_percent = median(
      het_percent_called,
      na.rm = TRUE
    ),
    
    min_het_percent = min(
      het_percent_called,
      na.rm = TRUE
    ),
    
    max_het_percent = max(
      het_percent_called,
      na.rm = TRUE
    ),
    
    mean_hom_alt_percent = mean(
      hom_alt_percent_called,
      na.rm = TRUE
    ),
    
    median_hom_alt_percent = median(
      hom_alt_percent_called,
      na.rm = TRUE
    ),
    
    .groups = "drop"
    
  )


cat("\n============================================================\n")
cat("COUNTRY SUMMARY\n")
cat("============================================================\n\n")

print(
  country_summary
)


# ============================================================
# 18. SUMMARY STATISTICS BY HOST
# ============================================================

host_summary <- plot_data %>%
  
  group_by(
    host
  ) %>%
  
  summarise(
    
    n = n(),
    
    mean_het_percent = mean(
      het_percent_called,
      na.rm = TRUE
    ),
    
    median_het_percent = median(
      het_percent_called,
      na.rm = TRUE
    ),
    
    min_het_percent = min(
      het_percent_called,
      na.rm = TRUE
    ),
    
    max_het_percent = max(
      het_percent_called,
      na.rm = TRUE
    ),
    
    mean_hom_alt_percent = mean(
      hom_alt_percent_called,
      na.rm = TRUE
    ),
    
    median_hom_alt_percent = median(
      hom_alt_percent_called,
      na.rm = TRUE
    ),
    
    .groups = "drop"
    
  )


cat("\n============================================================\n")
cat("HOST SUMMARY\n")
cat("============================================================\n\n")

print(
  host_summary
)


# ============================================================
# 19. GIVE BOTH PLOTS THE SAME Y-AXIS RANGE
#
# Observed overall heterozygous percentage range is
# approximately 5.95 - 17.53%, so both panels use 5.5 - 18%.
# ============================================================

p_country_final <- p_country +
  
  coord_cartesian(
    ylim = c(
      5.5,
      18
    )
  ) +
  
  theme(
    
    plot.tag = element_text(
      face = "bold",
      size = 14
    )
    
  )


p_host_final <- p_host +
  
  coord_cartesian(
    ylim = c(
      5.5,
      18
    )
  ) +
  
  theme(
    
    plot.tag = element_text(
      face = "bold",
      size = 14
    )
    
  )


# ============================================================
# 20. COMBINE COUNTRY + HOST PLOTS
#
# wrap_plots() is used instead of the "/" patchwork operator
# because the latter was not being recognised by this R
# session.
#
# A = country
# B = host
# ============================================================

combined_genotype_plot <- patchwork::wrap_plots(
  
  p_country_final,
  
  p_host_final,
  
  ncol = 1,
  
  heights = c(
    1.15,
    1
  )
  
) +
  
  patchwork::plot_annotation(
    tag_levels = "A"
  )


# ============================================================
# 21. DISPLAY COMBINED A/B FIGURE
# ============================================================

cat("\nDisplaying combined supplementary figure...\n")

print(
  combined_genotype_plot
)


# ============================================================
# 22. DEFINE LOCAL OUTPUT DIRECTORY
#
# S: can be read but is not writable from the current R
# session, so all figures are saved locally.
# ============================================================

USER_DIR <- Sys.getenv(
  "USERPROFILE"
)


if (
  is.na(USER_DIR) ||
  USER_DIR == "" ||
  !dir.exists(USER_DIR)
) {
  
  stop(
    "Could not identify the Windows user-profile directory."
  )
}


GENOTYPE_OUTDIR <- file.path(
  USER_DIR,
  "Genotype_composition_outputs"
)


if (!dir.exists(GENOTYPE_OUTDIR)) {
  
  created <- dir.create(
    
    GENOTYPE_OUTDIR,
    
    recursive = TRUE,
    
    showWarnings = TRUE
    
  )
  
  
  if (
    !created &&
    !dir.exists(GENOTYPE_OUTDIR)
  ) {
    
    stop(
      paste(
        "Could not create output directory:",
        GENOTYPE_OUTDIR
      )
    )
    
  }
}


# ============================================================
# 23. TEST WRITE PERMISSION
# ============================================================

write_status <- file.access(
  GENOTYPE_OUTDIR,
  mode = 2
)


cat("\n============================================================\n")
cat("OUTPUT DIRECTORY\n")
cat("============================================================\n\n")

cat(
  "Output directory:\n",
  GENOTYPE_OUTDIR,
  "\n\n"
)

cat(
  "Write-access status (0 = writable):\n"
)

print(
  write_status
)


if (write_status != 0) {
  
  stop(
    paste(
      "R cannot write to:",
      GENOTYPE_OUTDIR
    )
  )
}


# ============================================================
# 24. DEFINE OUTPUT FILENAMES
# ============================================================

COUNTRY_PDF <- file.path(
  GENOTYPE_OUTDIR,
  "Ascaris_heterozygous_genotypes_by_country.pdf"
)

COUNTRY_PNG <- file.path(
  GENOTYPE_OUTDIR,
  "Ascaris_heterozygous_genotypes_by_country.png"
)


HOST_PDF <- file.path(
  GENOTYPE_OUTDIR,
  "Ascaris_heterozygous_genotypes_by_host.pdf"
)

HOST_PNG <- file.path(
  GENOTYPE_OUTDIR,
  "Ascaris_heterozygous_genotypes_by_host.png"
)


COMBINED_PDF <- file.path(
  GENOTYPE_OUTDIR,
  "Ascaris_heterozygous_genotypes_combined.pdf"
)

COMBINED_PNG <- file.path(
  GENOTYPE_OUTDIR,
  "Ascaris_heterozygous_genotypes_combined.png"
)


COUNTRY_SUMMARY_FILE <- file.path(
  GENOTYPE_OUTDIR,
  "Ascaris_heterozygous_genotypes_country_summary.tsv"
)

HOST_SUMMARY_FILE <- file.path(
  GENOTYPE_OUTDIR,
  "Ascaris_heterozygous_genotypes_host_summary.tsv"
)


# ============================================================
# 25. SAVE COUNTRY PLOT
# ============================================================

cat("\nSaving country plot...\n")


ggsave(
  
  filename = COUNTRY_PDF,
  
  plot = p_country_final,
  
  width = 10,
  
  height = 5.5,
  
  units = "in"
  
)


ggsave(
  
  filename = COUNTRY_PNG,
  
  plot = p_country_final,
  
  width = 10,
  
  height = 5.5,
  
  units = "in",
  
  dpi = 600
  
)


# ============================================================
# 26. SAVE HOST PLOT
# ============================================================

cat("\nSaving host plot...\n")


ggsave(
  
  filename = HOST_PDF,
  
  plot = p_host_final,
  
  width = 7,
  
  height = 5.5,
  
  units = "in"
  
)


ggsave(
  
  filename = HOST_PNG,
  
  plot = p_host_final,
  
  width = 7,
  
  height = 5.5,
  
  units = "in",
  
  dpi = 600
  
)


# ============================================================
# 27. SAVE COMBINED A/B SUPPLEMENTARY FIGURE
# ============================================================

cat("\nSaving combined A/B figure...\n")


ggsave(
  
  filename = COMBINED_PDF,
  
  plot = combined_genotype_plot,
  
  width = 10,
  
  height = 10,
  
  units = "in"
  
)


ggsave(
  
  filename = COMBINED_PNG,
  
  plot = combined_genotype_plot,
  
  width = 10,
  
  height = 10,
  
  units = "in",
  
  dpi = 600
  
)


# ============================================================
# 28. SAVE SUMMARY TABLES
# ============================================================

write.table(
  
  country_summary,
  
  file = COUNTRY_SUMMARY_FILE,
  
  sep = "\t",
  
  quote = FALSE,
  
  row.names = FALSE
  
)


write.table(
  
  host_summary,
  
  file = HOST_SUMMARY_FILE,
  
  sep = "\t",
  
  quote = FALSE,
  
  row.names = FALSE
  
)


# ============================================================
# 29. CONFIRM OUTPUTS
# ============================================================

output_files <- c(
  
  COUNTRY_PDF,
  
  COUNTRY_PNG,
  
  HOST_PDF,
  
  HOST_PNG,
  
  COMBINED_PDF,
  
  COMBINED_PNG,
  
  COUNTRY_SUMMARY_FILE,
  
  HOST_SUMMARY_FILE
  
)


output_status <- file.exists(
  output_files
)


output_check <- data.frame(
  
  file = output_files,
  
  created = output_status,
  
  stringsAsFactors = FALSE
  
)


cat("\n============================================================\n")
cat("OUTPUT CHECK\n")
cat("============================================================\n\n")

print(
  output_check
)


if (!all(output_status)) {
  
  warning(
    "At least one genotype-composition output was not created."
  )
}


# ============================================================
# 30. FINAL ANALYSIS SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("FINAL GENOTYPE COMPOSITION SUMMARY\n")
cat("============================================================\n\n")


cat(
  "Samples analysed:",
  nrow(plot_data),
  "\n"
)


cat(
  "Minimum heterozygous genotype percentage:",
  min(
    plot_data$het_percent_called,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Maximum heterozygous genotype percentage:",
  max(
    plot_data$het_percent_called,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Mean heterozygous genotype percentage:",
  mean(
    plot_data$het_percent_called,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Minimum homozygous-alt genotype percentage:",
  min(
    plot_data$hom_alt_percent_called,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Maximum homozygous-alt genotype percentage:",
  max(
    plot_data$hom_alt_percent_called,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Mean homozygous-alt genotype percentage:",
  mean(
    plot_data$hom_alt_percent_called,
    na.rm = TRUE
  ),
  "\n"
)


cat("\nCountry sample counts:\n")

print(
  table(
    plot_data$country
  )
)


cat("\nHost sample counts:\n")

print(
  table(
    plot_data$host
  )
)


cat("\nFiles created:\n")

print(
  list.files(
    GENOTYPE_OUTDIR,
    full.names = TRUE
  )
)


cat("\n============================================================\n")
cat("GENOTYPE COMPOSITION ANALYSIS COMPLETE\n")
cat("============================================================\n")


cat(
  "\nFigures and summary tables are saved locally in:\n",
  GENOTYPE_OUTDIR,
  "\n"
)


cat(
  "\nThey can now be copied manually to:\n",
  INPUT_DIR,
  "\n"
)
```

## Outputs and QC checks

| Check | Expected/observed from source code |
|---|---|
| BAMs processed for ploidy | 156 `MarkDuplicates` BAMs |
| CN autosomes | 19 |
| CN 25-kb windows per sample | 7,111 expected |
| CN minimum mappability within windows | ≥90% |
| CN primary depth | ≥5× median autosomal depth |
| CN sensitivity depth | ≥10× |
| CN scaling baseline | Relative CN = 2 at sample autosomal median depth |
| Genotype count samples | 146 |
| Genotype-count fields | 8-column sample summary, `PSC` QC |

These are **pipeline parameters and validation conditions**, not proof that every original job succeeded; original SLURM logs and tables should be retained where available.

## Provenance and limitations

- The 25-kb BED (`Ascaris.autosomes19.25kb.unique90.bed`) exists as an analytical input; the exact generating `qc_25kb_mappability_070826.sub` script is not included in the current sources. Do not claim a documented executable BED-creation step until recovered.


- The chromosome CN heatmap reflects **relative sequencing depth** and is not a karyotype or formal confirmation of ploidy/aneuploidy.
- The copy-number dataset uses all 156 marked-duplicate BAMs before the final variant-QC exclusions; its ≥5× subset is distinct from the nuclear variant-calling n = 146 set.
- The mean-depth baseline is the median across included 25-kb windows, not necessarily coverage at every nucleotide.
- `sample_QC_list.csv` was **manually curated**; keep the actual file and describe its provenance.
- The original `SNP_genotype_counts` script is now embedded in full (04C); 04D independently checks its PSC record counts and percentage calculations. Retain the original SLURM log where available.
- **The PIXY all-sites correction is a separate pipeline**: `99_PIXY_allsites_correction.md`.
