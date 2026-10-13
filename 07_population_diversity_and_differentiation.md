# 07 — Windowed diversity and differentiation (historical thesis run)

## Purpose

Windowed diversity and differentiation (historical thesis run). This file archives the commands and analysis notes from the thesis/paper workflow.

## Inputs

Historical n=146 all-sites VCF, PIXY population lists. Source filenames and absolute HPC paths are retained in the original script blocks below.

## Workflow

See the original annotated analysis sequence below. Distinct script blocks are independent jobs/R sessions, **not one executable script**.

## Annotated scripts

> **Historical thesis pipeline — not the corrected PIXY input.** The all-sites VCF here was not restricted to uniquely mappable positions. The separately corrected VCF is documented in [`99`](99_PIXY_allsites_correction.md). 

### Historical all-sites joint-genotyping handoff 

This branch is distinct from the variant-only `GenotypeGVCFs` output in `05`. Earlier records describe a separate joint-genotyping job (`allsites_joint`) on the merged GVCFs with `--include-non-variant-sites`. 

```bash
REF=/parallel_scratch/lw01165/01_Reference/Ascaris_primary_scaffolds_plus_X.fasta
IN_GVCF=/parallel_scratch/lw01165/04_Variants/04a_merged/merged_all_samples.g.vcf.gz
OUT_VCF=/parallel_scratch/lw01165/04_Variants/04a_merged/merged_all_samples.ALLSITES.vcf.gz

gatk GenotypeGVCFs \
  -R "${REF}" \
  -V "${IN_GVCF}" \
  --include-non-variant-sites \
  -O "${OUT_VCF}"

gatk IndexFeatureFile -I "${OUT_VCF}"
```

The historical code below then subsets this output to the same 146 sample IDs and subsequently to 19 autosomes. The original all-sites branch **did not carry over the uniquely mappable-region restriction**; see `99` for the separate audit/correction.


## Comparisons between and within samples, stratified by host and geography

This annotated code section documents how the all-sites VCF was prepared and how windowed FST, dXY, π and Tajima’s D were generated. These statistics provide complementary views of differentiation and diversity:

- **FST**: relative differentiation between two groups.
- **dXY**: absolute nucleotide divergence between two groups.
- **π**: within-population nucleotide diversity.
- **Tajima’s D**: allele-frequency spectrum summary, useful for detecting shifts such as rare-allele excess or intermediate-frequency allele excess.

The key methodological point is that PIXY was run on an **all-sites VCF**, so invariant sites were retained for estimates such as dXY and π.


## Cluster code: 01 subset allsites to 146

#Purpose: prepare or subset the input VCF before PIXY.

```bash
#!/bin/bash
#SBATCH --job-name=allsites_146_subset
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=06:00:00
#SBATCH --output=allsites_146_subset.%j.out
#SBATCH --error=allsites_146_subset.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env

IN_ALLSITES=/parallel_scratch/lw01165/04_Variants/04a_merged/merged_all_samples.ALLSITES.vcf.gz
SNP146=/parallel_scratch/lw01165/05_analysis/05b_146_data/Ascaris.autosomes.SNPs.146.maxmiss0.90.vcf.gz
OUTDIR=/parallel_scratch/lw01165/05_analysis/05c_allsites_146

mkdir -p ${OUTDIR}

# extract sample list
bcftools query -l ${SNP146} > ${OUTDIR}/samples_146.txt

# subset all-sites VCF
bcftools view \
-S ${OUTDIR}/samples_146.txt \
-Oz \
-o ${OUTDIR}/Ascaris.autosomes.ALLSITES.146.vcf.gz \
${IN_ALLSITES}

# index
tabix -p vcf ${OUTDIR}/Ascaris.autosomes.ALLSITES.146.vcf.gz
```


## Cluster code: 02 restrict allsites to 19 autosomes

# Purpose: prepare or subset the input VCF before PIXY.

```bash
#!/bin/bash
#SBATCH --job-name=restrict_allsites_autosomes
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=08:00:00
#SBATCH --output=restrict_allsites_autosomes.%j.out
#SBATCH --error=restrict_allsites_autosomes.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate tools_env

WORKDIR=/parallel_scratch/lw01165/05_analysis/05c_allsites_146
cd ${WORKDIR}

OLDNAME=Ascaris.autosomes.ALLSITES.146.vcf.gz
OLDINDEX=${OLDNAME}.tbi

RENAMED=Ascaris.ALLSITES.146.vcf.gz
RENAMED_INDEX=${RENAMED}.tbi

AUTOSOME_VCF=Ascaris.autosomes19.ALLSITES.146.vcf.gz

AUTOSOME_LIST=${WORKDIR}/autosomes19.list

echo "===================================================="
echo "Renaming old file and restricting to 19 autosomes"
echo "Start: $(date)"
echo "Working directory: ${WORKDIR}"
echo "===================================================="

# check inputs exist
for f in "${OLDNAME}" "${OLDINDEX}"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: Missing required file: $f"
    exit 1
  fi
done

# rename original file so the name reflects content better
if [ ! -f "${RENAMED}" ]; then
  echo "Renaming ${OLDNAME} -> ${RENAMED}"
  mv "${OLDNAME}" "${RENAMED}"
  mv "${OLDINDEX}" "${RENAMED_INDEX}"
else
  echo "Renamed file already exists: ${RENAMED}"
  echo "Leaving original rename step unchanged"
fi

# write autosome contig list (first 19 chromosomes)
cat > "${AUTOSOME_LIST}" <<EOF
CM024168.1
CM024169.1
CM024170.1
CM024171.1
CM024172.1
CM024173.1
CM024174.1
CM024175.1
CM024176.1
CM024177.1
CM024178.1
CM024179.1
CM024180.1
CM024181.1
CM024182.1
CM024183.1
CM024184.1
CM024185.1
CM024186.1
EOF

echo "Autosome list written to ${AUTOSOME_LIST}"
cat "${AUTOSOME_LIST}"

# subset to autosomes only
echo "Subsetting all-sites VCF to the 19 autosomes..."
bcftools view \
  -r $(paste -sd, "${AUTOSOME_LIST}") \
  -Oz \
  -o "${AUTOSOME_VCF}" \
  "${RENAMED}"

# index autosome-only VCF
tabix -p vcf "${AUTOSOME_VCF}"

echo "Finished subsetting: $(date)"
echo "Output files:"
ls -lh "${RENAMED}" "${RENAMED_INDEX}" "${AUTOSOME_VCF}" "${AUTOSOME_VCF}.tbi"

echo "===================================================="
echo "Contigs present in autosome-only VCF:"
bcftools index -s "${AUTOSOME_VCF}"
echo "===================================================="
```


## Cluster code: 01 pixy fst pairwise 146

#Purpose: run or summarise PIXY windowed statistics.

```bash
#!/bin/bash
#SBATCH --job-name=pixy_fst_146
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=16
#SBATCH --mem=130G
#SBATCH --time=24:00:00
#SBATCH --output=pixy_fst_146.%j.out
#SBATCH --error=pixy_fst_146.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate pixy_env

VCF=/parallel_scratch/lw01165/05_analysis/05c_allsites_146/Ascaris.autosomes19.ALLSITES.146.vcf.gz
POPDIR=/parallel_scratch/lw01165/05_analysis/05i_pixy/popfiles
OUTDIR=/parallel_scratch/lw01165/05_analysis/05i_pixy/results

mkdir -p ${OUTDIR}

echo "===================================================="
echo "Running PIXY FST analyses"
echo "Start: $(date)"
echo "Host: $(hostname)"
echo "VCF: ${VCF}"
echo "Popfile directory: ${POPDIR}"
echo "Output directory: ${OUTDIR}"
echo "===================================================="

# check required files
for f in \
  "${VCF}" \
  "${VCF}.tbi" \
  "${POPDIR}/kenya_vs_ethiopia.txt" \
  "${POPDIR}/africa_vs_asia.txt" \
  "${POPDIR}/africa_vs_europe.txt" \
  "${POPDIR}/asia_vs_europe.txt" \
  "${POPDIR}/human_vs_pig.txt"
do
  if [ ! -f "$f" ]; then
    echo "ERROR: Missing required file: $f"
    exit 1
  fi
done

echo
echo "VCF sample count:"
bcftools query -l "${VCF}" | wc -l

echo
echo "Checking pixy is available:"
which pixy

for popfile in kenya_vs_ethiopia africa_vs_asia africa_vs_europe asia_vs_europe human_vs_pig
do
  echo
  echo "----------------------------------------------------"
  echo "Running PIXY for ${popfile}"
  echo "Start: $(date)"
  echo "----------------------------------------------------"

  pixy \
    --stats fst \
    --vcf "${VCF}" \
    --populations "${POPDIR}/${popfile}.txt" \
    --window_size 10000 \
    --n_cores "${SLURM_CPUS_PER_TASK}" \
    --output_folder "${OUTDIR}" \
    --output_prefix "${popfile}"

  echo "Finished ${popfile}: $(date)"
done

echo
echo "All PIXY runs complete: $(date)"
echo "Output files:"
ls -lh "${OUTDIR}"
echo "===================================================="
```


## Cluster code: 02 pixy dxy pairwise 146

# Purpose: run or summarise PIXY windowed statistics.

```bash
#!/bin/bash
#SBATCH --job-name=pixy_dxy_146
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=16
#SBATCH --mem=130G
#SBATCH --time=24:00:00
#SBATCH --output=pixy_dxy_146.%j.out
#SBATCH --error=pixy_dxy_146.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate pixy_env

VCF=/parallel_scratch/lw01165/05_analysis/05c_allsites_146/Ascaris.autosomes19.ALLSITES.146.vcf.gz
POPDIR=/parallel_scratch/lw01165/05_analysis/05i_pixy/popfiles
OUTDIR=/parallel_scratch/lw01165/05_analysis/05i_pixy/results_dxy

mkdir -p ${OUTDIR}

echo "===================================================="
echo "Running PIXY dXY analyses"
echo "Start: $(date)"
echo "Host: $(hostname)"
echo "VCF: ${VCF}"
echo "Popfile directory: ${POPDIR}"
echo "Output directory: ${OUTDIR}"
echo "===================================================="

# check required files
for f in \
  "${VCF}" \
  "${VCF}.tbi" \
  "${POPDIR}/kenya_vs_ethiopia.txt" \
  "${POPDIR}/africa_vs_asia.txt" \
  "${POPDIR}/africa_vs_europe.txt" \
  "${POPDIR}/asia_vs_europe.txt" \
  "${POPDIR}/africa_vs_nonafrica.txt" \
  "${POPDIR}/human_vs_pig.txt"
do
  if [ ! -f "$f" ]; then
    echo "ERROR: Missing required file: $f"
    exit 1
  fi
done

echo
echo "VCF sample count:"
bcftools query -l "${VCF}" | wc -l

echo
echo "Checking pixy is available:"
which pixy

for popfile in kenya_vs_ethiopia africa_vs_asia africa_vs_europe asia_vs_europe human_vs_pig
do
  echo
  echo "----------------------------------------------------"
  echo "Running PIXY for ${popfile}"
  echo "Start: $(date)"
  echo "----------------------------------------------------"

  pixy \
    --stats dxy \
    --vcf "${VCF}" \
    --populations "${POPDIR}/${popfile}.txt" \
    --window_size 10000 \
    --n_cores "${SLURM_CPUS_PER_TASK}" \
    --output_folder "${OUTDIR}" \
    --output_prefix "${popfile}"

  echo "Finished ${popfile}: $(date)"
done

echo
echo "All PIXY runs complete: $(date)"
echo "Output files:"
ls -lh "${OUTDIR}"
echo "===================================================="
```


## Cluster code: 03 pixy fst africa vs nonafrica

# Purpose: run or summarise PIXY windowed statistics.

```bash
#!/bin/bash
#SBATCH --job-name=pixy_fst_africa_nonafrica
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=16
#SBATCH --mem=130G
#SBATCH --time=24:00:00
#SBATCH --output=pixy_fst_africa_nonafrica.%j.out
#SBATCH --error=pixy_fst_africa_nonafrica.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate pixy_env

VCF=/parallel_scratch/lw01165/05_analysis/05c_allsites_146/Ascaris.autosomes19.ALLSITES.146.vcf.gz
POPFILE=/parallel_scratch/lw01165/05_analysis/05i_pixy/popfiles/africa_vs_nonafrica.txt
OUTDIR=/parallel_scratch/lw01165/05_analysis/05i_pixy/results

mkdir -p "${OUTDIR}"

echo "===================================================="
echo "Running PIXY FST: Africa vs non-Africa"
echo "Start: $(date)"
echo "Host: $(hostname)"
echo "VCF: ${VCF}"
echo "Population file: ${POPFILE}"
echo "Output directory: ${OUTDIR}"
echo "===================================================="

# check required files
for f in \
  "${VCF}" \
  "${VCF}.tbi" \
  "${POPFILE}"
do
  if [ ! -f "$f" ]; then
    echo "ERROR: Missing required file: $f"
    exit 1
  fi
done

echo
echo "VCF sample count:"
bcftools query -l "${VCF}" | wc -l

echo
echo "Preview of population file:"
head -5 "${POPFILE}"

echo
echo "Checking pixy is available:"
which pixy

echo
echo "----------------------------------------------------"
echo "Running PIXY for Africa vs non-Africa"
echo "Start: $(date)"
echo "----------------------------------------------------"

pixy \
  --stats fst \
  --vcf "${VCF}" \
  --populations "${POPFILE}" \
  --window_size 10000 \
  --n_cores "${SLURM_CPUS_PER_TASK}" \
  --output_folder "${OUTDIR}" \
  --output_prefix "africa_vs_nonafrica"

echo
echo "Finished: $(date)"
echo "Output files:"
ls -lh "${OUTDIR}" | grep africa_vs_nonafrica
echo "===================================================="
```


## Cluster code: 04 pixy dxy africa vs nonafrica

# Purpose: run or summarise PIXY windowed statistics.

```bash
#!/bin/bash
#SBATCH --job-name=pixy_dxy_africa_nonafrica
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=16
#SBATCH --mem=130G
#SBATCH --time=24:00:00
#SBATCH --output=pixy_dxy_africa_nonafrica.%j.out
#SBATCH --error=pixy_dxy_africa_nonafrica.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate pixy_env

VCF=/parallel_scratch/lw01165/05_analysis/05c_allsites_146/Ascaris.autosomes19.ALLSITES.146.vcf.gz
POPDIR=/parallel_scratch/lw01165/05_analysis/05i_pixy/popfiles
OUTDIR=/parallel_scratch/lw01165/05_analysis/05i_pixy/results_dxy_africa_nonafrica

mkdir -p "${OUTDIR}"

echo "===================================================="
echo "Running PIXY dXY analysis: Africa vs Non-Africa"
echo "Start: $(date)"
echo "Host: $(hostname)"
echo "VCF: ${VCF}"
echo "Popfile directory: ${POPDIR}"
echo "Output directory: ${OUTDIR}"
echo "===================================================="

for f in \
  "${VCF}" \
  "${VCF}.tbi" \
  "${POPDIR}/africa_vs_nonafrica.txt"
do
  if [ ! -f "$f" ]; then
    echo "ERROR: Missing required file: $f"
    exit 1
  fi
done

echo
echo "VCF sample count:"
bcftools query -l "${VCF}" | wc -l

echo
echo "Checking pixy is available:"
which pixy

echo
echo "Checking pop assignments:"
column -t "${POPDIR}/africa_vs_nonafrica.txt" | head

echo
echo "----------------------------------------------------"
echo "Running PIXY for africa_vs_nonafrica"
echo "Start: $(date)"
echo "----------------------------------------------------"

pixy \
  --stats dxy \
  --vcf "${VCF}" \
  --populations "${POPDIR}/africa_vs_nonafrica.txt" \
  --window_size 10000 \
  --n_cores "${SLURM_CPUS_PER_TASK}" \
  --output_folder "${OUTDIR}" \
  --output_prefix "africa_vs_nonafrica"

echo "Finished africa_vs_nonafrica: $(date)"

echo
echo "All PIXY runs complete: $(date)"
echo "Output files:"
ls -lh "${OUTDIR}"
echo "===================================================="
```


## Cluster code: 05 pixy diversity regions countries hosts 146

# Purpose: run or summarise PIXY windowed statistics.

```bash
#!/bin/bash
#SBATCH --job-name=pixy_diversity_146
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=pixy_diversity_146.%j.out
#SBATCH --error=pixy_diversity_146.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate pixy_env

VCF=/parallel_scratch/lw01165/05_analysis/05c_allsites_146/Ascaris.autosomes19.ALLSITES.146.vcf.gz
POPDIR=/parallel_scratch/lw01165/05_analysis/05i_pixy/popfiles_diversity
OUTDIR=/parallel_scratch/lw01165/05_analysis/05i_pixy/results_diversity

mkdir -p ${OUTDIR}

echo "===================================================="
echo "Running PIXY pi and Tajima's D analyses"
echo "Start: $(date)"
echo "Host: $(hostname)"
echo "VCF: ${VCF}"
echo "Popfile directory: ${POPDIR}"
echo "Output directory: ${OUTDIR}"
echo "===================================================="

# check required files
for f in \
  "${VCF}" \
  "${VCF}.tbi" \
  "${POPDIR}/africa.txt" \
  "${POPDIR}/asia.txt" \
  "${POPDIR}/europe.txt" \
  "${POPDIR}/kenya.txt" \
  "${POPDIR}/ethiopia.txt" \
  "${POPDIR}/human.txt" \
  "${POPDIR}/pig.txt"
do
  if [ ! -f "$f" ]; then
    echo "ERROR: Missing required file: $f"
    exit 1
  fi
done

echo
echo "Cleaning population files and checking formatting..."

for f in \
  "${POPDIR}/africa.txt" \
  "${POPDIR}/asia.txt" \
  "${POPDIR}/europe.txt" \
  "${POPDIR}/kenya.txt" \
  "${POPDIR}/ethiopia.txt" \
  "${POPDIR}/human.txt" \
  "${POPDIR}/pig.txt"
do
  # remove Windows carriage returns
  sed -i 's/\r$//' "$f"

  # strip trailing whitespace
  sed -i 's/[[:space:]]*$//' "$f"

  # replace one or more spaces with a tab
  perl -i -pe 's/ +/\t/g' "$f"

  echo
  echo "Preview: $(basename "$f")"
  head -3 "$f"
  echo "Field counts:"
  awk -F'\t' '{print NF}' "$f" | sort | uniq -c
done

echo
echo "VCF sample count:"
bcftools query -l "${VCF}" | wc -l

echo
echo "Checking pixy is available:"
which pixy

for popfile in africa asia europe kenya ethiopia human pig
do
  echo
  echo "----------------------------------------------------"
  echo "Running PIXY for ${popfile}"
  echo "Start: $(date)"
  echo "----------------------------------------------------"

  pixy \
    --stats pi tajima_d \
    --vcf "${VCF}" \
    --populations "${POPDIR}/${popfile}.txt" \
    --window_size 10000 \
    --n_cores "${SLURM_CPUS_PER_TASK}" \
    --output_folder "${OUTDIR}" \
    --output_prefix "${popfile}"

  echo "Finished ${popfile}: $(date)"
done

echo
echo "All PIXY diversity runs complete: $(date)"
echo "Output files:"
ls -lh "${OUTDIR}"
echo "===================================================="
```


## Cluster code: 06 pixy diversity nonafrica

# Purpose: run or summarise PIXY windowed statistics.

```bash
#!/bin/bash
#SBATCH --job-name=pixy_nonafrica
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=pixy_nonafrica.%j.out
#SBATCH --error=pixy_nonafrica.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate pixy_env

VCF=/parallel_scratch/lw01165/05_analysis/05c_allsites_146/Ascaris.autosomes19.ALLSITES.146.vcf.gz
POPDIR=/parallel_scratch/lw01165/05_analysis/05i_pixy/popfiles_diversity
OUTDIR=/parallel_scratch/lw01165/05_analysis/05i_pixy/results_diversity

mkdir -p ${OUTDIR}

echo "===================================================="
echo "Running PIXY for NON-AFRICA"
echo "Start: $(date)"
echo "Host: $(hostname)"
echo "===================================================="

#----------------------------------------------------
# Check required files
#----------------------------------------------------

for f in \
  "${VCF}" \
  "${VCF}.tbi" \
  "${POPDIR}/nonafrica.txt"
do
  if [ ! -f "$f" ]; then
    echo "ERROR: Missing required file: $f"
    exit 1
  fi
done

#----------------------------------------------------
# Clean and standardise nonafrica.txt
#----------------------------------------------------

echo
echo "Cleaning and standardising nonafrica.txt..."

# remove Windows carriage returns
sed -i 's/\r$//' "${POPDIR}/nonafrica.txt"

# strip trailing whitespace
sed -i 's/[[:space:]]*$//' "${POPDIR}/nonafrica.txt"

# replace spaces with tabs
perl -i -pe 's/ +/\t/g' "${POPDIR}/nonafrica.txt"

# enforce exactly two columns and set population label
awk -F'\t' '{
  if (NF >= 1) {
    print $1 "\tNonafrica"
  }
}' "${POPDIR}/nonafrica.txt" > "${POPDIR}/nonafrica.tmp"

mv "${POPDIR}/nonafrica.tmp" "${POPDIR}/nonafrica.txt"

echo
echo "Preview of cleaned nonafrica.txt:"
head -5 "${POPDIR}/nonafrica.txt"

echo "Field counts:"
awk -F'\t' '{print NF}' "${POPDIR}/nonafrica.txt" | sort | uniq -c

#----------------------------------------------------
# Run PIXY
#----------------------------------------------------

echo
echo "Running PIXY..."
echo "Start: $(date)"

pixy \
  --stats pi tajima_d \
  --vcf "${VCF}" \
  --populations "${POPDIR}/nonafrica.txt" \
  --window_size 10000 \
  --n_cores "${SLURM_CPUS_PER_TASK}" \
  --output_folder "${OUTDIR}" \
  --output_prefix "nonafrica"

echo "Finished: $(date)"

echo
echo "Output files:"
ls -lh "${OUTDIR}" | grep nonafrica

echo "===================================================="
```


## Cluster code: 07 pixy diversity all 146

# Purpose: run or summarise PIXY windowed statistics.

```bash
#!/bin/bash
#SBATCH --job-name=pixy_diversity_all146
#SBATCH --partition=high_mem
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=pixy_diversity_all146.%j.out
#SBATCH --error=pixy_diversity_all146.%j.err

set -eo pipefail

source /users/lw01165/miniconda3/etc/profile.d/conda.sh
conda activate pixy_env

VCF=/parallel_scratch/lw01165/05_analysis/05c_allsites_146/Ascaris.autosomes19.ALLSITES.146.vcf.gz
POPFILE=/parallel_scratch/lw01165/05_analysis/05i_pixy/popfiles_diversity/all_assigned_samples.txt
OUTDIR=/parallel_scratch/lw01165/05_analysis/05i_pixy/results_diversity

mkdir -p "${OUTDIR}"

echo "===================================================="
echo "Running PIXY for all 146 samples"
echo "Start: $(date)"
echo "===================================================="

pixy \
  --stats pi tajima_d \
  --vcf "${VCF}" \
  --populations "${POPFILE}" \
  --window_size 10000 \
  --n_cores "${SLURM_CPUS_PER_TASK}" \
  --output_folder "${OUTDIR}" \
  --output_prefix all_146

echo "Finished: $(date)"
echo "Output files:"
ls -lh "${OUTDIR}" | grep all_146
```


## Cluster code: 08 pixy summary awk

Purpose: run or summarise PIXY windowed statistics.

> **NB** This AWK summary was used as a quick check, but full 95% CI come from  R bootstrap script below.

```bash
#!/bin/bash
#SBATCH --job-name=pixy_summary
#SBATCH --output=pixy_summary.%j.out
#SBATCH --error=pixy_summary.%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G

set -eo pipefail

BASEDIR=/parallel_scratch/lw01165/05_analysis/05i_pixy

FST_DIR=${BASEDIR}/results
DIV_DIR=${BASEDIR}/results_diversity
OUTDIR=${BASEDIR}/pixy_all_summary

mkdir -p ${OUTDIR}

echo "Starting PIXY summary"
date

########################################
# 1. Genome-wide mean FST
########################################

echo "Summarising genome-wide FST"

for f in ${FST_DIR}/*fst.txt
do
    name=$(basename ${f} .txt)

    awk '
    NR>1 && $6!="NA" {
        sum+=$6
        n++
    }
    END{
        if(n>0){
            print "'$name'", sum/n
        }
    }' ${f}

done > ${OUTDIR}/fst_genomewide_means.txt


########################################
# 2. Windowed FST summary
########################################

echo "Summarising windowed FST"

for f in ${FST_DIR}/*fst.txt
do
    name=$(basename ${f} .txt)

    awk '
    NR>1 && $6!="NA" {print $6}
    ' ${f} | \
    awk -v pop=${name} '
    {
        sum+=$1
        n++
        if(NR==1 || $1<min){min=$1}
        if(NR==1 || $1>max){max=$1}
    }
    END{
        print pop, n, sum/n, min, max
    }'

done > ${OUTDIR}/fst_window_summary.txt


########################################
# 3. Nucleotide diversity (π)
########################################

echo "Summarising nucleotide diversity"

for f in ${DIV_DIR}/*pi.txt
do
    name=$(basename ${f} .txt)

    awk '
    NR>1 && $5!="NA" {
        sum+=$5
        n++
    }
    END{
        if(n>0){
            print "'$name'", sum/n
        }
    }' ${f}

done > ${OUTDIR}/pi_genomewide_means.txt


########################################
# 4. Tajima's D
########################################

echo "Summarising Tajimas D"

for f in ${DIV_DIR}/*tajima_d.txt
do
    name=$(basename ${f} .txt)

    awk '
    NR>1 && $5!="NA" {
        sum+=$5
        n++
    }
    END{
        if(n>0){
            print "'$name'", sum/n
        }
    }' ${f}

done > ${OUTDIR}/tajimasD_genomewide_means.txt


########################################
# 5. Window distributions for plotting
########################################

echo "Extracting window values for plotting"

cat ${FST_DIR}/*fst.txt | awk 'NR>1 {print $1,$2,$3,$6}' \
> ${OUTDIR}/fst_all_windows.txt

cat ${DIV_DIR}/*pi.txt | awk 'NR>1 {print $1,$2,$3,$5}' \
> ${OUTDIR}/pi_all_windows.txt

cat ${DIV_DIR}/*tajima_d.txt | awk 'NR>1 {print $1,$2,$3,$5}' \
> ${OUTDIR}/tajimasD_all_windows.txt


echo "Finished PIXY summary"
date
```


## RStudio code: genome-wide means and 10,000 bootstrap confidence intervals

This local R code summarises the 10 kb windowed outputs into genome-wide means and non-parametric bootstrap confidence intervals. 

> **REVIEW NOTE:** This is the code section to cite for the 10,000 bootstrap CIs in the methods. It should be kept distinct from the candidate-region scripts, which used 1,000 bootstrap replicates for internal candidate-analysis summaries.

```r

setwd("S:/Research/Lauren_PopGen/pi TajimasD")

base_dir <- "/parallel_scratch/lw01165/pi_TajimasD/text_files"
metadata_file <- "/parallel_scratch/lw01165/pi_TajimasD/sample_QC_list.csv"
outdir <- file.path(base_dir, "bootstrap_CI_outputs")

list.files("S:/Research/Lauren_PopGen/pi TajimasD/text files")

# ============================================================
# Bootstrap 95% CIs for genome-wide FST, dXY, pi and Tajima's D
# Using existing 10 kb window text files
# ============================================================

set.seed(123)

# Use 1000 first to test, then change to 10000 for final
n_boot <- 10000

# Input folders
base_dir <- "S:/Research/Lauren_PopGen/pi TajimasD/text files"
metadata_file <- "S:/Research/Lauren_PopGen/sample_QC_list.csv"

# Output folder
outdir <- file.path(base_dir, "bootstrap_CI_outputs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

read_table_auto <- function(file) {
  message("Reading: ", basename(file))
  
  x <- tryCatch(
    read.delim(file, header = TRUE, sep = "\t", check.names = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(x) || ncol(x) == 1) {
    x <- read.table(file, header = TRUE, check.names = FALSE)
  }
  
  return(x)
}

find_col <- function(dat, candidates) {
  nms <- names(dat)
  nms_clean <- tolower(gsub("[^a-z0-9]", "", nms))
  candidates_clean <- tolower(gsub("[^a-z0-9]", "", candidates))
  
  idx <- match(candidates_clean, nms_clean)
  idx <- idx[!is.na(idx)]
  
  if (length(idx) == 0) {
    stop(
      "Could not find expected column. Tried: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns are: ",
      paste(nms, collapse = ", ")
    )
  }
  
  return(nms[idx[1]])
}

bootstrap_mean_ci <- function(x, n_boot = 10000) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  observed_mean <- mean(x, na.rm = TRUE)
  
  boot_means <- replicate(
    n_boot,
    mean(sample(x, size = length(x), replace = TRUE), na.rm = TRUE)
  )
  
  data.frame(
    mean = observed_mean,
    ci_lower = unname(quantile(boot_means, 0.025, na.rm = TRUE)),
    ci_upper = unname(quantile(boot_means, 0.975, na.rm = TRUE)),
    n_windows = length(x)
  )
}

summarise_metric <- function(file, value_candidates, metric_name) {
  dat <- read_table_auto(file)
  value_col <- find_col(dat, value_candidates)
  
  message("  Using column: ", value_col, " for ", metric_name)
  
  result <- bootstrap_mean_ci(dat[[value_col]], n_boot = n_boot)
  result$metric <- metric_name
  result$value_column_used <- value_col
  
  return(result)
}

format_ci <- function(mean, lower, upper, digits = 4) {
  paste0(
    round(mean, digits),
    " (",
    round(lower, digits),
    "-",
    round(upper, digits),
    ")"
  )
}

# ------------------------------------------------------------
# Load metadata for sample sizes
# ------------------------------------------------------------

meta <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)

# Standardise metadata column names
names(meta) <- tolower(trimws(names(meta)))

# Check required columns exist
required_cols <- c("region", "country", "host")
missing_cols <- setdiff(required_cols, names(meta))

if (length(missing_cols) > 0) {
  stop(
    "Metadata file is missing these columns: ",
    paste(missing_cols, collapse = ", "),
    "\nAvailable columns are: ",
    paste(names(meta), collapse = ", ")
  )
}

# Clean metadata values
meta$region <- trimws(meta$region)
meta$country <- trimws(meta$country)
meta$host <- trimws(meta$host)

count_group <- function(group_type, group_value) {
  
  vals <- meta[[group_type]]
  vals_clean <- tolower(trimws(vals))
  group_clean <- tolower(trimws(group_value))
  
  # Special handling for NonAfrica
  if (group_type == "region" && group_clean %in% c("nonafrica", "non-africa", "non africa")) {
    return(sum(vals_clean %in% c("asia", "europe"), na.rm = TRUE))
  }
  
  sum(vals_clean == group_clean, na.rm = TRUE)
}

# ============================================================
# PART 1: Pairwise FST and dXY
# ============================================================

pairwise <- data.frame(
  comparison = c(
    "Africa vs Asia",
    "Africa vs Europe",
    "Africa vs non-Africa",
    "Asia vs Europe",
    "Kenya vs Ethiopia",
    "Human vs Pig"
  ),
  
  group1_type = c("region", "region", "region", "region", "country", "host"),
  group1_value = c("Africa", "Africa", "Africa", "Asia", "Kenya", "Human"),
  
  group2_type = c("region", "region", "region", "region", "country", "host"),
  group2_value = c("Asia", "Europe", "NonAfrica", "Europe", "Ethiopia", "Pig"),
  
  fst_file = c(
    "africa_vs_asia_fst.txt",
    "africa_vs_europe_fst.txt",
    "africa_vs_nonafrica_fst.txt",
    "asia_vs_europe_fst.txt",
    "kenya_vs_ethiopia_fst.txt",
    "human_vs_pig_fst.txt"
  ),
  
  dxy_file = c(
    "africa_vs_asia_dxy.txt",
    "africa_vs_europe_dxy.txt",
    "africa_vs_nonafrica_dxy.txt",
    "asia_vs_europe_dxy.txt",
    "kenya_vs_ethiopia_dxy.txt",
    "human_vs_pig_dxy.txt"
  ),
  
  stringsAsFactors = FALSE
)

pairwise_results <- list()

for (i in seq_len(nrow(pairwise))) {
  
  comp <- pairwise[i, ]
  
  fst_path <- file.path(base_dir, comp$fst_file)
  dxy_path <- file.path(base_dir, comp$dxy_file)
  
  fst_summary <- summarise_metric(
    fst_path,
    value_candidates = c("avg_wc_fst", "weighted_fst", "wc_fst", "fst", "FST"),
    metric_name = "FST"
  )
  
  dxy_summary <- summarise_metric(
    dxy_path,
    value_candidates = c("avg_dxy", "dxy", "DXY"),
    metric_name = "dXY"
  )
  
  combined <- data.frame(
    comparison = comp$comparison,
    group1 = comp$group1_value,
    n1 = count_group(comp$group1_type, comp$group1_value),
    group2 = comp$group2_value,
    n2 = count_group(comp$group2_type, comp$group2_value),
    
    mean_fst = fst_summary$mean,
    fst_ci_lower = fst_summary$ci_lower,
    fst_ci_upper = fst_summary$ci_upper,
    fst_windows = fst_summary$n_windows,
    fst_formatted = format_ci(
      fst_summary$mean,
      fst_summary$ci_lower,
      fst_summary$ci_upper
    ),
    
    mean_dxy = dxy_summary$mean,
    dxy_ci_lower = dxy_summary$ci_lower,
    dxy_ci_upper = dxy_summary$ci_upper,
    dxy_windows = dxy_summary$n_windows,
    dxy_formatted = format_ci(
      dxy_summary$mean,
      dxy_summary$ci_lower,
      dxy_summary$ci_upper
    )
  )
  
  pairwise_results[[i]] <- combined
}

pairwise_summary <- do.call(rbind, pairwise_results)

outdir <- "S:/Research/Lauren_PopGen/pi TajimasD/text files/bootstrap_CI_outputs"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

dir.exists(outdir)

write.csv(
  pairwise_summary,
  file.path(outdir, "pairwise_fst_dxy_bootstrap_summary.csv"),
  row.names = FALSE
)

# ============================================================
# PART 2: Within-population pi and Tajima's D
# ============================================================

within_pop <- data.frame(
  population = c(
    "Africa",
    "Asia",
    "Europe",
    "Kenya",
    "Ethiopia",
    "Human",
    "Pig",
    "non-Africa",
    "All 146"
  ),
  
  group_type = c(
    "region",
    "region",
    "region",
    "country",
    "country",
    "host",
    "host",
    "region",
    "all"
  ),
  
  group_value = c(
    "Africa",
    "Asia",
    "Europe",
    "Kenya",
    "Ethiopia",
    "Human",
    "Pig",
    "NonAfrica",
    "All"
  ),
  
  pi_file = c(
    "africa_pi.txt",
    "asia_pi.txt",
    "europe_pi.txt",
    "kenya_pi.txt",
    "ethiopia_pi.txt",
    "human_pi.txt",
    "pig_pi.txt",
    "nonafrica_pi.txt",
    "all_146_pi.txt"
  ),
  
  tajima_file = c(
    "africa_tajima_d.txt",
    "asia_tajima_d.txt",
    "europe_tajima_d.txt",
    "kenya_tajima_d.txt",
    "ethiopia_tajima_d.txt",
    "human_tajima_d.txt",
    "pig_tajima_d.txt",
    "nonafrica_tajima_d.txt",
    "all_146_tajima_d.txt"
  ),
  
  stringsAsFactors = FALSE
)

count_within_group <- function(group_type, group_value) {
  if (group_type == "all") {
    return(nrow(meta))
  } else {
    return(count_group(group_type, group_value))
  }
}

within_results <- list()

for (i in seq_len(nrow(within_pop))) {
  
  pop <- within_pop[i, ]
  
  pi_path <- file.path(base_dir, pop$pi_file)
  tajima_path <- file.path(base_dir, pop$tajima_file)
  
  pi_summary <- summarise_metric(
    pi_path,
    value_candidates = c("avg_pi", "pi", "PI", "window_pi"),
    metric_name = "pi"
  )
  
  tajima_summary <- summarise_metric(
    tajima_path,
    value_candidates = c("TajimaD", "tajima_d", "tajimas_d", "D", "tajima"),
    metric_name = "TajimaD"
  )
  
  combined <- data.frame(
    population = pop$population,
    n = count_within_group(pop$group_type, pop$group_value),
    
    mean_pi = pi_summary$mean,
    pi_ci_lower = pi_summary$ci_lower,
    pi_ci_upper = pi_summary$ci_upper,
    pi_windows = pi_summary$n_windows,
    pi_formatted = format_ci(
      pi_summary$mean,
      pi_summary$ci_lower,
      pi_summary$ci_upper
    ),
    
    mean_tajimas_d = tajima_summary$mean,
    tajimas_d_ci_lower = tajima_summary$ci_lower,
    tajimas_d_ci_upper = tajima_summary$ci_upper,
    tajimas_d_windows = tajima_summary$n_windows,
    tajimas_d_formatted = format_ci(
      tajima_summary$mean,
      tajima_summary$ci_lower,
      tajima_summary$ci_upper
    )
  )
  
  within_results[[i]] <- combined
}

within_summary <- do.call(rbind, within_results)

write.csv(
  within_summary,
  file.path(outdir, "within_population_pi_tajimasd_bootstrap_summary.csv"),
  row.names = FALSE
)

# ============================================================
# make  simplified tables
# ============================================================

pairwise_thesis_table <- pairwise_summary[, c(
  "comparison",
  "group1",
  "n1",
  "group2",
  "n2",
  "fst_formatted",
  "dxy_formatted",
  "fst_windows",
  "dxy_windows"
)]

names(pairwise_thesis_table) <- c(
  "Comparison",
  "Group_1",
  "n_1",
  "Group_2",
  "n_2",
  "FST_mean_95CI",
  "dXY_mean_95CI",
  "FST_windows",
  "dXY_windows"
)

within_thesis_table <- within_summary[, c(
  "population",
  "n",
  "pi_formatted",
  "tajimas_d_formatted",
  "pi_windows",
  "tajimas_d_windows"
)]

names(within_thesis_table) <- c(
  "Population",
  "n",
  "pi_mean_95CI",
  "TajimasD_mean_95CI",
  "pi_windows",
  "TajimasD_windows"
)

write.csv(
  pairwise_thesis_table,
  file.path(outdir, "THESIS_TABLE_pairwise_fst_dxy.csv"),
  row.names = FALSE
)

write.csv(
  within_thesis_table,
  file.path(outdir, "THESIS_TABLE_within_population_pi_tajimasd.csv"),
  row.names = FALSE
)

# ============================================================
# Print summaries
# ============================================================

cat("\nPairwise FST/dXY summary:\n")
print(pairwise_thesis_table)

cat("\nWithin-population pi/Tajima's D summary:\n")
print(within_thesis_table)

cat("\nDone.\n")
cat("Files written to:\n")
cat(outdir, "\n")
```

---

## Outputs and QC checks

10-kb FST, dXY, pi, Tajima’s D and bootstrap summaries. Consult recorded logs and intermediate tables to distinguish expected settings from actual execution success.

## Provenance and limitations

- Important: historical thesis PIXY branch; corrected all-sites VCF was NOT used to produce these results.
- Sample identifiers and host/country metadata were manually curated in `sample_QC_list.csv`; it is not generated by an algorithm in these files.
- For the discovery and correction of the PIXY all-sites input issue, see [`99_PIXY_allsites_correction.md`](99_PIXY_allsites_correction.md).
