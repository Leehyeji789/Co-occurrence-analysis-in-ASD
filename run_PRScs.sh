#!/usr/bin/env bash
#
# PRS-CS genome-wide polygenic score calculation script
# Usage:
#   bash run_PRScs.sh <sumstats> <bim_prefix> <trait> <N_GWAS> <phi> <out_dir>
#
# Example:
#   bash run_PRScs.sh \
#     /path/to/ldblk_ukbb_eur.ASD.tsv \
#     /path/to/genotype_prefix \
#     ASD \
#     58948 \
#     1e-2 \
#     /path/to/prs_output
#

set -euo pipefail

# -----------------------------
# Parse input arguments
# -----------------------------
SUMSTATS="$1"       # GWAS summary statistics (formatted for PRS-CS)
BIM_PREFIX="$2"     # PLINK prefix (without .bed/.bim/.fam)
TRAIT="$3"          # Trait name (for labeling outputs)
N_GWAS="$4"         # Effective GWAS sample size
PHI="${5:-1e-2}"    # Global shrinkage parameter (default: 1e-2)
OUT_DIR="$6"        # Output directory

mkdir -p "${OUT_DIR}/${TRAIT}"

# -----------------------------
# PRS-CS configuration
# -----------------------------
PRSCS_DIR="/data2/resources/PRS/PRScs"
REF_DIR="${PRSCS_DIR}/ldblk_ukbb_eur"
LOG_FILE="${OUT_DIR}/${TRAIT}_PRScs_log.txt"

echo "Running PRS-CS for trait: ${TRAIT}" > "${LOG_FILE}"

# -----------------------------
# Helper function
# -----------------------------
run_prscs_block() {
  local start_chr="$1"
  local end_chr="$2"
  local mkl_threads="$3"
  local label="$4"

  echo "[$(date)] Start PRS-CS for ${TRAIT}, chr${start_chr}-${end_chr} (${label})" | tee -a "${LOG_FILE}"
  local start_time
  start_time=$(date +%s)

  export MKL_NUM_THREADS="${mkl_threads}"
  export NUMEXPR_NUM_THREADS="${mkl_threads}"
  export OMP_NUM_THREADS="${mkl_threads}"

  for chr in $(seq "${start_chr}" "${end_chr}"); do
    python3 "${PRSCS_DIR}/PRScs.py" \
      --ref_dir="${REF_DIR}" \
      --bim_prefix="${BIM_PREFIX}" \
      --sst_file="${SUMSTATS}" \
      --n_gwas="${N_GWAS}" \
      --phi="${PHI}" \
      --chrom="${chr}" \
      --out_dir="${OUT_DIR}/${TRAIT}" \
      --seed=52 &
  done

  wait

  local end_time
  end_time=$(date +%s)
  echo "[$(date)] Finished PRS-CS for chr${start_chr}-${end_chr} (${label}); elapsed: $((end_time - start_time)) seconds" | tee -a "${LOG_FILE}"
}

# -----------------------------
# Run PRS-CS in three chromosome blocks
# (tuned for typical CPU/memory settings)
# -----------------------------

# Block 1: chr 1–5
run_prscs_block 1 5 8 "block1"

# Block 2: chr 6–13
run_prscs_block 6 13 5 "block2"

# Block 3: chr 14–22
run_prscs_block 14 22 4 "block3"

# -----------------------------
# Merge per-chromosome PRS-CS outputs
# -----------------------------
PST_PREFIX="${OUT_DIR}/${TRAIT}_pst_eff_a1_b0.5_phi${PHI}"

echo "[$(date)] Merging PRS-CS weight files into ${PST_PREFIX}.txt" | tee -a "${LOG_FILE}"

cat "${OUT_DIR}/${TRAIT}"/*_pst_eff_a1_b0.5_phi"${PHI}"_chr*.txt > "${PST_PREFIX}.txt"

# -----------------------------
# Compute polygenic scores with PLINK
# -----------------------------
echo "[$(date)] Computing PRS with PLINK" | tee -a "${LOG_FILE}"

plink \
  --bfile "${BIM_PREFIX}" \
  --score "${PST_PREFIX}.txt" 2 4 6 \
  --out "${OUT_DIR}/PRS.${TRAIT}_pst_eff_a1_b0.5_phi${PHI}"

echo "[$(date)] PRS-CS pipeline for ${TRAIT} completed." | tee -a "${LOG_FILE}"