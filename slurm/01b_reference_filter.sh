#!/bin/bash
#SBATCH --job-name=enamine_reffilter
#SBATCH --account=ACCOUNT          # <-- set your cluster account
#SBATCH --partition=standard       # <-- set your partition
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=logs/01b_reference_filter_%j.out
#SBATCH --error=logs/01b_reference_filter_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ymanasa@umich.edu

# ── Environment ───────────────────────────────────────────────────────────────
source ~/.bashrc
conda activate enamine-builder

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FP_DIR="${FP_DIR:-${SCRIPT_DIR}/results/fingerprints}"   # output of step 01
OUT_DIR="${OUT_DIR:-${SCRIPT_DIR}/results/reference_filter}"

mkdir -p "${OUT_DIR}" "${SCRIPT_DIR}/logs"

echo "SCRIPT_DIR : ${SCRIPT_DIR}"
echo "FP_DIR     : ${FP_DIR}"
echo "OUT_DIR    : ${OUT_DIR}"
echo "Host       : $(hostname)"
echo "Started    : $(date)"

# ── Run ───────────────────────────────────────────────────────────────────────
# --require_basic_n: keep compounds with a basic sp3 N (engages TgCPL Glu75/Asp218)
# --min_sim 0.30:    Tanimoto >= 0.30 to any of Cpd1/12/13/18
python3 "${SCRIPT_DIR}/1b_reference_filter.py" \
    --input   "${FP_DIR}/filtered_ligands.csv" \
    --fps     "${FP_DIR}/fingerprints.npy" \
    --out_dir "${OUT_DIR}" \
    --min_sim 0.30 \
    --require_basic_n

echo "Finished : $(date)"
