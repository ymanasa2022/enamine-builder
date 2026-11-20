#!/bin/bash
#SBATCH --job-name=enamine_fp
#SBATCH --account=ACCOUNT          # <-- set your cluster account
#SBATCH --partition=standard       # <-- set your partition
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=logs/01_fingerprint_%j.out
#SBATCH --error=logs/01_fingerprint_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ymanasa@umich.edu

# ── Environment ───────────────────────────────────────────────────────────────
# Activate conda env — adjust path/name to match your cluster setup
source ~/.bashrc
conda activate enamine-builder    # or: module load python/3.11; pip install -r requirements.txt

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${DATA_DIR:-${SCRIPT_DIR}/../../Data/ENAMINE/20261022}"   # override with: DATA_DIR=/path/to/data sbatch ...
OUT_DIR="${OUT_DIR:-${SCRIPT_DIR}/results/fingerprints}"

mkdir -p "${OUT_DIR}" "${SCRIPT_DIR}/logs"

echo "SCRIPT_DIR : ${SCRIPT_DIR}"
echo "DATA_DIR   : ${DATA_DIR}"
echo "OUT_DIR    : ${OUT_DIR}"
echo "Host       : $(hostname)"
echo "Started    : $(date)"

# ── Run ───────────────────────────────────────────────────────────────────────
python3 "${SCRIPT_DIR}/1_fingerprint.py" \
    --input   "${DATA_DIR}/enamine_ligands.csv" \
    --out_dir "${OUT_DIR}"

echo "Finished : $(date)"
