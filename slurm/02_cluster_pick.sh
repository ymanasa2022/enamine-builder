#!/bin/bash
#SBATCH --job-name=enamine_cluster
#SBATCH --account=ACCOUNT          # <-- set your cluster account
#SBATCH --partition=standard       # <-- set your partition
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8          # UMAP and KMeans use threading
#SBATCH --mem=64G
#SBATCH --output=logs/02_cluster_pick_%j.out
#SBATCH --error=logs/02_cluster_pick_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ymanasa@umich.edu

# ── Environment ───────────────────────────────────────────────────────────────
source ~/.bashrc
conda activate enamine-builder

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SIM_DIR="${SIM_DIR:-${SCRIPT_DIR}/results/reference_filter}"   # output of step 01b
OUT_DIR="${OUT_DIR:-${SCRIPT_DIR}/results/cluster_pick}"

mkdir -p "${OUT_DIR}" "${SCRIPT_DIR}/logs"

echo "SCRIPT_DIR : ${SCRIPT_DIR}"
echo "SIM_DIR    : ${SIM_DIR}"
echo "OUT_DIR    : ${OUT_DIR}"
echo "CPUS       : ${SLURM_CPUS_PER_TASK}"
echo "Host       : $(hostname)"
echo "Started    : $(date)"

# ── Run ───────────────────────────────────────────────────────────────────────
# --n_clusters 500: aim for ~500 diverse clusters; adjust based on how many
#                   compounds survive 1b_reference_filter (aim for ~1-2K picks)
python3 "${SCRIPT_DIR}/2_cluster_pick.py" \
    --fingerprints "${SIM_DIR}/similar_fps.npy" \
    --ligands      "${SIM_DIR}/similar_ligands.csv" \
    --out_dir      "${OUT_DIR}" \
    --n_clusters   500 \
    --seed         42

echo "Finished : $(date)"
