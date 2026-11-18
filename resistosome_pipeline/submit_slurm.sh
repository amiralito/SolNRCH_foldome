#!/bin/bash
#SBATCH --job-name=resistosome_analysis
#SBATCH --output=resistosome_%j.out
#SBATCH --error=resistosome_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=batch

# ── Configuration ──
INPUT_DIR="/path/to/your/structures"
OUTPUT_DIR="/path/to/output"
HMM="/path/to/MADA.hmm"
MHD_CSV="/path/to/NRCH_MHD_nlrexp_coords.csv"
PLOOP_CSV="/path/to/NRCH_ploop_nlrexp_coords.csv"
NBARC_CSV="/path/to/NRCH_clade_filtered_len_NBARC_filtered_95_coords.csv"
PIPELINE_DIR="/path/to/resistosome_pipeline"

# ── Environment ──
# module load python/3.12   # uncomment if needed
# source /path/to/venv/bin/activate  # uncomment if using virtualenv
# mkdssp and hmmsearch must be on PATH

# ── Run ──
cd "${PIPELINE_DIR}"

python run_pipeline.py \
    --input_dir "${INPUT_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --hmm "${HMM}" \
    --mhd_csv "${MHD_CSV}" \
    --ploop_csv "${PLOOP_CSV}" \
    --nbarc_csv "${NBARC_CSV}" \
    --contact_cutoff 4.0 \
    --pae_cutoff 10.0 \
    --k 1.96 \
    --workers "${SLURM_CPUS_PER_TASK:-8}" \
    --log_level INFO

echo "Pipeline finished at $(date)"
