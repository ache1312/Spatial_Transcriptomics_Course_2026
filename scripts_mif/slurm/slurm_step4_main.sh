#!/bin/bash
#SBATCH -J step4_interactions
#SBATCH -p main
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=logs/step4-main-%j.out
#SBATCH --error=logs/step4-main-%j.err

set -euo pipefail

PIPELINE_DIR="${PIPELINE_DIR:-/path/to/pipeline_leiden_validated}"
CONDA_ENV="${CONDA_ENV:-/path/to/.conda/envs/spatial-biology}"

cd "$PIPELINE_DIR"
mkdir -p logs

if [ ! -f "data/OUTPUT_imc_ERpos_with_graphs.h5ad" ]; then
  echo "ERROR: Missing data/OUTPUT_imc_ERpos_with_graphs.h5ad (run Step 3 first)" >&2
  exit 1
fi

NUM_THREADS="${NUM_THREADS:-${SLURM_CPUS_PER_TASK:-20}}"
export NUM_THREADS
export OMP_NUM_THREADS="$NUM_THREADS"
export OPENBLAS_NUM_THREADS="$NUM_THREADS"
export MKL_NUM_THREADS="$NUM_THREADS"
export NUMEXPR_NUM_THREADS="$NUM_THREADS"
export PYTHONUNBUFFERED=1
unset NUMBA_DISABLE_JIT
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-/tmp}"

conda run -p "$CONDA_ENV" python -u scripts/step4_interactions.py

