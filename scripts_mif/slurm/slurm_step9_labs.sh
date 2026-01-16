#!/bin/bash
#SBATCH -J step9_ripley
#SBATCH -p labs
#SBATCH --cpus-per-task=15
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/step9-labs-%j.out
#SBATCH --error=logs/step9-labs-%j.err

set -euo pipefail

PIPELINE_DIR="${PIPELINE_DIR:-/path/to/pipeline_leiden_validated}"
CONDA_ENV="${CONDA_ENV:-/path/to/.conda/envs/spatial-biology}"

cd "$PIPELINE_DIR"
mkdir -p logs

NUM_THREADS="${NUM_THREADS:-${SLURM_CPUS_PER_TASK:-15}}"
export NUM_THREADS
export OMP_NUM_THREADS="$NUM_THREADS"
export OPENBLAS_NUM_THREADS="$NUM_THREADS"
export MKL_NUM_THREADS="$NUM_THREADS"
export NUMEXPR_NUM_THREADS="$NUM_THREADS"
export PYTHONUNBUFFERED=1
unset NUMBA_DISABLE_JIT
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-/tmp}"

conda run -p "$CONDA_ENV" python -u scripts/step9_spatial_point_process.py
