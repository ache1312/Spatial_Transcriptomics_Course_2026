#!/bin/bash
#SBATCH -J step0_dataset
#SBATCH -p main
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=logs/step0-main-%j.out
#SBATCH --error=logs/step0-main-%j.err

set -euo pipefail

PIPELINE_DIR="${PIPELINE_DIR:-/path/to/pipeline_leiden_validated}"
CONDA_ENV="${CONDA_ENV:-/path/to/.conda/envs/spatial-biology}"

cd "$PIPELINE_DIR"
mkdir -p logs

NUM_THREADS="${NUM_THREADS:-${SLURM_CPUS_PER_TASK:-20}}"
export NUM_THREADS
export OMP_NUM_THREADS="$NUM_THREADS"
export OPENBLAS_NUM_THREADS="$NUM_THREADS"
export MKL_NUM_THREADS="$NUM_THREADS"
export NUMEXPR_NUM_THREADS="$NUM_THREADS"
export PYTHONUNBUFFERED=1
unset NUMBA_DISABLE_JIT
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-/tmp}"

conda run -p "$CONDA_ENV" python -u scripts/step0_describe_dataset.py

