#!/bin/bash
#SBATCH -J IMC_pipeline
#SBATCH -p main
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@example.com
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

# IMPORTANT: SLURM runs a spooled copy of this script, so `BASH_SOURCE[0]` is NOT
# the repo path. Use `PIPELINE_DIR` or fall back to `SLURM_SUBMIT_DIR`.
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-/home/courses/student75/mif_workshop}}"

# Prefer explicit CONDA_ENV, otherwise use the student env path (matches run_step_slurm.sh)
if [ -z "${CONDA_ENV:-}" ]; then
    CONDA_ENV="/home/courses/student75/.conda/envs/spatial-biology"
fi

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

bash slurm/run_pipeline_complete_v2.sh "$CONDA_ENV"
