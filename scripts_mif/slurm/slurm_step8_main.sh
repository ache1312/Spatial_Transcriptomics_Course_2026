#!/bin/bash
#SBATCH -J step8_survival
#SBATCH -p main
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/step8-main-%j.out
#SBATCH --error=logs/step8-main-%j.err

set -euo pipefail

PIPELINE_DIR="${PIPELINE_DIR:-/path/to/pipeline_leiden_validated}"
CONDA_ENV="${CONDA_ENV:-/path/to/.conda/envs/spatial-biology}"

cd "$PIPELINE_DIR"
mkdir -p logs

missing=0
for f in \
  analysis_imc/ERpos_directed_neighbors_r40_coarse.csv \
  analysis_imc/ERpos_nhood_enrichment_r25.csv \
  analysis_imc/ERpos_dispersion_metrics_by_patient.csv
do
  if [ ! -f "$f" ]; then
    echo "ERROR: Missing prerequisite file: $f" >&2
    missing=1
  fi
done
if [ "$missing" -ne 0 ]; then
  echo "ERROR: Step 8 prerequisites missing (run Steps 4, 5, 7 first)." >&2
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

conda run -p "$CONDA_ENV" python -u scripts/step8_survival.py

