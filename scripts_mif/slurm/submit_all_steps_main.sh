#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

PIPELINE_DIR="${PIPELINE_DIR:-$REPO_DIR}"
export PIPELINE_DIR

if [ -z "${CONDA_ENV:-}" ]; then
  echo "ERROR: Set CONDA_ENV (path to conda env) before submitting." >&2
  echo "Example: export CONDA_ENV=/path/to/.conda/envs/spatial-biology" >&2
  exit 1
fi
export CONDA_ENV

echo "Submitting workshop steps (MAIN) with dependencies..."
echo "PIPELINE_DIR=$PIPELINE_DIR"
echo "CONDA_ENV=$CONDA_ENV"

job_step0="$(sbatch --parsable "$SCRIPT_DIR/slurm_step0_main.sh")"
job_step1="$(sbatch --parsable "$SCRIPT_DIR/slurm_step1_main.sh")"
job_step2="$(sbatch --parsable "$SCRIPT_DIR/slurm_step2_main.sh")"

job_step3="$(sbatch --parsable "$SCRIPT_DIR/slurm_step3_main.sh")"
job_step4="$(sbatch --parsable --dependency=afterok:"$job_step3" "$SCRIPT_DIR/slurm_step4_main.sh")"

job_step5="$(sbatch --parsable "$SCRIPT_DIR/slurm_step5_main.sh")"
job_step6="$(sbatch --parsable "$SCRIPT_DIR/slurm_step6_main.sh")"
job_step7="$(sbatch --parsable "$SCRIPT_DIR/slurm_step7_main.sh")"

job_step8="$(sbatch --parsable --dependency=afterok:"$job_step4":"$job_step5":"$job_step7" "$SCRIPT_DIR/slurm_step8_main.sh")"

job_step9="$(sbatch --parsable "$SCRIPT_DIR/slurm_step9_main.sh")"
job_step10="$(sbatch --parsable "$SCRIPT_DIR/slurm_step10_main.sh")"
job_step11="$(sbatch --parsable "$SCRIPT_DIR/slurm_step11_main.sh")"
job_step12="$(sbatch --parsable "$SCRIPT_DIR/slurm_step12_main.sh")"
job_step13="$(sbatch --parsable "$SCRIPT_DIR/slurm_step13_main.sh")"

cat <<EOF
Submitted jobs:
  step0:  $job_step0
  step1:  $job_step1
  step2:  $job_step2
  step3:  $job_step3
  step4:  $job_step4   (afterok: step3)
  step5:  $job_step5
  step6:  $job_step6
  step7:  $job_step7
  step8:  $job_step8   (afterok: step4,step5,step7)
  step9:  $job_step9
  step10: $job_step10
  step11: $job_step11
  step12: $job_step12
  step13: $job_step13
EOF

