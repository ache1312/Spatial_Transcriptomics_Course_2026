#!/bin/bash
#---------------Script SBATCH - NLHPC ----------------
#SBATCH -J xenium_cellchat
#SBATCH -p labs
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=10668
#SBATCH --time=24:00:00
#SBATCH --mail-user=ahernandez@cienciavida.org
#SBATCH --mail-type=ALL

#--------------- Track start time ---------------------
start_time=$(date +%s)

#-----------------Modules------------------------------

ml purge
ml miniconda3
eval "$(conda shell.bash hook)"
conda activate SPATIAL


#--------------- Create folder for plots ---------------
PLOT_DIR="/home/courses/student75/xenium_plots"
mkdir -p "$PLOT_DIR"
REPORT_DIR="/home/courses/student75/xenium_reports"
mkdir -p "$REPORT_DIR"

# Report file (unique per job)
REPORT_FILE="${REPORT_DIR}/xenium_cellchat_report_${SLURM_JOB_ID}.log"

echo "Plot directory: $PLOT_DIR"
echo "Report will be saved to: $REPORT_FILE"

#--------------- Run R + capture console + time --------
{
  echo "========== Xenium CellChat report =========="
  echo "Job ID: ${SLURM_JOB_ID}"
  echo "Node: $(hostname)"
  echo "Start time: $(date -d "@${start_time}" '+%Y-%m-%d %H:%M:%S')"
  echo "Plot directory: ${PLOT_DIR}"
  echo
  echo "----- R console output -----"

  Rscript /home/courses/student75/xenium_cellchat.R "$PLOT_DIR"
  rstatus=$?

  echo "----- End of R console output -----"
  echo
  end_time=$(date +%s)
  runtime=$((end_time - start_time))

  echo "End time: $(date -d "@${end_time}" '+%Y-%m-%d %H:%M:%S')"
  echo "Rscript exit status: ${rstatus}"
  printf 'Total runtime: %02dh:%02dm:%02ds\n' \
    $((runtime/3600)) $((runtime%3600/60)) $((runtime%60))
} &> "$REPORT_FILE"

echo "Report written to: $REPORT_FILE"
