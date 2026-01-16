#!/bin/bash
#===============================================================================
# SPATIAL BIOLOGY PIPELINE - COMPLETE END-TO-END (VALIDATED VERSION)
#===============================================================================
#
# Description:
#   Complete ER+ spatial analysis pipeline (14 active steps) with full validation
#   Runs all steps in correct dependency order with time/RAM tracking
#
# Prerequisites:
#   - data/imc_ERpos_with_embeddings_and_full_names.h5ad (unified data file with all metadata)
#   - Conda environment with: scanpy, squidpy, lifelines, scikit-learn, scipy
#
# Pipeline Steps:
#   Steps 0-7, 9-13: Independent (can run in parallel)
#   Step 8: Survival (DEPENDS ON 4,5,7)
#   Steps 14-15: Archived (not part of active pipeline)
#
# Usage:
#   ./run_pipeline_complete_v2.sh [CONDA_ENV_PATH]
#
#   Options:
#     CONDA_ENV_PATH: Path to conda environment (optional, auto-detects if active)
#
# Time:
#   - Complete pipeline: varies by Step 5 runtime (14 active steps)
#
# Peak RAM: ~2-3 GB
#
# Output:
#   - analysis_imc/*.csv (35+ analysis files)
#   - figures/*/*.png (active steps only)
#   - logs/pipeline_YYYYMMDD_HHMMSS.log (detailed execution log)
#   - PIPELINE_EXECUTION_SUMMARY.txt (metrics summary)
#
# Reference: Keren et al. 2018, Eng et al. 2019
#===============================================================================

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------

# Get script directory (slurm/) and repo root (pipeline_leiden_validated/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_DIR"

# Conda environment
if [ -n "${1:-}" ]; then
    CONDA_ENV="$1"
elif [ -n "${CONDA_PREFIX:-}" ]; then
    CONDA_ENV="$CONDA_PREFIX"
else
    echo "ERROR: No conda environment specified and none currently active"
    echo "Usage: $0 [CONDA_ENV_PATH]"
    exit 1
fi

# Create logs directory
mkdir -p logs
TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
LOG_FILE="logs/pipeline_${TIMESTAMP}.log"

# Performance tracking arrays
declare -a STEP_TIMES
declare -a STEP_NAMES
declare -a STEP_STATUS

#-------------------------------------------------------------------------------
# Helper Functions
#-------------------------------------------------------------------------------

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

log_section() {
    log ""
    log "========================================================================"
    log "$1"
    log "========================================================================"
}

run_step() {
    local STEP_NUM=$1
    local STEP_NAME=$2
    local SCRIPT=$3

    log ""
    log "========================================================================"
    log "STEP $STEP_NUM: $STEP_NAME"
    log "Script: scripts/$SCRIPT"
    log "========================================================================"

    local START_TIME=$(date +%s)

    if [ ! -f "scripts/$SCRIPT" ]; then
        log "✗ ERROR: Script not found: scripts/$SCRIPT"
        STEP_STATUS+=("FAILED")
        STEP_TIMES+=(0)
        STEP_NAMES+=("Step $STEP_NUM: $STEP_NAME")
        exit 1
    fi

    # Run with time and resource tracking
    set +e  # Temporarily disable exit on error to capture failures

    /usr/bin/time -v conda run -p "$CONDA_ENV" \
        env NUMBA_CACHE_DIR=/tmp PYTHONUNBUFFERED=1 \
        python -u "scripts/$SCRIPT" 2>&1 | tee -a "$LOG_FILE"

    local EXIT_CODE=${PIPESTATUS[0]}
    set -e

    local END_TIME=$(date +%s)
    local DURATION=$((END_TIME - START_TIME))

    # Store metrics
    STEP_TIMES+=($DURATION)
    STEP_NAMES+=("Step $STEP_NUM: $STEP_NAME")

    if [ $EXIT_CODE -eq 0 ]; then
        log "✓ STEP $STEP_NUM COMPLETED in ${DURATION}s ($(($DURATION / 60))m $(($DURATION % 60))s)"
        STEP_STATUS+=("SUCCESS")
    else
        log "✗ STEP $STEP_NUM FAILED (exit code: $EXIT_CODE) after ${DURATION}s"
        STEP_STATUS+=("FAILED")
        exit 1
    fi
}

validate_file() {
    local FILE=$1
    local DESC=$2

    if [ -f "$FILE" ]; then
        local SIZE_MB=$(du -m "$FILE" | cut -f1)
        log "  ✓ $DESC: $FILE (${SIZE_MB} MB)"
        return 0
    else
        log "  ✗ MISSING: $DESC: $FILE"
        return 1
    fi
}

#-------------------------------------------------------------------------------
# Pre-flight Checks
#-------------------------------------------------------------------------------

log_section "SPATIAL BIOLOGY PIPELINE - END-TO-END VALIDATION"
log "Mode: COMPLETE (14 active steps; steps 14-15 archived)"
log "Working directory: $SCRIPT_DIR"
log "Conda environment: $CONDA_ENV"
log "Log file: $LOG_FILE"
log "Timestamp: $TIMESTAMP"

# Validate conda environment
if [ ! -d "$CONDA_ENV" ]; then
    log "✗ ERROR: Conda environment not found: $CONDA_ENV"
    exit 1
fi

log ""
log "Validating conda environment packages..."
conda run -p "$CONDA_ENV" python -c "
import sys
required = ['scanpy', 'squidpy', 'lifelines', 'sklearn', 'scipy', 'pandas', 'numpy']
missing = []
for pkg in required:
    try:
        __import__(pkg)
    except ImportError:
        missing.append(pkg)
if missing:
    print(f'✗ Missing packages: {missing}')
    sys.exit(1)
else:
    print('✓ All required packages installed')
" 2>&1 | tee -a "$LOG_FILE"

if [ ${PIPESTATUS[0]} -ne 0 ]; then
    log "✗ ERROR: Required packages missing in conda environment"
    exit 1
fi

# Validate input data
log ""
log "Validating input data files..."

ALL_FILES_OK=true

if ! validate_file "data/imc_ERpos_with_embeddings_and_full_names.h5ad" "Unified data file (all metadata)"; then
    ALL_FILES_OK=false
fi

if [ "$ALL_FILES_OK" != true ]; then
    log "✗ ERROR: Required input file missing"
    exit 1
fi

# Inspect main data file
log ""
log "Inspecting main data file..."
conda run -p "$CONDA_ENV" python - <<'PY' 2>&1 | tee -a "$LOG_FILE"
import scanpy as sc
import sys
adata = sc.read_h5ad('data/imc_ERpos_with_embeddings_and_full_names.h5ad')
print(f'  Cells: {adata.n_obs:,}')
print(f'  Markers: {adata.n_vars}')
print(f'  Patients: {adata.obs[\"Patient\"].nunique()}')
print(f'  Cell types: {adata.obs[\"celltype_full\"].nunique()}')
print(f'  Has spatial coords: {"spatial" in adata.obsm}')
print(f'  Has embeddings: {"X_umap" in adata.obsm}')

# Verify clinical metadata
required_clinical = ['Recurrence', 'Recurrence_time', 'Survival', 'Survival_time', 'age', 'tumor_size']
missing = [col for col in required_clinical if col not in adata.obs.columns]
if missing:
    print(f'  ✗ MISSING clinical columns: {missing}')
    sys.exit(1)
else:
    print(f'  ✓ All clinical metadata present')
PY

# Create output directories
log ""
log "Creating output directories..."
mkdir -p analysis_imc figures validation
log "  ✓ Directories created"

# Track total time
TOTAL_START=$(date +%s)

#===============================================================================
# INDEPENDENT STEPS (Steps 0-7, 9-13)
# Can run in any order - no dependencies
#===============================================================================

log_section "PHASE 1: INDEPENDENT ANALYSES"

# Step 0: Dataset Description
run_step 0 "Dataset Description" "step0_describe_dataset.py"

# Step 1: QC Visualization
run_step 1 "QC Visualization" "step1_visualize_imc_qc.py"

# Step 2: Embeddings QC
run_step 2 "Embeddings & QC" "step2_embeddings_and_qc.py"

# Step 3: Spatial Graphs
run_step 3 "Spatial Neighbor Graphs" "step3_neighbors.py"

# Step 4: Cell-Cell Interactions (REQUIRED FOR STEP 8)
run_step 4 "Cell-Cell Interactions" "step4_interactions.py"

# Step 5: Neighborhood Enrichment (REQUIRED FOR STEP 8)
run_step 5 "Neighborhood Enrichment" "step5_enrichment.py"

# Step 6: LDA Topics
run_step 6 "LDA Topic Modeling" "step6_lda.py"

# Step 7: Dispersion Metrics (REQUIRED FOR STEP 8)
run_step 7 "Dispersion & Isolation Metrics" "step7_dispersion.py"

# Step 9: Spatial Point Process
run_step 9 "Spatial Point Process (Ripley's L)" "step9_spatial_point_process.py"

# Step 10: Tumor-Immune Mixing
run_step 10 "Tumor-Immune Mixing Score" "step10_mixing_score.py"

# Step 11: Proliferation-Stratified Survival
run_step 11 "Proliferation-Stratified Survival" "step11_proliferation_stratified_survival.py"

# Step 12: T Cell Clustering
run_step 12 "T Cell Clustering Analysis" "step12_tcell_clustering_analysis.py"

# Step 13: Lineage Neighbors
run_step 13 "Lineage Neighbor Analysis" "step13_lineage_neighbors.py"

#===============================================================================
# DEPENDENT STEPS
# MUST run in specific order due to file dependencies
#===============================================================================

log_section "PHASE 2: DEPENDENT ANALYSES (Respecting Dependencies)"

# Step 8: Survival Analysis (DEPENDS ON: Steps 4, 5, 7)
log "Prerequisites for Step 8:"
validate_file "analysis_imc/ERpos_directed_neighbors_r40_coarse.csv" "Step 4 output"
validate_file "analysis_imc/ERpos_nhood_enrichment_r25.csv" "Step 5 output"
validate_file "analysis_imc/ERpos_dispersion_metrics_by_patient.csv" "Step 7 output"

run_step 8 "Survival Analysis (Univariate Cox PH)" "step8_survival.py"

#===============================================================================
# POST-EXECUTION VALIDATION
#===============================================================================

log_section "POST-EXECUTION VALIDATION"

# Count outputs
log "Counting generated outputs..."

CSV_COUNT=$(ls -1 analysis_imc/ERpos_*.csv 2>/dev/null | wc -l)
FIGURE_DIRS=$(find figures -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l)
FIGURE_COUNT=$(find figures -name "*.png" 2>/dev/null | wc -l)

log "  CSV files: $CSV_COUNT"
log "  Figure directories: $FIGURE_DIRS"
log "  Total figures: $FIGURE_COUNT"

# Validate key outputs exist
log ""
log "Validating key output files..."

CRITICAL_OUTPUTS=(
    "analysis_imc/ERpos_directed_neighbors_r40_coarse.csv:Step 4"
    "analysis_imc/ERpos_nhood_enrichment_r25.csv:Step 5"
    "analysis_imc/ERpos_dispersion_metrics_by_patient.csv:Step 7"
    "analysis_imc/ERpos_survival_cox_univariate.csv:Step 8"
)

VALIDATION_OK=true
for entry in "${CRITICAL_OUTPUTS[@]}"; do
    FILE="${entry%%:*}"
    DESC="${entry##*:}"
    if [ ! -f "$FILE" ]; then
        log "  ✗ MISSING: $FILE ($DESC)"
        VALIDATION_OK=false
    else
        ROWS=$(wc -l < "$FILE")
        log "  ✓ $FILE ($DESC) - $ROWS rows"
    fi
done

if [ "$VALIDATION_OK" != true ]; then
    log "⚠ WARNING: Some critical outputs are missing"
fi

#===============================================================================
# GENERATE EXECUTION SUMMARY
#===============================================================================

TOTAL_END=$(date +%s)
TOTAL_DURATION=$((TOTAL_END - TOTAL_START))
TOTAL_MINUTES=$((TOTAL_DURATION / 60))
TOTAL_SECONDS=$((TOTAL_DURATION % 60))

log_section "GENERATING EXECUTION SUMMARY"

SUMMARY_FILE="PIPELINE_EXECUTION_SUMMARY.txt"

{
    echo "================================================================================"
    echo "PIPELINE EXECUTION SUMMARY"
    echo "================================================================================"
    echo ""
    echo "Execution Details:"
    echo "  Timestamp: $TIMESTAMP"
    echo "  Mode: COMPLETE (14 active steps; steps 14-15 archived)"
    echo "  Total time: ${TOTAL_MINUTES}m ${TOTAL_SECONDS}s (${TOTAL_DURATION}s)"
    echo "  Log file: $LOG_FILE"
    echo ""
    echo "Environment:"
    echo "  Working directory: $SCRIPT_DIR"
    echo "  Conda environment: $CONDA_ENV"
    echo ""
    echo "Outputs Generated:"
    echo "  CSV files: $CSV_COUNT"
    echo "  Figure directories: $FIGURE_DIRS"
    echo "  Total figures: $FIGURE_COUNT"
    echo ""
    echo "================================================================================"
    echo "STEP-BY-STEP PERFORMANCE"
    echo "================================================================================"
    echo ""
    printf "%-5s %-50s %-10s %-10s\n" "Step" "Name" "Time (s)" "Status"
    echo "--------------------------------------------------------------------------------"

    for i in "${!STEP_NAMES[@]}"; do
        printf "%-5s %-50s %-10s %-10s\n" \
            "${STEP_NAMES[$i]%%:*}" \
            "${STEP_NAMES[$i]#*: }" \
            "${STEP_TIMES[$i]}" \
            "${STEP_STATUS[$i]}"
    done

    echo ""
    echo "================================================================================"
    echo "KEY FINDINGS"
    echo "================================================================================"
    echo ""

    if [ -f "analysis_imc/ERpos_survival_cox_univariate.csv" ]; then
        echo "Survival Analysis (Step 8):"
        SURVIVAL_VARS=$(wc -l < "analysis_imc/ERpos_survival_cox_univariate.csv")
        echo "  Variables tested: $((SURVIVAL_VARS - 1))"
        SIG_VARS=$(awk -F',' 'NR>1 && $5<0.05 {count++} END {print count+0}' "analysis_imc/ERpos_survival_cox_univariate.csv")
        echo "  Significant (p<0.05): $SIG_VARS"
        echo ""
    fi

    echo "Archived steps (not run):"
    echo "  - Step 14 Multivariable Cox"
    echo "  - Step 15 Hierarchical Subtyping"
    echo ""

    echo "================================================================================"
    echo "PIPELINE STATUS: $([ "$VALIDATION_OK" == true ] && echo "✅ SUCCESS" || echo "⚠ COMPLETED WITH WARNINGS")"
    echo "================================================================================"

} | tee "$SUMMARY_FILE"

log ""
log "✓ Summary saved to: $SUMMARY_FILE"

#===============================================================================
# FINAL REPORT
#===============================================================================

log_section "PIPELINE EXECUTION COMPLETE"

log "Total execution time: ${TOTAL_MINUTES}m ${TOTAL_SECONDS}s"
log "Steps completed: ${#STEP_NAMES[@]}"
log "Success rate: $(echo "${STEP_STATUS[@]}" | grep -o "SUCCESS" | wc -l)/${#STEP_STATUS[@]}"
log ""
log "Key outputs:"
log "  1. analysis_imc/*.csv ($CSV_COUNT files)"
log "  2. figures/*/ ($FIGURE_COUNT figures across $FIGURE_DIRS directories)"
log "  3. $SUMMARY_FILE (execution metrics)"
log "  4. $LOG_FILE (detailed log)"
log ""

if [ "$VALIDATION_OK" == true ]; then
    log "✅ ALL VALIDATIONS PASSED - Pipeline ready for workshop!"
else
    log "⚠ Some validations failed - Please review log file"
fi

log ""
log "Next steps:"
log "  1. Review $SUMMARY_FILE for performance metrics"
log "  2. Check figures/ for all visualizations"
log "  3. Analyze CSVs in analysis_imc/"
log "  4. Review docs/STEP*_GUIDE.md for interpretations"

log_section "END OF PIPELINE"

exit 0
