#!/usr/bin/env python
"""
Step 7: Dispersion & Isolation Metrics

ISOLATED PIPELINE - Run from pipeline_leiden_validated/ directory:
  cd /path/to/pipeline_leiden_validated
  conda run -p /path/to/squidpy-env python scripts/step7_dispersion.py

Purpose:
  Calculates spatial dispersion and isolation metrics (Wortman et al.):
  1. k-NN distances (k=5, k=10) - Measures lymphocyte dispersion
  2. Isolation metric (<5 neighbors within 20μm) - Identifies isolated cells
  3. Nearest neighbor distances (min distance to same cell type)

  As mentioned in Keren et al.:
  "Wortman et al. developed metrics for lymphocyte isolation and spatial dispersion
  which were linked to longer OS... Higher values of both metrics are indicative of
  spatial dispersion while lower values are associated with clustering."

Input:
  - data/imc_ERpos_with_embeddings_and_full_names.h5ad (with spatial coordinates)

Outputs (CSVs):
  - analysis_imc/ERpos_dispersion_metrics_by_patient.csv
  - analysis_imc/ERpos_dispersion_metrics_per_cell.csv

Outputs (Figures):
  - figures/step7_dispersion/01_dispersion_by_celltype.png
  - figures/step7_dispersion/02_isolation_by_celltype.png
  - figures/step7_dispersion/03_dispersion_vs_survival.png

Parameters:
  - K_VALUES = [5, 10] (k-nearest neighbors)
  - ISOLATION_RADIUS_UM = 20 (isolation radius in microns)
  - ISOLATION_THRESHOLD = 5 (<5 neighbors = isolated)

Notes:
  - Processes per patient/ROI to avoid cross-patient neighbors
  - Uses full cell type names (celltype_full) for better interpretability
  - Excludes cluster "20" and "Unclassified Cells" from analysis
  - Higher dispersion = more spatial mixing (linked to better survival)
"""

import os
import sys
import time
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import cKDTree
from pathlib import Path

# ============================================================================
# PARAMETERS
# ============================================================================
K_VALUES = [5, 10]  # k-nearest neighbors to compute
ISOLATION_RADIUS_UM = 20  # Isolation radius in microns
ISOLATION_THRESHOLD = 5  # <5 neighbors = isolated

# Input/Output paths
INPUT_H5AD = "data/imc_ERpos_with_embeddings_and_full_names.h5ad"
OUTPUT_DIR = Path("analysis_imc")
FIGURE_DIR = Path("figures/step7_dispersion")

# Create output directories
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
FIGURE_DIR.mkdir(exist_ok=True, parents=True)

# Set matplotlib parameters
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    print("=" * 80)
    print("STEP 7: DISPERSION & ISOLATION METRICS (Wortman et al.)")
    print("=" * 80)
    print(f"Parameters:")
    print(f"  k-NN values: {K_VALUES}")
    print(f"  Isolation: <{ISOLATION_THRESHOLD} neighbors within {ISOLATION_RADIUS_UM}μm\n")

    # ========================================================================
    # 1. Load data
    # ========================================================================
    print("[1/5] Loading ER+ dataset...")
    t0 = time.time()
    adata = sc.read_h5ad(INPUT_H5AD)

    # Use full cell type names (celltype_full column)
    if 'celltype_full' in adata.obs.columns:
        adata.obs['celltype'] = adata.obs['celltype_full'].astype('category')
    else:
        # Fallback to leiden if celltype_full not available
        adata.obs['celltype'] = adata.obs['leiden'].astype('category')

    # Filter out cluster "20" or "Unclassified" (case-insensitive)
    n_before = adata.n_obs
    celltype_str = adata.obs['celltype'].astype(str).str.lower()
    mask_unclassified = (
        celltype_str.str.contains('unclassified', na=False) |
        celltype_str.isin(['20', 'nan'])
    )
    adata = adata[~mask_unclassified].copy()
    n_after = adata.n_obs
    n_filtered = n_before - n_after

    print(f"   Loaded: {n_before:,} cells")
    print(f"   Filtered: {n_filtered:,} cells (cluster 20/Unclassified)")
    print(f"   Remaining: {n_after:,} cells × {adata.shape[1]} markers")
    print(f"   Patients: {adata.obs['Patient'].nunique()}")
    print(f"   Cell types: {adata.obs['celltype'].nunique()} unique types")
    print(f"   Time: {time.time() - t0:.1f}s")

    # Get spatial coordinates in microns
    if 'spatial' not in adata.obsm:
        print("\n   ERROR: No spatial coordinates found in adata.obsm['spatial']")
        sys.exit(1)

    coords_pixels = adata.obsm['spatial']

    # Convert to microns (1.0 microns/pixel for processed coordinates)
    MICRONS_PER_PIXEL = 1.0
    coords_um = coords_pixels * MICRONS_PER_PIXEL

    print(f"\n   Spatial coordinates range:")
    print(f"     X: {coords_um[:, 0].min():.1f} - {coords_um[:, 0].max():.1f} μm")
    print(f"     Y: {coords_um[:, 1].min():.1f} - {coords_um[:, 1].max():.1f} μm")

    # ========================================================================
    # 2. Calculate dispersion metrics
    # ========================================================================
    print("\n[2/5] Calculating dispersion metrics (k-NN distances)...")
    t0 = time.time()

    # Initialize result columns
    for k in K_VALUES:
        adata.obs[f'knn_dist_k{k}'] = 0.0

    adata.obs['isolated'] = 0
    adata.obs['nn_same_celltype_dist'] = 0.0

    # Process per patient to avoid cross-patient neighbors
    patients = adata.obs['Patient'].unique()
    total_patients = len(patients)

    for i, patient in enumerate(patients, 1):
        if i % 20 == 0 or i == total_patients:
            print(f"   Processing patient {i}/{total_patients}...", flush=True)

        patient_mask = adata.obs['Patient'].values == patient
        patient_idx = np.where(patient_mask)[0]
        patient_index = adata.obs.index[patient_idx]
        patient_coords = coords_um[patient_idx, :]

        if len(patient_idx) < max(K_VALUES) + 1:
            continue  # Skip if not enough cells

        # Build KD-tree for this patient
        tree = cKDTree(patient_coords)

        # k-NN distances
        max_k = max(K_VALUES) + 1
        distances, _ = tree.query(patient_coords, k=max_k)

        for k in K_VALUES:
            # distances[:, k] = distance to k-th neighbor (0-indexed)
            adata.obs.loc[patient_index, f'knn_dist_k{k}'] = distances[:, k]

        # Isolation metric: cells with <threshold neighbors within radius
        for local_i, global_i in enumerate(patient_idx):
            neighbors = tree.query_ball_point(patient_coords[local_i], r=ISOLATION_RADIUS_UM)
            n_neighbors = len(neighbors) - 1  # Exclude self

            if n_neighbors < ISOLATION_THRESHOLD:
                adata.obs.loc[patient_index[local_i], 'isolated'] = 1

        # Nearest neighbor distance to same celltype
        celltypes = adata.obs.loc[patient_index, 'celltype'].values

        for celltype in np.unique(celltypes):
            ct_mask_local = celltypes == celltype
            ct_idx_local = np.where(ct_mask_local)[0]

            if len(ct_idx_local) < 2:
                continue  # Skip if only 1 cell of this type

            ct_coords = patient_coords[ct_mask_local, :]
            ct_tree = cKDTree(ct_coords)

            # Query 2 nearest neighbors (self + 1 other)
            dist_ct, _ = ct_tree.query(ct_coords, k=2)

            # dist_ct[:, 1] = distance to nearest same-celltype neighbor
            global_ct_idx = patient_index[ct_idx_local]
            adata.obs.loc[global_ct_idx, 'nn_same_celltype_dist'] = dist_ct[:, 1]

    print(f"   Time: {time.time() - t0:.1f}s")

    # ========================================================================
    # 3. Aggregate to patient level
    # ========================================================================
    print("\n[3/5] Aggregating to patient level...")

    patient_metrics = []

    for patient in patients:
        patient_mask = adata.obs['Patient'] == patient
        patient_cells = adata.obs.loc[patient_mask]

        metrics = {'Patient': patient}

        # Total cells
        metrics['total_cells'] = patient_mask.sum()

        # Cell type fractions
        for celltype in adata.obs['celltype'].unique():
            n_cells = (patient_cells['celltype'] == celltype).sum()
            metrics[f'n_{celltype}'] = n_cells
            metrics[f'frac_{celltype}'] = n_cells / metrics['total_cells'] * 100 if metrics['total_cells'] > 0 else 0

        # Define immune cell types (adjust based on your leiden labels)
        IMMUNE_CELLTYPES = ['CD3+ T lymphocytes', 'CD20+ B lymphocytes', 'Macrophages']

        # k-NN dispersion metrics (mean across all cells)
        for k in K_VALUES:
            col = f'knn_dist_k{k}'
            metrics[f'mean_{col}'] = patient_cells[col].mean()
            metrics[f'median_{col}'] = patient_cells[col].median()

            # Also calculate for immune cells specifically
            immune_mask = patient_cells['celltype'].isin(IMMUNE_CELLTYPES)
            if immune_mask.sum() > 0:
                metrics[f'immune_mean_{col}'] = patient_cells.loc[immune_mask, col].mean()
                metrics[f'immune_median_{col}'] = patient_cells.loc[immune_mask, col].median()

        # Isolation metrics
        n_isolated = patient_cells['isolated'].sum()
        metrics['n_isolated'] = n_isolated
        metrics['frac_isolated'] = n_isolated / metrics['total_cells'] * 100 if metrics['total_cells'] > 0 else 0

        # Immune-specific isolation
        immune_mask = patient_cells['celltype'].isin(IMMUNE_CELLTYPES)
        if immune_mask.sum() > 0:
            immune_isolated = patient_cells.loc[immune_mask, 'isolated'].sum()
            metrics['immune_n_isolated'] = immune_isolated
            metrics['immune_frac_isolated'] = immune_isolated / immune_mask.sum() * 100

        # Nearest neighbor distances (by celltype)
        for celltype in adata.obs['celltype'].unique():
            ct_mask = patient_cells['celltype'] == celltype
            if ct_mask.sum() > 1:  # Need at least 2 cells
                metrics[f'{celltype}_nn_dist_mean'] = patient_cells.loc[ct_mask, 'nn_same_celltype_dist'].mean()
                metrics[f'{celltype}_nn_dist_median'] = patient_cells.loc[ct_mask, 'nn_same_celltype_dist'].median()

        patient_metrics.append(metrics)

    df_patient = pd.DataFrame(patient_metrics)

    # Save patient-level metrics
    output_file_patient = OUTPUT_DIR / "ERpos_dispersion_metrics_by_patient.csv"
    df_patient.to_csv(output_file_patient, index=False)
    print(f"   ✓ Saved: {output_file_patient}")
    print(f"     Shape: {df_patient.shape[0]} patients × {df_patient.shape[1]} metrics")

    # ========================================================================
    # 4. Save per-cell metrics
    # ========================================================================
    print("\n[4/5] Saving per-cell metrics...")

    # Save per-cell dispersion metrics (subset of columns)
    cell_metrics = adata.obs[[
        'Patient', 'slide_scene', 'celltype',
        'knn_dist_k5', 'knn_dist_k10',
        'isolated', 'nn_same_celltype_dist'
    ]].copy()

    # Double-check: filter out cluster 20/Unclassified (case-insensitive)
    celltype_str = cell_metrics['celltype'].astype(str).str.lower()
    mask_unclassified = (
        celltype_str.str.contains('unclassified', na=False) |
        celltype_str.isin(['20', 'nan'])
    )
    cell_metrics = cell_metrics[~mask_unclassified].copy()

    output_file_cells = OUTPUT_DIR / "ERpos_dispersion_metrics_per_cell.csv"
    cell_metrics.to_csv(output_file_cells)
    print(f"   ✓ Saved: {output_file_cells}")
    print(f"     Shape: {cell_metrics.shape[0]} cells × {cell_metrics.shape[1]} metrics")

    # Print summary statistics
    print("\n" + "=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)

    print("\nDispersion metrics (k-NN distances) - Patient averages:")
    for k in K_VALUES:
        col = f'mean_knn_dist_k{k}'
        print(f"  k={k}: mean={df_patient[col].mean():.2f}μm, median={df_patient[col].median():.2f}μm")

    print(f"\nIsolation (<{ISOLATION_THRESHOLD} neighbors @ {ISOLATION_RADIUS_UM}μm):")
    print(f"  Mean % isolated cells: {df_patient['frac_isolated'].mean():.2f}%")
    if 'immune_frac_isolated' in df_patient.columns:
        print(f"  Mean % isolated immune: {df_patient['immune_frac_isolated'].mean():.2f}%")

    # ========================================================================
    # 5. Generate visualizations
    # ========================================================================
    print("\n[5/5] Generating visualizations...")

    # ------------------------------------------------------------------------
    # Plot 1: Dispersion metrics by cell type (boxplot)
    # ------------------------------------------------------------------------
    print("   Generating Plot 1: Dispersion metrics by cell type...")

    # Prepare data for plotting (cell-level k-NN distances)
    df_plot = cell_metrics[['celltype', 'knn_dist_k5', 'knn_dist_k10']].copy()

    # Filter out cluster 20/Unclassified and rare cell types
    # Use higher threshold (5000) to ensure meaningful boxplots
    MIN_CELLS_FOR_PLOT = 5000
    celltype_counts = df_plot['celltype'].value_counts()
    valid_celltypes = celltype_counts[celltype_counts >= MIN_CELLS_FOR_PLOT].index

    celltype_str = df_plot['celltype'].astype(str).str.lower()
    mask_valid = (
        df_plot['celltype'].isin(valid_celltypes) &
        ~celltype_str.str.contains('unclassified', na=False) &
        ~celltype_str.isin(['20', 'nan'])
    )
    df_plot = df_plot[mask_valid].copy()

    # Select top 12 most abundant cell types for clarity
    top_celltypes = df_plot['celltype'].value_counts().head(12).index.tolist()
    df_plot_top = df_plot[df_plot['celltype'].isin(top_celltypes)].copy()

    # CRITICAL: Remove unused categories from categorical dtype
    # This prevents seaborn from showing empty categories
    if df_plot_top['celltype'].dtype.name == 'category':
        df_plot_top['celltype'] = df_plot_top['celltype'].cat.remove_unused_categories()

    print(f"     Filtered to {len(valid_celltypes)} cell types (>= {MIN_CELLS_FOR_PLOT} cells)")
    print(f"     Plotting top {len(top_celltypes)} most abundant cell types")
    print(f"     Actual types in plot data: {df_plot_top['celltype'].nunique()}")
    print(f"     Cell counts per type in plot:")
    for ct, count in df_plot_top['celltype'].value_counts().items():
        print(f"       - {ct}: {count}")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=300)

    # k=5
    sns.boxplot(
        data=df_plot_top,
        x='celltype',
        y='knn_dist_k5',
        palette='Set2',
        ax=axes[0]
    )
    axes[0].set_xlabel('Cell Type', fontsize=11)
    axes[0].set_ylabel('k-NN Distance (μm)', fontsize=11)
    axes[0].set_title('Dispersion: 5-NN Distance by Cell Type', fontsize=12, fontweight='bold')
    axes[0].grid(True, alpha=0.3, axis='y')
    # Rotate labels for full cell type names
    plt.setp(axes[0].xaxis.get_majorticklabels(), rotation=60, ha='right', fontsize=7)

    # k=10
    sns.boxplot(
        data=df_plot_top,
        x='celltype',
        y='knn_dist_k10',
        palette='Set2',
        ax=axes[1]
    )
    axes[1].set_xlabel('Cell Type', fontsize=11)
    axes[1].set_ylabel('k-NN Distance (μm)', fontsize=11)
    axes[1].set_title('Dispersion: 10-NN Distance by Cell Type', fontsize=12, fontweight='bold')
    axes[1].grid(True, alpha=0.3, axis='y')
    # Rotate labels for full cell type names
    plt.setp(axes[1].xaxis.get_majorticklabels(), rotation=60, ha='right', fontsize=7)

    plt.tight_layout()
    plt.savefig(FIGURE_DIR / '01_dispersion_by_celltype.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("     ✓ Saved: 01_dispersion_by_celltype.png")

    # ------------------------------------------------------------------------
    # Plot 2: Isolation metrics by cell type (barplot)
    # ------------------------------------------------------------------------
    print("   Generating Plot 2: Isolation metrics by cell type...")

    # Calculate isolation fraction per celltype (use same top_celltypes from Plot 1)
    isolation_by_ct = []

    # Use filtered cell_metrics (already excludes Unclassified and rare types)
    cell_metrics_filtered = cell_metrics[cell_metrics['celltype'].isin(top_celltypes)].copy()

    for celltype in top_celltypes:
        ct_cells = cell_metrics_filtered[cell_metrics_filtered['celltype'] == celltype]
        n_total = len(ct_cells)
        n_isolated = ct_cells['isolated'].sum()
        frac_isolated = (n_isolated / n_total * 100) if n_total > 0 else 0

        isolation_by_ct.append({
            'celltype': celltype,
            'frac_isolated': frac_isolated,
            'n_total': n_total,
            'n_isolated': n_isolated
        })

    df_isolation = pd.DataFrame(isolation_by_ct).sort_values('frac_isolated', ascending=False)

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

    bars = ax.barh(range(len(df_isolation)), df_isolation['frac_isolated'].values,
                   color='steelblue', alpha=0.8)

    ax.set_yticks(range(len(df_isolation)))
    # Use smaller font for full cell type names
    ax.set_yticklabels(df_isolation['celltype'].values, fontsize=8)
    ax.set_xlabel(f'% Isolated Cells (<{ISOLATION_THRESHOLD} neighbors @ {ISOLATION_RADIUS_UM}μm)', fontsize=11)
    ax.set_title('Cell Isolation by Cell Type', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')

    # Add value labels
    for i, (idx, row) in enumerate(df_isolation.iterrows()):
        ax.text(row['frac_isolated'] + 0.5, i, f"{row['frac_isolated']:.1f}%",
                va='center', fontsize=7)

    plt.tight_layout()
    plt.savefig(FIGURE_DIR / '02_isolation_by_celltype.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("     ✓ Saved: 02_isolation_by_celltype.png")

    # ------------------------------------------------------------------------
    # Plot 3: Dispersion vs survival (scatter)
    # ------------------------------------------------------------------------
    print("   Generating Plot 3: Dispersion vs survival...")

    # Merge dispersion metrics with survival data
    patient_obs = adata.obs[['Patient', 'Survival_time', 'Survival', 'Recurrence_time', 'Recurrence']].drop_duplicates(subset='Patient')
    df_merged = df_patient.merge(patient_obs, on='Patient', how='left')

    # Create scatter plots for k=5 and k=10 vs OS and RFS
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=300)

    # OS vs k=5
    axes[0, 0].scatter(df_merged['mean_knn_dist_k5'], df_merged['Survival_time'],
                       c=df_merged['Survival'], cmap='RdYlGn', alpha=0.6, s=30)
    axes[0, 0].set_xlabel('Mean k-NN Distance (k=5, μm)', fontsize=10)
    axes[0, 0].set_ylabel('Overall Survival Time (months)', fontsize=10)
    axes[0, 0].set_title('Dispersion vs Overall Survival (k=5)', fontsize=11, fontweight='bold')
    axes[0, 0].grid(True, alpha=0.3)

    # OS vs k=10
    axes[0, 1].scatter(df_merged['mean_knn_dist_k10'], df_merged['Survival_time'],
                       c=df_merged['Survival'], cmap='RdYlGn', alpha=0.6, s=30)
    axes[0, 1].set_xlabel('Mean k-NN Distance (k=10, μm)', fontsize=10)
    axes[0, 1].set_ylabel('Overall Survival Time (months)', fontsize=10)
    axes[0, 1].set_title('Dispersion vs Overall Survival (k=10)', fontsize=11, fontweight='bold')
    axes[0, 1].grid(True, alpha=0.3)

    # RFS vs k=5
    axes[1, 0].scatter(df_merged['mean_knn_dist_k5'], df_merged['Recurrence_time'],
                       c=df_merged['Recurrence'], cmap='RdYlGn', alpha=0.6, s=30)
    axes[1, 0].set_xlabel('Mean k-NN Distance (k=5, μm)', fontsize=10)
    axes[1, 0].set_ylabel('Recurrence-Free Survival (months)', fontsize=10)
    axes[1, 0].set_title('Dispersion vs Recurrence-Free Survival (k=5)', fontsize=11, fontweight='bold')
    axes[1, 0].grid(True, alpha=0.3)

    # RFS vs k=10
    axes[1, 1].scatter(df_merged['mean_knn_dist_k10'], df_merged['Recurrence_time'],
                       c=df_merged['Recurrence'], cmap='RdYlGn', alpha=0.6, s=30)
    axes[1, 1].set_xlabel('Mean k-NN Distance (k=10, μm)', fontsize=10)
    axes[1, 1].set_ylabel('Recurrence-Free Survival (months)', fontsize=10)
    axes[1, 1].set_title('Dispersion vs Recurrence-Free Survival (k=10)', fontsize=11, fontweight='bold')
    axes[1, 1].grid(True, alpha=0.3)

    # Add colorbar legend
    from matplotlib import cm
    from matplotlib.colors import Normalize

    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    norm = Normalize(vmin=0, vmax=1)
    sm = cm.ScalarMappable(cmap='RdYlGn', norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Event Status\n(0=alive/no recur,\n1=dead/recur)', fontsize=9)

    plt.tight_layout(rect=[0, 0, 0.9, 1])
    plt.savefig(FIGURE_DIR / '03_dispersion_vs_survival.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("     ✓ Saved: 03_dispersion_vs_survival.png")

    # ========================================================================
    # Summary
    # ========================================================================
    print("\n" + "=" * 80)
    print("STEP 7 COMPLETE")
    print("=" * 80)
    print(f"\nOutputs:")
    print(f"  CSVs:")
    print(f"    - ERpos_dispersion_metrics_by_patient.csv ({df_patient.shape[0]} patients × {df_patient.shape[1]} metrics)")
    print(f"    - ERpos_dispersion_metrics_per_cell.csv ({cell_metrics.shape[0]} cells × {cell_metrics.shape[1]} metrics)")
    print(f"  Figures:")
    print(f"    - 01_dispersion_by_celltype.png")
    print(f"    - 02_isolation_by_celltype.png")
    print(f"    - 03_dispersion_vs_survival.png")
    print("\nNext Step: Step 8 (Survival Analysis)")
    print("=" * 80)


if __name__ == "__main__":
    main()
