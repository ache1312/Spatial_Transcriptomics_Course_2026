#!/usr/bin/env python
"""
STEP 10: Tumor-Immune Mixing (neighbor-fraction index)

ISOLATED PIPELINE - Run from pipeline_leiden_validated/ directory:
  cd /path/to/pipeline_leiden_validated
  conda run -p /path/to/squidpy-env python scripts/step10_mixing_score.py

Implements a practical tumor–immune mixing index inspired by Keren et al.
(Cell, 2018), but not an exact reproduction of their contact-based definition.

This mixing index quantifies the degree of spatial intermixing between
tumor (epithelial) and immune compartments:
- High mixing: Tumor and immune cells are spatially interleaved
- Low mixing: Tumor and immune cells are segregated into compartments

Methods:
1. Local mixing index: per-cell fraction of heterotypic neighbors in a radius-based neighbor graph
2. Compartment boundaries: edge density between tumor-immune
3. Occupancy metrics: % of each compartment occupied by the other

Input:  data/imc_ERpos_with_embeddings_and_full_names.h5ad (with full cell type names)
Output: analysis_imc/ERpos_mixing_score_by_patient.csv
        analysis_imc/ERpos_step10_cox_overall_mixing_per0p1.csv
        analysis_imc/ERpos_step10_cox_celltype_mixing_per0p1.csv
        analysis_imc/ERpos_step10_cox_overall_mixing_stratified_by_prolif.csv
        figures/step10_mixing_score/*.png (01–05 + spatial examples)
"""

import os
import time
import warnings
import pandas as pd
import numpy as np
import scanpy as sc
import squidpy as sq
from scipy.spatial import Voronoi, voronoi_plot_2d  # pylint: disable=no-name-in-module
from pathlib import Path

warnings.filterwarnings('ignore')

# Settings
RADIUS_UM = 40  # Neighborhood radius
MICRONS_PER_PIXEL = 1.0

OUTPUT_DIR = Path("analysis_imc")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

print("=" * 80)
print("STEP 10: TUMOR-IMMUNE MIXING SCORE (Keren et al. Fig 4I)")
print("=" * 80)
print(f"Radius: {RADIUS_UM}μm")
print()

# Load data
print("[1/4] Loading ER+ dataset...")
t0 = time.time()
adata = sc.read_h5ad("data/imc_ERpos_with_embeddings_and_full_names.h5ad")

# Use full cell type names (celltype_full column)
if 'celltype_full' in adata.obs.columns:
    adata.obs['celltype'] = adata.obs['celltype_full'].astype('category')
    # Exclude "Unclassified cells"
    adata = adata[adata.obs['celltype'] != 'Unclassified cells'].copy()
else:
    raise ValueError("Column 'celltype_full' not found in h5ad file")

print(f"Loaded: {adata.shape[0]:,} cells × {adata.shape[1]} markers")
print(f"Patients: {adata.obs['Patient'].nunique()}")
print(f"Cell types: {adata.obs['celltype'].nunique()}")
print(f"Time: {time.time() - t0:.1f}s")

# Get spatial coordinates
coords_pixels = adata.obsm['spatial']
coords_um = coords_pixels * MICRONS_PER_PIXEL

# Simplify to tumor (epithelial) vs immune
# Map full cell type names to compartments
adata.obs['compartment'] = 'other'

# Tumor: All epithelial/tumor phenotypes (full names)
tumor_types = [
    'CD44+ tumor cells',
    'Basal tumor cells',
    'Luminal tumor cells',
    'Luminal ER+ tumor cells',
    'Cytokeratin-low tumor cells',
    'Proliferating tumor cells',
    'Ki67+ proliferative tumor cells',
    'HER2+ Ki67+ proliferative tumor cells',
    'HER2+ tumor cells',
    'HER2+ ER+ tumor cells',
    'HER2+ ER+ PR+ tumor cells',
    'HER2+ Ki67+ tumor cells',
    'EGFR+ tumor cells',
    'Myoepithelial cells'
]
adata.obs.loc[adata.obs['celltype'].isin(tumor_types), 'compartment'] = 'tumor'

# Immune: T cells, B cells, Macrophages (full names)
immune_types = [
    'CD3+ T lymphocytes',
    'CD20+ B lymphocytes',
    'CD68+ macrophages'
]
adata.obs.loc[adata.obs['celltype'].isin(immune_types), 'compartment'] = 'immune'

print(f"\nCompartments:")
print(f"  Tumor (epithelial): {(adata.obs['compartment'] == 'tumor').sum():,} cells")
print(f"  Immune: {(adata.obs['compartment'] == 'immune').sum():,} cells")
print(f"  Other (stromal, endothelial): {(adata.obs['compartment'] == 'other').sum():,} cells")


def calculate_mixing_index(adata_roi, radius=40):
    """
    Calculate local mixing index for tumor-immune

    Mixing index = (heterotypic neighbors) / (total neighbors)

    High mixing: Cells have many neighbors of opposite type
    Low mixing: Cells have few neighbors of opposite type

    Args:
        adata_roi: AnnData for single ROI
        radius: Neighborhood radius (microns)

    Returns:
        dict with mixing metrics
    """
    def _mean_neighbor_fraction(conn_mat, src_mask, tgt_mask):
        src_idx = np.where(src_mask)[0]
        if src_idx.size == 0:
            return 0.0, 0
        sub = conn_mat[src_idx]
        total = np.asarray(sub.sum(axis=1)).ravel()
        cross = np.asarray(sub[:, tgt_mask].sum(axis=1)).ravel()
        valid = total > 0
        if not valid.any():
            return 0.0, int(valid.sum())
        return float(np.mean(cross[valid] / total[valid])), int(valid.sum())

    immune_subtypes = {
        'cd3': 'CD3+ T lymphocytes',
        'cd68': 'CD68+ macrophages',
        'cd20': 'CD20+ B lymphocytes',
    }

    # Build spatial graph
    sq.gr.spatial_neighbors(
        adata_roi,
        coord_type="generic",
        radius=radius,
        delaunay=False
    )

    conn = adata_roi.obsp['spatial_connectivities']
    compartments = adata_roi.obs['compartment'].values
    celltypes = adata_roi.obs['celltype'].astype(str).values

    # For tumor cells, count immune neighbors
    tumor_mask = compartments == 'tumor'
    immune_mask = compartments == 'immune'

    if tumor_mask.sum() == 0 or immune_mask.sum() == 0:
        out = {
            'tumor_immune_mixing': 0,
            'immune_tumor_mixing': 0,
            'overall_mixing': 0,
            'n_tumor': tumor_mask.sum(),
            'n_immune': immune_mask.sum(),
            'n_tumor_with_immune_neighbors': 0,
            'n_immune_with_tumor_neighbors': 0,
        }
        for key, label in immune_subtypes.items():
            subtype_mask = celltypes == label
            out[f'n_{key}'] = int(subtype_mask.sum())
            out[f'tumor_{key}_mix'] = 0.0
            out[f'{key}_tumor_mix'] = 0.0
            out[f'overall_{key}_mix'] = 0.0
        return out

    tumor_mixing_score, n_tumor_with_neighbors = _mean_neighbor_fraction(conn, tumor_mask, immune_mask)
    immune_mixing_score, n_immune_with_neighbors = _mean_neighbor_fraction(conn, immune_mask, tumor_mask)
    overall_mixing = (tumor_mixing_score + immune_mixing_score) / 2

    out = {
        'tumor_immune_mixing': tumor_mixing_score,
        'immune_tumor_mixing': immune_mixing_score,
        'overall_mixing': overall_mixing,
        'n_tumor': tumor_mask.sum(),
        'n_immune': immune_mask.sum(),
        'n_tumor_with_immune_neighbors': n_tumor_with_neighbors,
        'n_immune_with_tumor_neighbors': n_immune_with_neighbors,
    }

    # Cell-type-specific mixing (Tumor × CD3/CD68/CD20)
    for key, label in immune_subtypes.items():
        subtype_mask = celltypes == label
        out[f'n_{key}'] = int(subtype_mask.sum())

        t2s, _n_t_with, = _mean_neighbor_fraction(conn, tumor_mask, subtype_mask)
        s2t, _n_s_with, = _mean_neighbor_fraction(conn, subtype_mask, tumor_mask)
        out[f'tumor_{key}_mix'] = t2s
        out[f'{key}_tumor_mix'] = s2t
        out[f'overall_{key}_mix'] = (t2s + s2t) / 2.0

    return out


def calculate_boundary_density(adata_roi, radius=40):
    """
    Calculate density of tumor-immune boundaries

    Boundary density = # tumor-immune edges / max possible edges

    High boundary density: Many contact points between compartments
    Low boundary density: Few contact points (segregation)

    Args:
        adata_roi: AnnData for single ROI
        radius: Neighborhood radius (microns)

    Returns:
        float: boundary density [0, 1]
    """
    if 'spatial_connectivities' not in adata_roi.obsp:
        sq.gr.spatial_neighbors(
            adata_roi,
            coord_type="generic",
            radius=radius,
            delaunay=False
        )

    conn = adata_roi.obsp['spatial_connectivities']
    compartments = adata_roi.obs['compartment'].values

    tumor_mask = compartments == 'tumor'
    immune_mask = compartments == 'immune'

    n_tumor = tumor_mask.sum()
    n_immune = immune_mask.sum()

    if n_tumor == 0 or n_immune == 0:
        return 0

    # Count tumor-immune edges
    n_boundary_edges = 0

    for idx in np.where(tumor_mask)[0]:
        neighbors_idx = conn[idx].nonzero()[1]
        n_boundary_edges += immune_mask[neighbors_idx].sum()

    # Max possible edges
    max_edges = n_tumor * n_immune

    # Normalize by max
    boundary_density = n_boundary_edges / max_edges if max_edges > 0 else 0

    return boundary_density


def calculate_occupancy_metrics(adata_roi):
    """
    Calculate compartment occupancy metrics

    Measures spatial extent and overlap of compartments using
    Voronoi tessellation or simple bounding boxes

    Args:
        adata_roi: AnnData for single ROI

    Returns:
        dict with occupancy metrics
    """
    coords = adata_roi.obsm['spatial'] * MICRONS_PER_PIXEL
    compartments = adata_roi.obs['compartment'].values

    tumor_mask = compartments == 'tumor'
    immune_mask = compartments == 'immune'

    if tumor_mask.sum() == 0 or immune_mask.sum() == 0:
        return {
            'tumor_area_frac': 0,
            'immune_area_frac': 0,
            'overlap_area_frac': 0
        }

    # Simple bounding box approach
    tumor_coords = coords[tumor_mask, :]
    immune_coords = coords[immune_mask, :]

    # Bounding boxes
    tumor_bbox = [
        tumor_coords[:, 0].min(), tumor_coords[:, 0].max(),
        tumor_coords[:, 1].min(), tumor_coords[:, 1].max()
    ]

    immune_bbox = [
        immune_coords[:, 0].min(), immune_coords[:, 0].max(),
        immune_coords[:, 1].min(), immune_coords[:, 1].max()
    ]

    # Total tissue area (all cells)
    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()
    total_area = (x_max - x_min) * (y_max - y_min)

    if total_area == 0:
        return {
            'tumor_area_frac': 0,
            'immune_area_frac': 0,
            'overlap_area_frac': 0
        }

    # Areas
    tumor_area = (tumor_bbox[1] - tumor_bbox[0]) * (tumor_bbox[3] - tumor_bbox[2])
    immune_area = (immune_bbox[1] - immune_bbox[0]) * (immune_bbox[3] - immune_bbox[2])

    # Overlap (intersection of bboxes)
    x_overlap = max(0, min(tumor_bbox[1], immune_bbox[1]) - max(tumor_bbox[0], immune_bbox[0]))
    y_overlap = max(0, min(tumor_bbox[3], immune_bbox[3]) - max(tumor_bbox[2], immune_bbox[2]))
    overlap_area = x_overlap * y_overlap

    return {
        'tumor_area_frac': tumor_area / total_area,
        'immune_area_frac': immune_area / total_area,
        'overlap_area_frac': overlap_area / total_area
    }


print("\n[2/4] Calculating mixing metrics per ROI...")
t0 = time.time()

scenes = adata.obs['slide_scene'].unique()
total = len(scenes)

all_metrics = []

for i, scene in enumerate(scenes, 1):
    if i % 50 == 0 or i == total:
        print(f"  Processing ROI {i}/{total}...", flush=True)

    scene_mask = adata.obs['slide_scene'] == scene
    adata_scene = adata[scene_mask].copy()

    if adata_scene.n_obs < 10:
        continue

    patient = adata_scene.obs['Patient'].iloc[0]

    # Calculate mixing metrics
    mixing = calculate_mixing_index(adata_scene, radius=RADIUS_UM)
    boundary_dens = calculate_boundary_density(adata_scene, radius=RADIUS_UM)
    occupancy = calculate_occupancy_metrics(adata_scene)

    metrics = {
        'slide_scene': scene,
        'Patient': patient,
        'n_cells': adata_scene.n_obs,
        **mixing,
        'boundary_density': boundary_dens,
        **occupancy
    }

    all_metrics.append(metrics)

print(f"Time: {time.time() - t0:.1f}s")

print("\n[3/4] Aggregating to patient level...")
df = pd.DataFrame(all_metrics)

# Aggregate to patient
patient_metrics = df.groupby('Patient').agg({
    'n_cells': 'sum',
    'n_tumor': 'sum',
    'n_immune': 'sum',
    'n_cd3': 'sum',
    'n_cd68': 'sum',
    'n_cd20': 'sum',
    'tumor_immune_mixing': 'mean',
    'immune_tumor_mixing': 'mean',
    'overall_mixing': 'mean',
    'tumor_cd3_mix': 'mean',
    'cd3_tumor_mix': 'mean',
    'overall_cd3_mix': 'mean',
    'tumor_cd68_mix': 'mean',
    'cd68_tumor_mix': 'mean',
    'overall_cd68_mix': 'mean',
    'tumor_cd20_mix': 'mean',
    'cd20_tumor_mix': 'mean',
    'overall_cd20_mix': 'mean',
    'boundary_density': 'mean',
    'tumor_area_frac': 'mean',
    'immune_area_frac': 'mean',
    'overlap_area_frac': 'mean',
    'n_tumor_with_immune_neighbors': 'sum',
    'n_immune_with_tumor_neighbors': 'sum'
}).reset_index()

# Calculate derived metrics
patient_metrics['tumor_immune_ratio'] = patient_metrics['n_tumor'] / patient_metrics['n_immune']
patient_metrics['frac_tumor'] = patient_metrics['n_tumor'] / patient_metrics['n_cells'] * 100
patient_metrics['frac_immune'] = patient_metrics['n_immune'] / patient_metrics['n_cells'] * 100

print("\n[4/4] Saving results...")
output_file = OUTPUT_DIR / "ERpos_mixing_score_by_patient.csv"
patient_metrics.to_csv(output_file, index=False)
print(f"Saved: {output_file}")
print(f"Shape: {patient_metrics.shape[0]} patients × {patient_metrics.shape[1]} metrics")

print("\n" + "=" * 80)
print("SUMMARY STATISTICS")
print("=" * 80)

print("\nTumor-Immune Mixing Scores (Patient Averages):")
print(f"  Overall mixing score: {patient_metrics['overall_mixing'].mean():.3f} ± {patient_metrics['overall_mixing'].std():.3f}")
print(f"  Tumor-to-immune mixing: {patient_metrics['tumor_immune_mixing'].mean():.3f}")
print(f"  Immune-to-tumor mixing: {patient_metrics['immune_tumor_mixing'].mean():.3f}")
print(f"  Boundary density: {patient_metrics['boundary_density'].mean():.3f}")

print("\nCompartment Occupancy:")
print(f"  Tumor area fraction: {patient_metrics['tumor_area_frac'].mean():.2f}")
print(f"  Immune area fraction: {patient_metrics['immune_area_frac'].mean():.2f}")
print(f"  Overlap area fraction: {patient_metrics['overlap_area_frac'].mean():.2f}")

print("\nInterpretation:")
print("  - Mixing score > 0.5: High intermixing of tumor and immune")
print("  - Mixing score < 0.3: Segregated compartments")
print("  - Boundary density: Quantifies contact points")
print("  - Overlap area: Spatial co-localization")

print("\nNote on original definition (Keren et al., Cell 2018):")
print("  The paper defines neighbors via direct cell contact and uses an interaction-ratio mixing score.")
print("  This script uses a radius-based neighbor graph and a neighbor-fraction mixing index for IMC.")

# ============================================================================
# VISUALIZATIONS
# ============================================================================
print("\n[GENERATING VISUALIZATIONS]")

import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

# Set matplotlib parameters
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

# Create figures directory
FIGURES_DIR = Path("figures/step10_mixing_score")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ========================================================================
# PLOT 1: Mixing score distribution by patient
# ========================================================================
print("Generating Plot 1: Mixing score distribution...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=300)

# Plot 1a: Overall mixing score
ax = axes[0, 0]
sns.histplot(
    data=patient_metrics,
    x='overall_mixing',
    bins=30,
    kde=True,
    color='steelblue',
    ax=ax
)
ax.axvline(patient_metrics['overall_mixing'].median(), color='red',
           linestyle='--', linewidth=1.5, label=f"Median = {patient_metrics['overall_mixing'].median():.3f}")
ax.set_xlabel('Overall Mixing Score', fontsize=11)
ax.set_ylabel('Number of Patients', fontsize=11)
ax.set_title('Overall Tumor-Immune Mixing Score Distribution', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Plot 1b: Tumor-to-immune mixing
ax = axes[0, 1]
sns.histplot(
    data=patient_metrics,
    x='tumor_immune_mixing',
    bins=30,
    kde=True,
    color='coral',
    ax=ax
)
ax.axvline(patient_metrics['tumor_immune_mixing'].median(), color='red',
           linestyle='--', linewidth=1.5, label=f"Median = {patient_metrics['tumor_immune_mixing'].median():.3f}")
ax.set_xlabel('Tumor→Immune Mixing Score', fontsize=11)
ax.set_ylabel('Number of Patients', fontsize=11)
ax.set_title('Tumor Cell Mixing (Immune Neighbors)', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Plot 1c: Boundary density
ax = axes[1, 0]
sns.histplot(
    data=patient_metrics,
    x='boundary_density',
    bins=30,
    kde=True,
    color='seagreen',
    ax=ax
)
ax.axvline(patient_metrics['boundary_density'].median(), color='red',
           linestyle='--', linewidth=1.5, label=f"Median = {patient_metrics['boundary_density'].median():.4f}")
ax.set_xlabel('Boundary Density', fontsize=11)
ax.set_ylabel('Number of Patients', fontsize=11)
ax.set_title('Tumor-Immune Boundary Density', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Plot 1d: Overlap area fraction
ax = axes[1, 1]
sns.histplot(
    data=patient_metrics,
    x='overlap_area_frac',
    bins=30,
    kde=True,
    color='mediumpurple',
    ax=ax
)
ax.axvline(patient_metrics['overlap_area_frac'].median(), color='red',
           linestyle='--', linewidth=1.5, label=f"Median = {patient_metrics['overlap_area_frac'].median():.3f}")
ax.set_xlabel('Overlap Area Fraction', fontsize=11)
ax.set_ylabel('Number of Patients', fontsize=11)
ax.set_title('Tumor-Immune Spatial Overlap', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

plt.suptitle('Tumor-Immune Mixing Metrics Distribution', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig(FIGURES_DIR / '01_mixing_distribution.png', dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: 01_mixing_distribution.png")

# ========================================================================
# PLOT 2: Mixing score vs survival (KM curves by quartile)
# ========================================================================
print("Generating Plot 2: Mixing score vs survival...")

# Merge with survival + clinical covariates
patient_obs = adata.obs[
    ['Patient', 'Survival_time', 'Survival', 'Recurrence_time', 'Recurrence', 'age', 'grade', 'Stage', 'N_Stage', 'cohort']
].drop_duplicates(subset='Patient')
df_merged = patient_metrics.merge(patient_obs, on='Patient', how='left')


def _to_months(series):
    s = pd.to_numeric(series, errors='coerce')
    if s.dropna().empty:
        return s
    return s / 30.44 if s.max() > 500 else s


# Normalize times (these are in days in this dataset; convert to months for plotting/modeling)
df_merged['OS_months'] = _to_months(df_merged['Survival_time'])
df_merged['OS_event'] = pd.to_numeric(df_merged['Survival'], errors='coerce')
df_merged['RFS_months'] = _to_months(df_merged['Recurrence_time'])
df_merged['RFS_event'] = pd.to_numeric(df_merged['Recurrence'], errors='coerce')

# Remove missing survival data
df_merged_os = df_merged.dropna(subset=['OS_months', 'OS_event'])
df_merged_rfs = df_merged.dropna(subset=['RFS_months', 'RFS_event'])

if len(df_merged_os) > 0:
    # Stratify by mixing score quartiles (rank-based to avoid qcut edge duplicates)
    quartile_labels = ['Q1 (Low)', 'Q2', 'Q3', 'Q4 (High)']
    df_merged_os['mixing_quartile'] = pd.qcut(
        df_merged_os['overall_mixing'].rank(method='first'),
        q=4,
        labels=quartile_labels
    )

    if len(df_merged_rfs) > 0:
        df_merged_rfs['mixing_quartile'] = pd.qcut(
            df_merged_rfs['overall_mixing'].rank(method='first'),
            q=4,
            labels=quartile_labels
        )
    else:
        df_merged_rfs['mixing_quartile'] = pd.Series(index=df_merged_rfs.index, dtype='object')

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=300)

    # Plot 2a: Overall Survival
    ax = axes[0]
    kmf = KaplanMeierFitter()

    for quartile in quartile_labels:
        mask = df_merged_os['mixing_quartile'] == quartile
        if mask.sum() > 0:
            kmf.fit(
                df_merged_os.loc[mask, 'OS_months'],
                df_merged_os.loc[mask, 'OS_event'],
                label=f'{quartile} (n={mask.sum()})'
            )
            kmf.plot_survival_function(ax=ax, ci_show=False)

    ax.set_xlabel('Time (months)', fontsize=11)
    ax.set_ylabel('Survival Probability', fontsize=11)
    ax.set_title('Overall Survival by Mixing Score Quartile', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9, loc='best')
    ax.grid(True, alpha=0.3)

    try:
        os_p = multivariate_logrank_test(
            event_durations=df_merged_os['OS_months'],
            groups=df_merged_os['mixing_quartile'],
            event_observed=df_merged_os['OS_event']
        ).p_value
        ax.text(
            0.02, 0.02,
            f'Log-rank p = {os_p:.3g}',
            transform=ax.transAxes,
            fontsize=10,
            ha='left',
            va='bottom',
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white', alpha=0.75, edgecolor='#cbd5e1'),
        )
    except Exception:  # noqa: BLE001
        pass

    # Plot 2b: Recurrence-Free Survival
    ax = axes[1]
    kmf = KaplanMeierFitter()

    if len(df_merged_rfs) > 0:
        for quartile in quartile_labels:
            mask = df_merged_rfs['mixing_quartile'] == quartile
            if mask.sum() > 0:
                kmf.fit(
                    df_merged_rfs.loc[mask, 'RFS_months'],
                    df_merged_rfs.loc[mask, 'RFS_event'],
                    label=f'{quartile} (n={mask.sum()})'
                )
                kmf.plot_survival_function(ax=ax, ci_show=False)

        try:
            rfs_p = multivariate_logrank_test(
                event_durations=df_merged_rfs['RFS_months'],
                groups=df_merged_rfs['mixing_quartile'],
                event_observed=df_merged_rfs['RFS_event']
            ).p_value
            ax.text(
                0.02, 0.02,
                f'Log-rank p = {rfs_p:.3g}',
                transform=ax.transAxes,
                fontsize=10,
                ha='left',
                va='bottom',
                bbox=dict(boxstyle='round,pad=0.25', facecolor='white', alpha=0.75, edgecolor='#cbd5e1'),
            )
        except Exception:  # noqa: BLE001
            pass
    else:
        ax.text(0.5, 0.5, 'RFS unavailable', transform=ax.transAxes, ha='center', va='center')

    ax.set_xlabel('Time (months)', fontsize=11)
    ax.set_ylabel('Recurrence-Free Probability', fontsize=11)
    ax.set_title('Recurrence-Free Survival by Mixing Score Quartile', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9, loc='best')
    ax.grid(True, alpha=0.3)

    plt.suptitle('Survival Analysis: Tumor-Immune Mixing Score (Quartiles)', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / '02_mixing_vs_survival.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: 02_mixing_vs_survival.png")
else:
    print("  Skipped: No survival data available")

# ========================================================================
# PLOT 2b: Multivariable Cox models (overall + cell-type-specific mixing)
# ========================================================================
print("Generating Plot 2b: Multivariable Cox models...")

from lifelines import CoxPHFitter


def _zscore(series):
    s = pd.to_numeric(series, errors='coerce')
    mu = s.mean()
    sd = s.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return s * np.nan
    return (s - mu) / sd


def _fit_cox_table(df_in, duration_col, event_col, covariates, penalizer=0.1, label=''):
    model_df = df_in[[duration_col, event_col, *covariates]].copy()
    model_df = model_df.dropna()
    model_df = model_df[model_df[duration_col] > 0].copy()
    model_df[event_col] = pd.to_numeric(model_df[event_col], errors='coerce').astype(int)

    # Drop constant covariates in this subset
    kept = []
    dropped = []
    for c in covariates:
        if model_df[c].nunique(dropna=True) <= 1:
            dropped.append(c)
            continue
        kept.append(c)

    if dropped:
        print(f"  Dropped non-informative covariates ({label}): {dropped}")

    model_df = model_df[[duration_col, event_col, *kept]].copy()

    # Robust fitting: retry with stronger ridge penalty if convergence fails.
    last_err = None
    cph = None
    for pen in [penalizer, max(1.0, penalizer * 10), max(10.0, penalizer * 100)]:
        try:
            cph = CoxPHFitter(penalizer=pen, l1_ratio=0.0)
            cph.fit(model_df, duration_col=duration_col, event_col=event_col)
            break
        except Exception as e:  # noqa: BLE001
            last_err = e
            cph = None
            continue

    if cph is None:
        print(f"  Cox fit failed ({label}): {type(last_err).__name__}: {last_err}")
        return pd.DataFrame()

    summ = cph.summary.reset_index()
    if 'index' in summ.columns:
        summ = summ.rename(columns={'index': 'covariate'})
    elif 'covariate' not in summ.columns:
        summ = summ.rename(columns={summ.columns[0]: 'covariate'})

    summ['HR'] = np.exp(summ['coef'])
    summ['HR_lower_95'] = np.exp(summ['coef lower 95%'])
    summ['HR_upper_95'] = np.exp(summ['coef upper 95%'])
    summ['n'] = len(model_df)
    return summ


def _forest_plot_overall(mix_rows, outpath):
    if mix_rows.empty:
        return

    # Step 8-style forest plot (point + CI + annotation boxes)
    fig, ax = plt.subplots(figsize=(10, max(3.2, len(mix_rows) * 0.55)), dpi=250)
    df_plot = mix_rows.copy().reset_index(drop=True)
    df_plot['var_display'] = df_plot['outcome'].astype(str)

    y_pos = np.arange(len(df_plot))
    colors = ['red' if hr > 1 else 'blue' for hr in df_plot['HR']]
    ax.scatter(df_plot['HR'], y_pos, s=90, c=colors, zorder=3, alpha=0.8)

    for idx, row in df_plot.iterrows():
        ax.plot([row['HR_lower_95'], row['HR_upper_95']], [idx, idx], c=colors[idx], linewidth=2, alpha=0.6, zorder=2)

    ax.axvline(x=1.0, color='black', linestyle='--', linewidth=1, alpha=0.5, zorder=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df_plot['var_display'], fontsize=11)
    ax.set_xlabel('Hazard Ratio (95% CI)', fontsize=12, fontweight='bold')
    ax.set_title('Cox Model (Adjusted): Overall Mixing per +0.1', fontsize=14, fontweight='bold', pad=16)
    ax.grid(axis='x', alpha=0.3, zorder=0)

    for idx, row in df_plot.iterrows():
        sig_marker = " *" if row['p'] < 0.05 else ""
        hr_text = f"HR={row['HR']:.2f}\np={row['p']:.4f}{sig_marker}"
        x_pos = max(row['HR_upper_95'], ax.get_xlim()[1] * 0.7)
        bg_color = 'lightyellow' if row['p'] < 0.05 else 'white'
        edge_color = 'orange' if row['p'] < 0.05 else 'none'
        ax.text(
            x_pos,
            idx,
            hr_text,
            va='center',
            fontsize=9,
            bbox=dict(boxstyle='round,pad=0.3', facecolor=bg_color, alpha=0.8, edgecolor=edge_color, linewidth=1.5),
        )

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='red', alpha=0.8, label='Increased risk (HR > 1)'),
        Patch(facecolor='blue', alpha=0.8, label='Decreased risk (HR < 1)'),
        Patch(facecolor='lightyellow', alpha=0.8, edgecolor='orange', linewidth=1.5, label='Statistically significant (p < 0.05)'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=10)

    ax.text(
        0.02,
        -0.22,
        "Adjusted for: age, grade, stage, N stage, cohort (when varying), proliferation proxy (Ki67-like)\n"
        "Ridge penalizer=0.1; HR per +0.1 increase in overall mixing.",
        transform=ax.transAxes,
        fontsize=10,
        ha='left',
        va='top',
    )

    fig.tight_layout()
    fig.savefig(outpath, bbox_inches='tight')
    plt.close(fig)

def _forest_plot_multivariable_full(os_summ, rfs_summ, covar_order, outpath):
    """
    Step 8-style forest plot for the *full* adjusted Cox model terms (mixing + covariates),
    shown separately for OS and RFS.
    """
    def _display_name(cov):
        base = {
            'overall_mixing_per_0p1': 'Overall mixing',
            'age_z': 'Age',
            'grade_z': 'Grade',
            'stage_z': 'Stage',
            'n_stage_z': 'N stage',
            'prolif_frac_z': 'Proliferation proxy',
        }
        if cov in base:
            return base[cov]
        if cov.startswith('cohort_'):
            return f"Cohort: {cov.replace('cohort_', '')}"
        return cov

    outcomes = [('OS', os_summ), ('RFS', rfs_summ)]
    dfs = []
    for outcome, summ in outcomes:
        if summ is None or summ.empty:
            dfs.append(pd.DataFrame())
            continue
        summ = summ.copy()
        summ['covariate'] = summ['covariate'].astype(str)
        keep = [c for c in covar_order if c in set(summ['covariate'])]
        summ = summ[summ['covariate'].isin(keep)].copy()
        summ['order'] = summ['covariate'].map({c: i for i, c in enumerate(keep)})
        summ = summ.sort_values('order')
        summ['var_display'] = summ['covariate'].map(_display_name)
        summ['outcome'] = outcome
        dfs.append(summ)

    if all(d.empty for d in dfs):
        return

    # Determine common x-limits
    all_lo = np.concatenate([d['HR_lower_95'].to_numpy() for d in dfs if not d.empty])
    all_hi = np.concatenate([d['HR_upper_95'].to_numpy() for d in dfs if not d.empty])
    xmin = max(0.05, float(np.nanmin(all_lo)) * 0.9)
    xmax = max(1.2, float(np.nanmax(all_hi)) * 1.1)

    heights = [max(2.6, len(d) * 0.55) if not d.empty else 2.0 for d in dfs]
    fig, axes = plt.subplots(2, 1, figsize=(10.5, sum(heights)), dpi=250, sharex=True)

    for ax, (outcome, df_plot) in zip(axes, zip(['OS', 'RFS'], dfs)):
        if df_plot is None or df_plot.empty:
            ax.axis('off')
            continue
        df_plot = df_plot.reset_index(drop=True)
        y_pos = np.arange(len(df_plot))
        colors = ['red' if hr > 1 else 'blue' for hr in df_plot['HR']]
        ax.scatter(df_plot['HR'], y_pos, s=90, c=colors, zorder=3, alpha=0.8)
        for idx, row in df_plot.iterrows():
            ax.plot([row['HR_lower_95'], row['HR_upper_95']], [idx, idx], c=colors[idx], linewidth=2, alpha=0.6, zorder=2)
        ax.axvline(x=1.0, color='black', linestyle='--', linewidth=1, alpha=0.5, zorder=1)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(df_plot['var_display'], fontsize=11)
        ax.set_title(outcome, fontsize=12, fontweight='bold', pad=10)
        ax.grid(axis='x', alpha=0.3, zorder=0)
        ax.set_xlim(xmin, xmax)

        for idx, row in df_plot.iterrows():
            p = float(row['p'])
            sig_marker = " *" if p < 0.05 else ""
            hr_text = f"HR={row['HR']:.2f}\np={p:.4f}{sig_marker}"
            x_pos = max(float(row['HR_upper_95']), ax.get_xlim()[1] * 0.7)
            bg_color = 'lightyellow' if p < 0.05 else 'white'
            edge_color = 'orange' if p < 0.05 else 'none'
            ax.text(
                x_pos,
                idx,
                hr_text,
                va='center',
                fontsize=9,
                bbox=dict(boxstyle='round,pad=0.3', facecolor=bg_color, alpha=0.8, edgecolor=edge_color, linewidth=1.5),
            )

    axes[-1].set_xlabel('Hazard Ratio (95% CI)', fontsize=12, fontweight='bold')
    fig.suptitle('Adjusted Cox Model: Overall Mixing and Covariates', fontsize=14, fontweight='bold', y=0.98)

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='red', alpha=0.8, label='Increased risk (HR > 1)'),
        Patch(facecolor='blue', alpha=0.8, label='Decreased risk (HR < 1)'),
        Patch(facecolor='lightyellow', alpha=0.8, edgecolor='orange', linewidth=1.5, label='Statistically significant (p < 0.05)'),
    ]
    # Legend in the upper-left, without shrinking the axes excessively.
    fig.legend(
        handles=legend_elements,
        loc='upper left',
        bbox_to_anchor=(0.02, 0.98),
        fontsize=8,
        frameon=True,
        borderpad=0.4,
        labelspacing=0.3,
        handlelength=1.4,
        handletextpad=0.6,
    )
    fig.tight_layout(rect=[0, 0.02, 1, 0.90])
    fig.savefig(outpath, bbox_inches='tight')
    plt.close(fig)


def _forest_plot_celltypes(df_rows, outpath):
    if df_rows.empty:
        return
    order = ['mix_cd3', 'mix_cd68', 'mix_cd20']
    label_map = {'mix_cd3': 'Tumor×CD3', 'mix_cd68': 'Tumor×CD68', 'mix_cd20': 'Tumor×CD20'}

    # Step 8-style forest plot with two panels (OS/RFS)
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 7.0), dpi=250, sharex=True)
    for ax, outcome in zip(axes, ['OS', 'RFS']):
        sub = df_rows[df_rows['outcome'] == outcome].copy()
        if sub.empty:
            ax.axis('off')
            continue
        sub['term'] = pd.Categorical(sub['term'], categories=order, ordered=True)
        sub = sub.sort_values('term')

        y_pos = np.arange(len(sub))
        hr = sub['HR'].to_numpy()
        lo = sub['HR_lower_95'].to_numpy()
        hi = sub['HR_upper_95'].to_numpy()
        pvals = sub['p'].to_numpy()

        colors = ['red' if h > 1 else 'blue' for h in hr]
        ax.scatter(hr, y_pos, s=90, c=colors, zorder=3, alpha=0.8)
        for idx, (l, u) in enumerate(zip(lo, hi)):
            ax.plot([l, u], [idx, idx], c=colors[idx], linewidth=2, alpha=0.6, zorder=2)

        ax.axvline(x=1.0, color='black', linestyle='--', linewidth=1, alpha=0.5, zorder=1)
        ax.set_yticks(y_pos)
        ax.set_yticklabels([label_map.get(t, t) for t in sub['term'].astype(str)], fontsize=11)
        ax.set_title(outcome, fontsize=12, fontweight='bold', pad=10)
        ax.grid(axis='x', alpha=0.3, zorder=0)

        for idx, (h, l, u, p) in enumerate(zip(hr, lo, hi, pvals)):
            sig_marker = " *" if p < 0.05 else ""
            hr_text = f"HR={h:.2f}\np={p:.4f}{sig_marker}"
            x_pos = max(u, ax.get_xlim()[1] * 0.7)
            bg_color = 'lightyellow' if p < 0.05 else 'white'
            edge_color = 'orange' if p < 0.05 else 'none'
            ax.text(
                x_pos,
                idx,
                hr_text,
                va='center',
                fontsize=9,
                bbox=dict(boxstyle='round,pad=0.3', facecolor=bg_color, alpha=0.8, edgecolor=edge_color, linewidth=1.5),
            )

    axes[-1].set_xlabel('Hazard Ratio (95% CI)', fontsize=12, fontweight='bold')
    xmax = float(np.nanmax(df_rows['HR_upper_95']))
    xmin = float(np.nanmin(df_rows['HR_lower_95']))
    axes[-1].set_xlim(max(0.05, xmin * 0.9), max(1.2, xmax * 1.1))
    fig.suptitle('Adjusted Cox Model: Cell-type-specific Mixing', fontsize=14, fontweight='bold', y=0.98)
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='red', alpha=0.8, label='Increased risk (HR > 1)'),
        Patch(facecolor='blue', alpha=0.8, label='Decreased risk (HR < 1)'),
        Patch(facecolor='lightyellow', alpha=0.8, edgecolor='orange', linewidth=1.5, label='Statistically significant (p < 0.05)'),
    ]
    fig.legend(
        handles=legend_elements,
        loc='upper left',
        bbox_to_anchor=(0.02, 0.98),
        fontsize=8,
        frameon=True,
        borderpad=0.4,
        labelspacing=0.3,
        handlelength=1.4,
        handletextpad=0.6,
    )
    fig.tight_layout(rect=[0, 0.02, 1, 0.90])
    fig.savefig(outpath, bbox_inches='tight')
    plt.close(fig)


# Proliferation proxy (Ki67-like) per patient
prolif_types = {'Ki67+ proliferative tumor cells', 'HER2+ Ki67+ proliferative tumor cells'}
prolif_df = pd.DataFrame({
    'Patient': adata.obs['Patient'].astype(str).values,
    'is_prolif': adata.obs['celltype'].astype(str).isin(prolif_types).values,
})
prolif_df = prolif_df.groupby('Patient', as_index=False).agg(
    n_cells=('is_prolif', 'size'),
    n_prolif=('is_prolif', 'sum')
)
prolif_df['prolif_frac'] = prolif_df['n_prolif'] / prolif_df['n_cells']

cox_df = df_merged.copy()
cox_df['Patient'] = cox_df['Patient'].astype(str)
cox_df = cox_df.merge(prolif_df[['Patient', 'prolif_frac']], on='Patient', how='left')

# Encode cohort if variable
if 'cohort' in cox_df.columns:
    cox_df['cohort'] = cox_df['cohort'].astype('category')
    cohort_dummies = pd.get_dummies(cox_df['cohort'], prefix='cohort', drop_first=True)
    cox_df = pd.concat([cox_df, cohort_dummies], axis=1)

# Standard covariates
cox_df['age_z'] = _zscore(cox_df['age']) if 'age' in cox_df.columns else np.nan
cox_df['grade_z'] = _zscore(cox_df['grade']) if 'grade' in cox_df.columns else np.nan
cox_df['stage_z'] = _zscore(cox_df['Stage']) if 'Stage' in cox_df.columns else np.nan
cox_df['n_stage_z'] = _zscore(cox_df['N_Stage']) if 'N_Stage' in cox_df.columns else np.nan
cox_df['prolif_frac_z'] = _zscore(cox_df['prolif_frac'])

# Mixing terms: per +0.1 increase
cox_df['overall_mixing_per_0p1'] = pd.to_numeric(cox_df['overall_mixing'], errors='coerce') / 0.1
for term in ['overall_cd3_mix', 'overall_cd68_mix', 'overall_cd20_mix']:
    if term in cox_df.columns:
        cox_df[f'{term}_per_0p1'] = pd.to_numeric(cox_df[term], errors='coerce') / 0.1

base_covars = ['age_z', 'grade_z', 'stage_z', 'n_stage_z', 'prolif_frac_z']
base_covars += [c for c in cox_df.columns if c.startswith('cohort_')]

# Fit overall mixing Cox
overall_rows = []
overall_term = 'overall_mixing_per_0p1'
covars_overall = [overall_term, *base_covars]

os_summ = pd.DataFrame()
rfs_summ = pd.DataFrame()

if cox_df['OS_months'].notna().any():
    os_summ = _fit_cox_table(
        cox_df.dropna(subset=['OS_months', 'OS_event']),
        'OS_months',
        'OS_event',
        covars_overall,
        penalizer=0.1,
        label='OS'
    )
    if not os_summ.empty:
        os_summ.to_csv(OUTPUT_DIR / 'ERpos_step10_cox_OS_full_overall_mixing_per0p1.csv', index=False)
    row = os_summ[os_summ['covariate'] == overall_term].copy()
    if not row.empty:
        row['outcome'] = 'OS'
        overall_rows.append(row)

if cox_df['RFS_months'].notna().any():
    rfs_summ = _fit_cox_table(
        cox_df.dropna(subset=['RFS_months', 'RFS_event']),
        'RFS_months',
        'RFS_event',
        covars_overall,
        penalizer=0.1,
        label='RFS'
    )
    if not rfs_summ.empty:
        rfs_summ.to_csv(OUTPUT_DIR / 'ERpos_step10_cox_RFS_full_overall_mixing_per0p1.csv', index=False)
    row = rfs_summ[rfs_summ['covariate'] == overall_term].copy()
    if not row.empty:
        row['outcome'] = 'RFS'
        overall_rows.append(row)

overall_out = pd.concat(overall_rows, ignore_index=True) if overall_rows else pd.DataFrame()
if not overall_out.empty:
    overall_out.to_csv(OUTPUT_DIR / 'ERpos_step10_cox_overall_mixing_per0p1.csv', index=False)
    _forest_plot_multivariable_full(
        os_summ,
        rfs_summ,
        [*base_covars, overall_term],
        FIGURES_DIR / '04_cox_forest_overall_mixing.png',
    )
    print("  Saved: 04_cox_forest_overall_mixing.png")

# Fit cell-type-specific mixing Cox (separate models per term)
celltype_terms = [
    ('mix_cd3', 'overall_cd3_mix_per_0p1'),
    ('mix_cd68', 'overall_cd68_mix_per_0p1'),
    ('mix_cd20', 'overall_cd20_mix_per_0p1'),
]

cell_rows = []
cell_full_rows = []
for outcome, time_col, event_col in [('OS', 'OS_months', 'OS_event'), ('RFS', 'RFS_months', 'RFS_event')]:
    if time_col not in cox_df.columns or event_col not in cox_df.columns:
        continue
    dfx = cox_df.dropna(subset=[time_col, event_col]).copy()
    for term_name, term_col in celltype_terms:
        if term_col not in dfx.columns:
            continue
        covs = [term_col, *base_covars]
        summ = _fit_cox_table(dfx, time_col, event_col, covs, penalizer=0.1, label=f'{outcome}:{term_name}')
        if not summ.empty:
            summ2 = summ.copy()
            summ2['outcome'] = outcome
            summ2['term'] = term_name
            cell_full_rows.append(summ2)
        row = summ[summ['covariate'] == term_col].copy()
        if row.empty:
            continue
        row['outcome'] = outcome
        row['term'] = term_name
        cell_rows.append(row[['outcome', 'term', 'coef', 'HR', 'HR_lower_95', 'HR_upper_95', 'p', 'n']])

cell_out = pd.concat(cell_rows, ignore_index=True) if cell_rows else pd.DataFrame()
if not cell_out.empty:
    cell_out.to_csv(OUTPUT_DIR / 'ERpos_step10_cox_celltype_mixing_per0p1.csv', index=False)
    _forest_plot_celltypes(cell_out, FIGURES_DIR / '05_cox_forest_celltype_mixing.png')
    print("  Saved: 05_cox_forest_celltype_mixing.png")

cell_full_out = pd.concat(cell_full_rows, ignore_index=True) if cell_full_rows else pd.DataFrame()
if not cell_full_out.empty:
    cell_full_out.to_csv(OUTPUT_DIR / 'ERpos_step10_cox_celltype_mixing_full_per0p1.csv', index=False)

# Proliferation stratification + interaction (median split)
print("Generating Plot 2c: Proliferation stratification + interaction...")

if 'prolif_frac' in cox_df.columns and cox_df['prolif_frac'].notna().any():
    prolif_median = float(cox_df['prolif_frac'].median())
    cox_df['prolif_high'] = (cox_df['prolif_frac'] >= prolif_median).astype(int)
    cox_df['mix_x_prolif'] = cox_df['overall_mixing_per_0p1'] * cox_df['prolif_high']

    strat_rows = []
    for outcome, time_col, event_col in [('OS', 'OS_months', 'OS_event'), ('RFS', 'RFS_months', 'RFS_event')]:
        if time_col not in cox_df.columns or event_col not in cox_df.columns:
            continue
        dfx = cox_df.dropna(subset=[time_col, event_col]).copy()

        for group_name, group_mask in [('prolif_low', dfx['prolif_high'] == 0), ('prolif_high', dfx['prolif_high'] == 1)]:
            df_g = dfx[group_mask].copy()
            if df_g.empty:
                continue
            summ = _fit_cox_table(
                df_g,
                time_col,
                event_col,
                [overall_term, *base_covars],
                penalizer=0.1,
                label=f'{outcome}:{group_name}',
            )
            if summ.empty:
                continue
            row = summ[summ['covariate'] == overall_term].copy()
            if row.empty:
                continue
            row['outcome'] = outcome
            row['group'] = group_name
            row['term'] = 'overall_mixing_per_0p1'
            strat_rows.append(row[['outcome', 'group', 'term', 'coef', 'HR', 'HR_lower_95', 'HR_upper_95', 'p', 'n']])

        # Interaction model (omit prolif_frac_z to reduce redundancy with prolif_high)
        interaction_covars = ['overall_mixing_per_0p1', 'prolif_high', 'mix_x_prolif', 'age_z', 'grade_z', 'stage_z', 'n_stage_z']
        interaction_covars += [c for c in dfx.columns if c.startswith('cohort_')]
        summ_int = _fit_cox_table(
            dfx,
            time_col,
            event_col,
            interaction_covars,
            penalizer=0.1,
            label=f'{outcome}:interaction',
        )
        if not summ_int.empty and 'mix_x_prolif' in summ_int['covariate'].values:
            row = summ_int[summ_int['covariate'] == 'mix_x_prolif'].copy()
            row['outcome'] = outcome
            row['group'] = 'interaction'
            row['term'] = 'mix_x_prolif'
            strat_rows.append(row[['outcome', 'group', 'term', 'coef', 'HR', 'HR_lower_95', 'HR_upper_95', 'p', 'n']])

    strat_out = pd.concat(strat_rows, ignore_index=True) if strat_rows else pd.DataFrame()
    if not strat_out.empty:
        strat_out.to_csv(OUTPUT_DIR / 'ERpos_step10_cox_overall_mixing_stratified_by_prolif.csv', index=False)
        print("  Saved: ERpos_step10_cox_overall_mixing_stratified_by_prolif.csv")
else:
    print("  Skipped: proliferation proxy unavailable")

# ========================================================================
# PLOT 3: Spatial visualization of mixing (example ROIs)
# ========================================================================
print("Generating Plot 3: Spatial visualization of mixing...")

def _select_example_scenes(df_in, n_pairs=2):
    """
    Select example ROIs for visualization.

    Prefer ROIs with both tumor and immune cells present so that "mixing" is
    visually interpretable (avoid trivial cases with ~0 immune cells).

    Additionally, match low-mixing examples to high-mixing examples by
    (approximate) tumor/immune cell counts to make the comparison fair.
    """
    threshold_schedule = [
        (200, 20),
        (100, 10),
        (50, 5),
        (25, 2),
        (0, 1),
        (0, 0),
    ]

    mix_col = 'overall_mixing'

    for min_tumor, min_immune in threshold_schedule:
        candidates = df_in[
            (df_in['n_tumor'] >= min_tumor)
            & (df_in['n_immune'] >= min_immune)
            & df_in[mix_col].notna()
        ].copy()

        if len(candidates) < (2 * n_pairs):
            continue

        high_cut = candidates[mix_col].quantile(0.80)
        low_cut = candidates[mix_col].quantile(0.20)
        high_pool = candidates[candidates[mix_col] >= high_cut].sort_values(mix_col, ascending=False)
        low_pool = candidates[candidates[mix_col] <= low_cut].sort_values(mix_col, ascending=True)

        if len(high_pool) < n_pairs or len(low_pool) < n_pairs:
            continue

        # High mixing examples:
        # - Left panel (primary): enforce minimum counts for interpretability.
        # - Right panel (partner): keep a strong high-mixing example (highest mixing remaining).
        high_candidates = high_pool.head(max(50, 10 * n_pairs)).copy()

        primary_scene = None
        # Requested constraint: at least 2000 tumor and 100 immune (if available).
        for min_immune_high in [100, 75, 50, 40, 30, 20, 10, min_immune]:
            pool = high_candidates[
                (high_candidates['n_tumor'] >= 2000)
                & (high_candidates['n_immune'] >= max(min_immune, min_immune_high))
            ]
            if pool.empty:
                continue
            primary_scene = pool.sort_values(mix_col, ascending=False).iloc[0]['slide_scene']
            break

        if primary_scene is None:
            # Fallback: keep pipeline running if dataset has no such ROI.
            # Prefer tumor-rich + immune-rich within the high-mixing pool.
            primary_scene = high_candidates.sort_values(
                ['n_tumor', 'n_immune', 'overall_mixing'],
                ascending=False
            ).iloc[0]['slide_scene']

        remaining_high = high_candidates[high_candidates['slide_scene'] != primary_scene].copy()
        if remaining_high.empty:
            continue

        partner_scene = None
        for min_immune_partner in [100, 75, 50, 40, 30, 20, 10, min_immune]:
            partner_pool = remaining_high[remaining_high['n_immune'] >= max(min_immune, min_immune_partner)]
            if partner_pool.empty:
                continue
            partner_scene = partner_pool.sort_values(mix_col, ascending=False).iloc[0]['slide_scene']
            break

        if partner_scene is None:
            # Fallback: keep pipeline running if immune-rich high-mixing example is not available.
            partner_scene = remaining_high.sort_values(mix_col, ascending=False).iloc[0]['slide_scene']
            print(
                "  Warning: could not find a second high-mixing ROI with the desired minimum immune count; "
                "using the best available high-mixing ROI."
            )

        # Keep ordering stable: primary first (left), partner second (right)
        highs = candidates[candidates['slide_scene'].isin([primary_scene, partner_scene])].copy()
        highs = highs.set_index('slide_scene').loc[[primary_scene, partner_scene]].reset_index()
        lows = []
        used = set(highs['slide_scene'].tolist())

        for _, high_row in highs.iterrows():
            pool = low_pool[~low_pool['slide_scene'].isin(used)]
            if pool.empty:
                break

            target = np.array([
                np.log1p(high_row['n_tumor']),
                np.log1p(high_row['n_immune']),
                np.log1p(high_row['n_cells']),
            ])
            pool_mat = np.vstack([
                np.log1p(pool['n_tumor'].to_numpy()),
                np.log1p(pool['n_immune'].to_numpy()),
                np.log1p(pool['n_cells'].to_numpy()),
            ]).T

            dist = np.linalg.norm(pool_mat - target, axis=1)
            best_idx = int(np.argmin(dist))
            best_scene = pool.iloc[best_idx]['slide_scene']
            lows.append(best_scene)
            used.add(best_scene)

        if len(lows) == n_pairs:
            return highs['slide_scene'].tolist(), lows

    # Fallback: extremes from all scenes (keeps pipeline running even if immune is rare)
    high = df_in[df_in['overall_mixing'].notna()].nlargest(n_pairs, 'overall_mixing')['slide_scene'].tolist()
    remaining = df_in[~df_in['slide_scene'].isin(high)]
    low = remaining[remaining['overall_mixing'].notna()].nsmallest(n_pairs, 'overall_mixing')['slide_scene'].tolist()
    return high, low


# Select 4 example ROIs: 2 with high mixing, 2 with low mixing (matched counts)
high_mixing_scenes, low_mixing_scenes = _select_example_scenes(df, n_pairs=2)
print(f"  Selected high-mixing scenes: {high_mixing_scenes}")
print(f"  Selected low-mixing scenes:  {low_mixing_scenes}")

selected_scenes = list(high_mixing_scenes) + list(low_mixing_scenes)

# Standardize axis ranges across panels (same X/Y dimensions in μm)
scene_windows = {}
global_range = 0.0
for scene in selected_scenes:
    scene_mask = adata.obs['slide_scene'] == scene
    scene_coords = coords_um[scene_mask, :]
    x_min, x_max = float(scene_coords[:, 0].min()), float(scene_coords[:, 0].max())
    y_min, y_max = float(scene_coords[:, 1].min()), float(scene_coords[:, 1].max())
    x_c = (x_min + x_max) / 2.0
    y_c = (y_min + y_max) / 2.0
    side = max(x_max - x_min, y_max - y_min)
    scene_windows[scene] = (x_c, y_c, side)
    global_range = max(global_range, side)

global_half_range = (global_range * 1.05) / 2.0  # small padding

fig, axes = plt.subplots(2, 2, figsize=(14, 12), dpi=300)

# Plot high mixing examples
for idx, scene in enumerate(high_mixing_scenes):
    ax = axes[0, idx]

    scene_mask = adata.obs['slide_scene'] == scene
    scene_data = adata.obs.loc[scene_mask]
    scene_coords = coords_um[scene_mask, :]

    compartment = scene_data['compartment'].values
    colors = {'tumor': 'red', 'immune': 'blue', 'other': 'green'}
    for comp, size, alpha in [('other', 4, 0.35), ('tumor', 5, 0.65), ('immune', 12, 0.9)]:
        comp_mask = compartment == comp
        if comp_mask.any():
            ax.scatter(
                scene_coords[comp_mask, 0],
                scene_coords[comp_mask, 1],
                c=colors[comp],
                s=size,
                alpha=alpha,
                linewidths=0,
            )

    mixing_score = df[df['slide_scene'] == scene]['overall_mixing'].iloc[0]
    n_tumor = df[df['slide_scene'] == scene]['n_tumor'].iloc[0]
    n_immune = df[df['slide_scene'] == scene]['n_immune'].iloc[0]

    ax.set_title(f'High Mixing (score={mixing_score:.3f})\nTumor: {n_tumor}, Immune: {n_immune}',
                 fontsize=10, fontweight='bold')
    ax.set_xlabel('X (μm)', fontsize=9)
    ax.set_ylabel('Y (μm)', fontsize=9)
    x_c, y_c, _side = scene_windows[scene]
    ax.set_xlim(x_c - global_half_range, x_c + global_half_range)
    ax.set_ylim(y_c - global_half_range, y_c + global_half_range)
    ax.set_aspect('equal', adjustable='box')

# Plot low mixing examples
for idx, scene in enumerate(low_mixing_scenes):
    ax = axes[1, idx]

    scene_mask = adata.obs['slide_scene'] == scene
    scene_data = adata.obs.loc[scene_mask]
    scene_coords = coords_um[scene_mask, :]

    compartment = scene_data['compartment'].values
    colors = {'tumor': 'red', 'immune': 'blue', 'other': 'green'}
    for comp, size, alpha in [('other', 4, 0.35), ('tumor', 5, 0.65), ('immune', 12, 0.9)]:
        comp_mask = compartment == comp
        if comp_mask.any():
            ax.scatter(
                scene_coords[comp_mask, 0],
                scene_coords[comp_mask, 1],
                c=colors[comp],
                s=size,
                alpha=alpha,
                linewidths=0,
            )

    mixing_score = df[df['slide_scene'] == scene]['overall_mixing'].iloc[0]
    n_tumor = df[df['slide_scene'] == scene]['n_tumor'].iloc[0]
    n_immune = df[df['slide_scene'] == scene]['n_immune'].iloc[0]

    ax.set_title(f'Low Mixing (score={mixing_score:.3f})\nTumor: {n_tumor}, Immune: {n_immune}',
                 fontsize=10, fontweight='bold')
    ax.set_xlabel('X (μm)', fontsize=9)
    ax.set_ylabel('Y (μm)', fontsize=9)
    x_c, y_c, _side = scene_windows[scene]
    ax.set_xlim(x_c - global_half_range, x_c + global_half_range)
    ax.set_ylim(y_c - global_half_range, y_c + global_half_range)
    ax.set_aspect('equal', adjustable='box')

# Add legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='red', alpha=0.6, label='Tumor (epithelial)'),
    Patch(facecolor='blue', alpha=0.6, label='Immune'),
    Patch(facecolor='green', alpha=0.6, label='Other (stromal)')
]
fig.legend(
    handles=legend_elements,
    loc='lower center',
    bbox_to_anchor=(0.5, 0.01),
    ncol=3,
    fontsize=10,
    frameon=True,
)

plt.suptitle('Spatial Visualization of Tumor-Immune Mixing\n(Examples of High vs Low Mixing)',
             fontsize=13, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0.08, 1, 0.96])
plt.savefig(FIGURES_DIR / '03_mixing_spatial_examples.png', dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: 03_mixing_spatial_examples.png")

print(f"\n✅ Step 10 visualizations complete! (figures/step10_mixing_score/)")

print("\n" + "=" * 80)
print("✓ STEP 10 COMPLETE")
print("=" * 80)
print()
