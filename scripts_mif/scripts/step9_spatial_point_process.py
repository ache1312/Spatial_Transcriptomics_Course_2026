#!/usr/bin/env python
"""
STEP 9: Spatial Point Process Metrics

ISOLATED PIPELINE - Run from pipeline_leiden_validated/ directory:
  cd /path/to/pipeline_leiden_validated
  conda run -p /path/to/squidpy-env python scripts/step9_spatial_point_process.py

Implements spatial point process metrics mentioned in Keren et al. Supplemental Figure 11:
1. Ripley's L function (tests for clustering/dispersion)
2. K-cross function (cross-type spatial correlation)
3. G-cross function (nearest-neighbor cross-type distribution)

As noted in the paper:
"Additional spatial metrics, including Ripley's L, Kcross, and Gcross, did not yield
any significant biomarkers in the validation cohort (Supplemental Figure 11A)."

However, these metrics demonstrate methodological rigor and completeness.

Input:  data/imc_ERpos_with_embeddings_and_full_names.h5ad (with full cell type names)
Output: analysis_imc/ERpos_spatial_point_process_by_patient.csv
        figures/step9_spatial_point_process/*.png (4 visualization plots)
"""

import os
import time
import warnings
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.spatial import distance_matrix
from pathlib import Path

warnings.filterwarnings('ignore')

# Settings
RADII = np.linspace(0, 100, 21)  # 0, 5, 10, ..., 100 microns
MICRONS_PER_PIXEL = 1.0

OUTPUT_DIR = Path("analysis_imc")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

print("=" * 80)
print("STEP 9: SPATIAL POINT PROCESS METRICS (Ripley's L, K-cross, G-cross)")
print("=" * 80)
print(f"Radii: 0 - 100 μm ({len(RADII)} steps)")
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

# Get spatial coordinates in microns
coords_pixels = adata.obsm['spatial']
coords_um = coords_pixels * MICRONS_PER_PIXEL


def ripleys_L(coords, radii, area):
    """
    Calculate Ripley's L function for a point pattern

    L(r) = sqrt(K(r) / pi) - r

    where K(r) = (Area / n^2) * sum(I(d_ij < r))

    L(r) ≈ 0: random (Poisson)
    L(r) > 0: clustering
    L(r) < 0: dispersion/regularity

    Args:
        coords: Nx2 array of coordinates (microns)
        radii: Array of radii to compute (microns)
        area: Total area (square microns)

    Returns:
        Array of L(r) values
    """
    n = len(coords)

    if n < 2:
        return np.zeros(len(radii))

    # Compute pairwise distances
    dists = distance_matrix(coords, coords)

    # For each radius, count pairs within radius
    L_values = []

    for r in radii:
        if r == 0:
            L_values.append(0)
            continue

        # Count pairs within radius (excluding self, i≠j)
        count = np.sum((dists < r) & (dists > 0)) / 2  # Divide by 2 to avoid double counting

        # K(r) estimation (edge correction ignored for simplicity)
        K_r = (area / (n * (n - 1))) * 2 * count

        # L(r) transformation
        L_r = np.sqrt(K_r / np.pi) - r

        L_values.append(L_r)

    return np.array(L_values)


def k_cross(coords_i, coords_j, radii, area):
    """
    Calculate K-cross function between two point types

    K_ij(r) measures spatial correlation between type i and type j

    K_ij(r) > K_Poisson(r): attraction between types
    K_ij(r) < K_Poisson(r): repulsion between types

    Args:
        coords_i: Nx2 array for type i
        coords_j: Mx2 array for type j
        radii: Array of radii
        area: Total area

    Returns:
        Array of K_cross(r) values
    """
    n_i = len(coords_i)
    n_j = len(coords_j)

    if n_i == 0 or n_j == 0:
        return np.zeros(len(radii))

    # Cross-type distances
    dists = distance_matrix(coords_i, coords_j)

    K_values = []

    for r in radii:
        if r == 0:
            K_values.append(0)
            continue

        # Count type j points within r of each type i point
        count = np.sum(dists < r)

        # K_cross estimation
        K_ij = (area / (n_i * n_j)) * count

        K_values.append(K_ij)

    return np.array(K_values)


def g_cross(coords_i, coords_j):
    """
    Calculate G-cross function: distribution of nearest-neighbor distances
    from type i to type j

    G_ij(r) = proportion of type i points with nearest type j neighbor at distance < r

    Args:
        coords_i: Nx2 array for type i
        coords_j: Mx2 array for type j

    Returns:
        mean nearest-neighbor distance from i to j
    """
    n_i = len(coords_i)
    n_j = len(coords_j)

    if n_i == 0 or n_j == 0:
        return np.nan

    # Cross-type distances
    dists = distance_matrix(coords_i, coords_j)

    # For each type i point, find nearest type j point
    nn_dists = np.min(dists, axis=1)

    # Return mean nearest-neighbor distance
    return np.mean(nn_dists)


print("\n[2/4] Calculating spatial point process metrics per ROI...")
t0 = time.time()

scenes = adata.obs['slide_scene'].unique()
total = len(scenes)

all_metrics = []

for i, scene in enumerate(scenes, 1):
    if i % 50 == 0 or i == total:
        print(f"  Processing ROI {i}/{total}...", flush=True)

    scene_mask = adata.obs['slide_scene'] == scene
    scene_data = adata.obs.loc[scene_mask]
    scene_coords = coords_um[scene_mask, :]

    if len(scene_coords) < 10:
        continue

    patient = scene_data['Patient'].iloc[0]

    # Calculate bounding box area
    x_min, y_min = scene_coords.min(axis=0)
    x_max, y_max = scene_coords.max(axis=0)
    area = (x_max - x_min) * (y_max - y_min)

    # Ripley's L for each cell type
    for celltype in scene_data['celltype'].unique():
        ct_mask = scene_data['celltype'] == celltype
        ct_coords = scene_coords[ct_mask.values, :]

        if len(ct_coords) >= 2:
            L_values = ripleys_L(ct_coords, RADII, area)

            # Store key summary statistics
            all_metrics.append({
                'slide_scene': scene,
                'Patient': patient,
                'metric': 'ripleys_L',
                'celltype': celltype,
                'celltype_pair': celltype,  # For consistency with K-cross
                'L_mean': np.mean(L_values),
                'L_max': np.max(L_values),
                'L_at_40um': L_values[8] if len(L_values) > 8 else np.nan,  # r=40 is index 8
                'n_cells': len(ct_coords)
            })

    # K-cross and G-cross for pairs of cell types
    celltypes = scene_data['celltype'].unique()

    for i, ct_i in enumerate(celltypes):
        ct_i_mask = scene_data['celltype'] == ct_i
        coords_i = scene_coords[ct_i_mask.values, :]

        for ct_j in celltypes[i:]:  # Only upper triangle (avoid redundancy)
            if ct_i == ct_j:
                continue

            ct_j_mask = scene_data['celltype'] == ct_j
            coords_j = scene_coords[ct_j_mask.values, :]

            if len(coords_i) > 0 and len(coords_j) > 0:
                # K-cross
                K_values = k_cross(coords_i, coords_j, RADII, area)

                all_metrics.append({
                    'slide_scene': scene,
                    'Patient': patient,
                    'metric': 'K_cross',
                    'celltype': ct_i,
                    'celltype_pair': f'{ct_i}-{ct_j}',
                    'K_mean': np.mean(K_values),
                    'K_at_40um': K_values[8] if len(K_values) > 8 else np.nan,
                    'n_cells': len(coords_i)
                })

                # G-cross (mean nearest-neighbor distance)
                g_ij = g_cross(coords_i, coords_j)
                g_ji = g_cross(coords_j, coords_i)

                all_metrics.append({
                    'slide_scene': scene,
                    'Patient': patient,
                    'metric': 'G_cross',
                    'celltype': ct_i,
                    'celltype_pair': f'{ct_i}-{ct_j}',
                    'G_mean_nn_dist': g_ij,
                    'G_mean_nn_dist_reverse': g_ji,
                    'n_cells': len(coords_i)
                })

print(f"Time: {time.time() - t0:.1f}s")

print("\n[3/4] Aggregating to patient level...")
df = pd.DataFrame(all_metrics)

# Aggregate to patient level
patient_metrics = []

for patient in df['Patient'].unique():
    patient_df = df[df['Patient'] == patient]

    # Ripley's L summary
    for celltype in patient_df['celltype'].unique():
        L_df = patient_df[
            (patient_df['metric'] == 'ripleys_L') &
            (patient_df['celltype'] == celltype)
        ]

        if len(L_df) > 0:
            patient_metrics.append({
                'Patient': patient,
                'metric_type': 'ripleys_L',
                'celltype': celltype,
                'celltype_pair': celltype,  # FIX: Add celltype_pair for pivot compatibility
                'value': L_df['L_mean'].mean(),
                'value_at_40um': L_df['L_at_40um'].mean()
            })

    # K-cross summary
    for pair in patient_df[patient_df['metric'] == 'K_cross']['celltype_pair'].unique():
        K_df = patient_df[
            (patient_df['metric'] == 'K_cross') &
            (patient_df['celltype_pair'] == pair)
        ]

        if len(K_df) > 0:
            patient_metrics.append({
                'Patient': patient,
                'metric_type': 'K_cross',
                'celltype': pair.split('-')[0],
                'celltype_pair': pair,
                'value': K_df['K_mean'].mean(),
                'value_at_40um': K_df['K_at_40um'].mean()
            })

    # G-cross summary
    for pair in patient_df[patient_df['metric'] == 'G_cross']['celltype_pair'].unique():
        G_df = patient_df[
            (patient_df['metric'] == 'G_cross') &
            (patient_df['celltype_pair'] == pair)
        ]

        if len(G_df) > 0:
            patient_metrics.append({
                'Patient': patient,
                'metric_type': 'G_cross',
                'celltype': pair.split('-')[0],
                'celltype_pair': pair,
                'value': G_df['G_mean_nn_dist'].mean()
            })

df_patient = pd.DataFrame(patient_metrics)

# Pivot to wide format
df_wide = df_patient.pivot_table(
    index='Patient',
    columns=['metric_type', 'celltype_pair'],
    values='value',
    fill_value=0
)

# Flatten column names
df_wide.columns = [f'{metric}_{pair}' for metric, pair in df_wide.columns]
df_wide = df_wide.reset_index()

print("\n[4/4] Saving results...")
output_file = OUTPUT_DIR / "ERpos_spatial_point_process_by_patient.csv"
df_wide.to_csv(output_file, index=False)
print(f"Saved: {output_file}")
print(f"Shape: {df_wide.shape[0]} patients × {df_wide.shape[1]-1} metrics")

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

print("\nMetrics calculated:")
print(f"  - Ripley's L (per celltype): Clustering indicator")
print(f"  - K-cross (celltype pairs): Spatial correlation")
print(f"  - G-cross (celltype pairs): Nearest-neighbor distances")

print("\nInterpretation:")
print("  - Ripley's L > 0: Clustering (cells aggregate)")
print("  - Ripley's L < 0: Dispersion (cells avoid each other)")
print("  - K-cross high: Types co-locate")
print("  - K-cross low: Types segregate")
print("  - G-cross low: Types are close neighbors")

print("\nNote from Keren et al.:")
print("  'Additional spatial metrics, including Ripley's L, Kcross, and Gcross,")
print("   did not yield any significant biomarkers in the validation cohort'")
print("  However, these metrics provide comprehensive spatial characterization.")

# ============================================================================
# VISUALIZATIONS
# ============================================================================
print("\n[GENERATING VISUALIZATIONS]")

import matplotlib.pyplot as plt
import seaborn as sns

# Set matplotlib parameters
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

# Create figures directory
FIGURES_DIR = Path("figures/step9_spatial_point_process")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ========================================================================
# PLOT 1: Ripley's L function curves for example cell types
# ========================================================================
print("Generating Plot 1: Ripley's L function curves...")

# Select 4 most abundant cell types for demonstration
top_celltypes = adata.obs['celltype'].value_counts().head(4).index.tolist()

# For each cell type, calculate average L(r) curve across patients
fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=300)
axes = axes.flatten()

for idx, celltype in enumerate(top_celltypes):
    ax = axes[idx]

    # Get all ROIs for this celltype
    ct_L_data = df[
        (df['metric'] == 'ripleys_L') &
        (df['celltype'] == celltype)
    ]

    if len(ct_L_data) > 0:
        # For simplicity, calculate L for a few example ROIs
        example_scenes = ct_L_data['slide_scene'].unique()[:5]

        for scene in example_scenes:
            scene_mask = adata.obs['slide_scene'] == scene
            scene_coords = coords_um[scene_mask, :]
            scene_data = adata.obs.loc[scene_mask]

            ct_mask = scene_data['celltype'] == celltype
            ct_coords = scene_coords[ct_mask.values, :]

            if len(ct_coords) >= 2:
                # Calculate area
                x_min, y_min = scene_coords.min(axis=0)
                x_max, y_max = scene_coords.max(axis=0)
                area = (x_max - x_min) * (y_max - y_min)

                # Calculate L(r)
                L_values = ripleys_L(ct_coords, RADII, area)

                ax.plot(RADII, L_values, alpha=0.3, color='steelblue', linewidth=0.8)

        # Add CSR reference line (L=0)
        ax.axhline(0, color='red', linestyle='--', linewidth=1.5, label='CSR (L=0)')

        ax.set_xlabel('Radius (μm)', fontsize=10)
        ax.set_ylabel("Ripley's L(r)", fontsize=10)
        ax.set_title(f"{celltype}\n(n={ct_L_data['n_cells'].sum():,} cells)", fontsize=10, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)

plt.suptitle("Ripley's L Function: Spatial Clustering Analysis\n(L>0: clustering, L<0: dispersion, L≈0: random)",
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig(FIGURES_DIR / '01_ripleys_L_curves.png', dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: 01_ripleys_L_curves.png")

# ========================================================================
# PLOT 2: G-cross function (nearest-neighbor distances for key pairs)
# ========================================================================
print("Generating Plot 2: G-cross nearest-neighbor distances...")

# Select top 10 pairs with most data
g_cross_data = df_patient[df_patient['metric_type'] == 'G_cross'].copy()

if len(g_cross_data) > 0:
    # Count how many patients have each pair
    pair_counts = g_cross_data.groupby('celltype_pair').size().sort_values(ascending=False)
    top_pairs = pair_counts.head(10).index.tolist()

    # Get mean values for these pairs
    pair_values = []
    for pair in top_pairs:
        pair_df = g_cross_data[g_cross_data['celltype_pair'] == pair]
        mean_val = pair_df['value'].mean()
        pair_values.append({'pair': pair, 'mean_nn_dist': mean_val})

    df_gcross_plot = pd.DataFrame(pair_values).sort_values('mean_nn_dist')

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

    bars = ax.barh(range(len(df_gcross_plot)), df_gcross_plot['mean_nn_dist'].values,
                   color='coral', alpha=0.8)

    ax.set_yticks(range(len(df_gcross_plot)))
    ax.set_yticklabels(df_gcross_plot['pair'].values, fontsize=9)
    ax.set_xlabel('Mean Nearest-Neighbor Distance (μm)', fontsize=11)
    ax.set_title('G-cross Function: Cross-Type Nearest-Neighbor Distances\nTop 10 Cell Type Pairs',
                 fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')

    # Add value labels
    for i, value in enumerate(df_gcross_plot['mean_nn_dist'].values):
        ax.text(value + 1, i, f'{value:.1f}μm', va='center', fontsize=8)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / '02_gcross_nn_distances.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: 02_gcross_nn_distances.png")
else:
    print("  Skipped: No G-cross data available")

# ========================================================================
# PLOT 3: Spatial clustering scores by patient (Ripley's L at 40μm)
# ========================================================================
print("Generating Plot 3: Spatial clustering scores by patient...")

# Extract L values at 40μm for top cell types
L_data_40um = df_patient[
    (df_patient['metric_type'] == 'ripleys_L') &
    (df_patient['celltype'].isin(top_celltypes))
].copy()

if len(L_data_40um) > 0 and 'value_at_40um' in L_data_40um.columns:
    # Pivot to patient × celltype matrix
    L_matrix = L_data_40um.pivot_table(
        index='Patient',
        columns='celltype',
        values='value_at_40um',
        fill_value=0
    )

    # Select patients with most data (top 30)
    patient_counts = L_matrix.count(axis=1).sort_values(ascending=False)
    top_patients = patient_counts.head(30).index.tolist()
    L_matrix_top = L_matrix.loc[top_patients, top_celltypes]

    fig, ax = plt.subplots(figsize=(10, 12), dpi=300)

    sns.heatmap(
        L_matrix_top,
        cmap='RdBu_r',
        center=0,
        robust=True,
        cbar_kws={'label': "Ripley's L(40μm)", 'shrink': 0.5},
        yticklabels=False,  # Too many patients
        xticklabels=True,
        ax=ax
    )

    ax.set_xlabel('Cell Type', fontsize=11)
    ax.set_ylabel(f'Patients (n={len(top_patients)})', fontsize=11)
    ax.set_title("Spatial Clustering Scores by Patient\n(Ripley's L at r=40μm)",
                 fontsize=12, fontweight='bold')
    plt.xticks(rotation=45, ha='right', fontsize=9)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / '03_clustering_by_patient.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: 03_clustering_by_patient.png")
else:
    print("  Skipped: No L(40μm) data available")

# ========================================================================
# PLOT 4: CSR comparison (distribution of L values)
# ========================================================================
print("Generating Plot 4: CSR comparison...")

# Compare observed L values against CSR expectation (L=0)
L_all_data = df_patient[df_patient['metric_type'] == 'ripleys_L'].copy()

if len(L_all_data) > 0:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=300)

    # Plot 1: Distribution of L values
    ax = axes[0]
    sns.histplot(
        data=L_all_data,
        x='value',
        bins=50,
        kde=True,
        color='steelblue',
        ax=ax
    )

    ax.axvline(0, color='red', linestyle='--', linewidth=2, label='CSR (L=0)')
    ax.set_xlabel("Ripley's L", fontsize=11)
    ax.set_ylabel('Frequency', fontsize=11)
    ax.set_title('Distribution of Spatial Clustering Scores\n(All cell types, all patients)',
                 fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Plot 2: L values by cell type (boxplot)
    ax = axes[1]

    L_top_ct = L_all_data[L_all_data['celltype'].isin(top_celltypes)]

    sns.boxplot(
        data=L_top_ct,
        x='celltype',
        y='value',
        palette='Set2',
        ax=ax
    )

    ax.axhline(0, color='red', linestyle='--', linewidth=1.5, label='CSR (L=0)')
    ax.set_xlabel('Cell Type', fontsize=11)
    ax.set_ylabel("Ripley's L", fontsize=11)
    ax.set_title('Spatial Clustering by Cell Type\n(vs Complete Spatial Randomness)',
                 fontsize=12, fontweight='bold')
    ax.tick_params(axis='x', rotation=45)
    ax.grid(True, alpha=0.3, axis='y')
    ax.legend(fontsize=9)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right', fontsize=9)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / '04_csr_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: 04_csr_comparison.png")
else:
    print("  Skipped: No Ripley's L data for CSR comparison")

print(f"\n✅ Step 9 visualizations complete! (figures/step9_spatial_point_process/)")

print("\n" + "=" * 80)
print("✓ STEP 9 COMPLETE")
print("=" * 80)
print()
