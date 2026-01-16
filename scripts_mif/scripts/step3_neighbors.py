#!/usr/bin/env python3
"""
STEP 3: Spatial Neighbor Graphs & Composition Analysis (Squidpy)

Purpose:
    Compute neighborhood compositions at multiple spatial scales using squidpy
    for both COARSE (5 types) and FINE (20 types) granularities.

Key Changes from Previous Version:
    - Processes per ROI (slide_scene) to avoid cross-ROI edges
    - Calculates for both granularities in a single script
    - Includes internal validation against author's pre-computed data

What it does:
    1. Loads ER+ IMC dataset
    2. Builds spatial graphs per ROI using squidpy (radius-based neighbors)
    3. Counts neighbors by cell type for each cell
    4. Aggregates by patient and cell type
    5. Validates against author's CSVs (internal check)
    6. Generates 3 core visualizations (+ optional validation plot)

Inputs:
    - data/imc_ERpos_with_embeddings_and_full_names.h5ad
    - data/raw/cycif_tma/data/*NeighborhoodCounts*.csv (for validation only)

Outputs:
    - analysis_imc/ERpos_neighbor_composition_r{25,40,100}_coarse.csv (3 files)
    - analysis_imc/ERpos_neighbor_composition_r{25,40}_fine.csv (2 files)
    - data/OUTPUT_imc_ERpos_with_graphs.h5ad (with spatial graphs)
    - figures/step3_neighbors/*.png (3 core plots + optional validation plot)

Radii by Granularity:
    - COARSE (5 types): r=25, 40, 100 μm
    - FINE (20 types): r=25, 40 μm (r=100 excluded: too many types, sparse neighbors)

"""

import os
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.sparse import csr_matrix, lil_matrix
import warnings
warnings.filterwarnings('ignore')

import anndata as ad
import squidpy as sq

# ============================================================================
# CONFIGURATION
# ============================================================================

# Cell type categories
COARSE_TYPES = ["endothelial", "epithelial", "fibroblast", "immune", "stromal"]

# Radii by granularity
RADII_COARSE = [25, 40, 100]
RADII_FINE = [25, 40]

# Input/Output paths
INPUT_H5AD = "data/imc_ERpos_with_embeddings_and_full_names.h5ad"
OUTPUT_H5AD = "data/OUTPUT_imc_ERpos_with_graphs.h5ad"
OUTPUT_DIR = Path("analysis_imc")
FIGURES_DIR = Path("figures/step3_neighbors")

# Matplotlib settings
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'

# Labels to drop from fine-grained heatmaps
UNCLASSIFIED_LABELS = {
    "20",
    "Unclassified",
    "Unclassified cells",
    "unclassified cells",
}


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def print_section(title):
    """Print formatted section header"""
    print(f"\n{'=' * 80}")
    print(f"{title}")
    print(f"{'=' * 80}")


# ============================================================================
# CORE FUNCTIONS: SPATIAL GRAPH CONSTRUCTION
# ============================================================================

def build_spatial_graph_per_roi(adata, radius, verbose=True):
    """
    Build spatial graph per ROI to avoid cross-ROI edges.

    Args:
        adata: AnnData with 'slide_scene' column and spatial coordinates
        radius: Neighborhood radius in microns
        verbose: Print progress updates

    Returns:
        Sparse connectivity matrix (n_cells x n_cells)
    """
    n_cells = adata.n_obs
    connectivity = lil_matrix((n_cells, n_cells))

    scenes = adata.obs['slide_scene'].unique()
    total_scenes = len(scenes)

    if verbose:
        print(f"   Processing {total_scenes} ROIs with radius={radius}μm...")

    for idx, scene in enumerate(scenes, 1):
        # Get indices for this ROI
        roi_mask = adata.obs['slide_scene'] == scene
        roi_indices = np.where(roi_mask)[0]

        if len(roi_indices) < 2:
            continue

        # Subset adata for this ROI
        adata_roi = adata[roi_mask].copy()

        # Build spatial graph with squidpy
        sq.gr.spatial_neighbors(
            adata_roi,
            coord_type="generic",
            radius=radius,
            delaunay=False,
            spatial_key='spatial'
        )

        # Extract connectivity
        roi_conn = adata_roi.obsp['spatial_connectivities']

        # Map back to global indices
        for i_local in range(len(roi_indices)):
            i_global = roi_indices[i_local]
            neighbors_local = roi_conn[i_local].nonzero()[1]

            for j_local in neighbors_local:
                j_global = roi_indices[j_local]
                connectivity[i_global, j_global] = 1

        # Progress update
        if verbose and (idx % 50 == 0 or idx == total_scenes):
            print(f"      ROI {idx}/{total_scenes} ({100*idx/total_scenes:.1f}%)", flush=True)

    return csr_matrix(connectivity)


def count_neighbors_by_celltype(adata, conn_matrix, celltype_col, celltypes):
    """
    Count neighbors of each cell type for every cell.

    Args:
        adata: AnnData object
        conn_matrix: Sparse connectivity matrix
        celltype_col: Column name for cell type annotation
        celltypes: List of cell types to count

    Returns:
        DataFrame with neighbor counts per cell
    """
    neighbor_counts = {}

    for ct in celltypes:
        ct_mask = (adata.obs[celltype_col] == ct).values
        counts = np.array(conn_matrix[:, ct_mask].sum(axis=1)).flatten()
        neighbor_counts[ct] = counts

    df = pd.DataFrame(neighbor_counts, index=adata.obs.index)
    df['Patient'] = adata.obs['Patient'].values
    df['celltype'] = adata.obs[celltype_col].values
    df['slide_scene'] = adata.obs['slide_scene'].values

    return df


# ============================================================================
# VALIDATION
# ============================================================================

def validate_against_authors(squidpy_df, author_csv_path, radius, celltypes):
    """
    Compare squidpy-computed counts against author's pre-computed CSVs.

    Args:
        squidpy_df: DataFrame with squidpy neighbor counts (cell-level)
        author_csv_path: Path to author's CSV file
        radius: Radius used
        celltypes: List of cell types to compare

    Returns:
        Dictionary with correlation metrics
    """
    if not Path(author_csv_path).exists():
        print(f"   ⚠ Author CSV not found: {author_csv_path}")
        return None

    print(f"   Validating squidpy vs author's data (r={radius})...")

    # Load author's data
    author_df = pd.read_csv(author_csv_path, low_memory=False).rename(
        columns={"Unnamed: 0": "cell_id"}
    )

    # Align by cell_id
    squidpy_df_indexed = squidpy_df.reset_index()
    squidpy_df_indexed.columns = ['cell_id'] + list(squidpy_df.columns)

    merged = squidpy_df_indexed.merge(
        author_df[['cell_id'] + celltypes],
        on='cell_id',
        how='inner',
        suffixes=('_squidpy', '_author')
    )

    if len(merged) == 0:
        print("   ⚠ No common cells found for validation")
        return None

    # Calculate correlations
    correlations = {}
    for ct in celltypes:
        sq_col = f'{ct}_squidpy' if f'{ct}_squidpy' in merged.columns else ct
        auth_col = f'{ct}_author' if f'{ct}_author' in merged.columns else ct

        # Handle case where suffix might not be added if column names don't conflict
        if ct in merged.columns and f'{ct}_author' not in merged.columns:
            # Need to check which is which
            continue

        if sq_col in merged.columns and auth_col in merged.columns:
            corr = merged[sq_col].corr(merged[auth_col])
            mae = (merged[sq_col] - merged[auth_col]).abs().mean()
            correlations[ct] = {'correlation': corr, 'MAE': mae}

    # Print summary
    print(f"\n   Validation Results (r={radius}):")
    print(f"   {'Cell Type':<15} {'Correlation':>12} {'MAE':>10}")
    print(f"   {'-'*40}")
    for ct, metrics in correlations.items():
        status = '✓' if metrics['correlation'] > 0.95 else '⚠'
        print(f"   {status} {ct:<14} {metrics['correlation']:>11.4f} {metrics['MAE']:>9.2f}")

    avg_corr = np.mean([m['correlation'] for m in correlations.values()])
    print(f"   {'-'*40}")
    print(f"   {'Average':<15} {avg_corr:>11.4f}")
    print(f"   {'Status':<15} {'PASS ✓' if avg_corr > 0.95 else 'REVIEW ⚠':>12}\n")

    return correlations


# ============================================================================
# MAIN PROCESSING PIPELINE
# ============================================================================

def process_granularity(adata, granularity, radii, celltype_col, celltypes):
    """
    Process a specific granularity (coarse or fine).

    Args:
        adata: AnnData object
        granularity: 'coarse' or 'fine'
        radii: List of radii to process
        celltype_col: Column name for cell type
        celltypes: List of cell types

    Returns:
        Dictionary mapping radius -> (connectivity_matrix, neighbor_df)
    """
    print_section(f"PROCESSING {granularity.upper()} GRANULARITY")
    print(f"Cell types: {len(celltypes)} ({', '.join(celltypes[:3])}...)")
    print(f"Radii: {radii} μm\n")

    results = {}

    for radius in radii:
        t0 = time.time()
        print(f"  [{granularity}] Radius = {radius} μm")

        # Build spatial graph
        conn_matrix = build_spatial_graph_per_roi(adata, radius, verbose=True)

        # Count neighbors
        print(f"   Counting neighbors by cell type...")
        neighbor_df = count_neighbors_by_celltype(adata, conn_matrix, celltype_col, celltypes)

        # Aggregate by patient
        print(f"   Aggregating by patient...")
        agg_patient_celltype = neighbor_df.groupby(['Patient', 'celltype'])[celltypes].mean()
        agg_patient = neighbor_df.groupby('Patient')[celltypes].mean()

        # Save CSVs
        output_file_pc = OUTPUT_DIR / f"ERpos_neighbor_composition_r{radius}_{granularity}.csv"
        agg_patient.to_csv(output_file_pc)
        print(f"   ✓ Saved: {output_file_pc}")

        # Store results
        results[radius] = {
            'connectivity': conn_matrix,
            'neighbor_df': neighbor_df,
            'agg_patient': agg_patient,
            'agg_patient_celltype': agg_patient_celltype
        }

        elapsed = time.time() - t0
        print(f"   Time: {elapsed:.1f}s\n")

    return results


def main():
    print_section("STEP 3: SPATIAL NEIGHBOR GRAPHS WITH SQUIDPY")
    print("\n📊 Computing neighborhood compositions with squidpy:")
    print("   • COARSE (5 types): r=25, 40, 100 μm")
    print("   • FINE (20 types): r=25, 40 μm")
    print("\nMethod: Graph per ROI to avoid cross-ROI edges")

    t0_total = time.time()

    # Create output directories
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # ========================================================================
    # 1. LOAD DATA
    # ========================================================================
    print_section("LOADING DATA")
    print(f"Input: {INPUT_H5AD}")

    adata = ad.read_h5ad(INPUT_H5AD)
    print(f"✓ Loaded: {adata.n_obs:,} cells × {adata.n_vars} markers")
    print(f"  Patients: {adata.obs['Patient'].nunique()}")
    print(f"  ROIs: {adata.obs['slide_scene'].nunique()}")

    # Ensure categorical
    adata.obs['leidencelltype5'] = adata.obs['leidencelltype5'].astype('category')
    adata.obs['leiden'] = adata.obs['leiden'].astype('category')

    # Get cell type lists
    fine_types = sorted(adata.obs['leiden'].cat.categories.tolist())
    fine_type_map = None
    if 'celltype_full' in adata.obs.columns:
        mapping_df = adata.obs[['leiden', 'celltype_full']].drop_duplicates()
        fine_type_map = dict(zip(mapping_df['leiden'].astype(str),
                                 mapping_df['celltype_full'].astype(str)))

    print(f"\nCell type breakdown:")
    print(f"  COARSE: {len(COARSE_TYPES)} types")
    print(f"  FINE: {len(fine_types)} types")

    # ========================================================================
    # 2. PROCESS COARSE GRANULARITY (5 types, r=25,40,100)
    # ========================================================================
    results_coarse = process_granularity(
        adata,
        granularity='coarse',
        radii=RADII_COARSE,
        celltype_col='leidencelltype5',
        celltypes=COARSE_TYPES
    )

    # ========================================================================
    # 3. PROCESS FINE GRANULARITY (20 types, r=25,40)
    # ========================================================================
    results_fine = process_granularity(
        adata,
        granularity='fine',
        radii=RADII_FINE,
        celltype_col='leiden',
        celltypes=fine_types
    )

    # ========================================================================
    # 4. VALIDATION AGAINST AUTHORS (COARSE ONLY)
    # ========================================================================
    print_section("VALIDATION: SQUIDPY VS AUTHORS")

    validation_results = {}

    for radius in [25, 40]:
        author_csv = f"data/raw/cycif_tma/data/20220420_JP-TMAs_IMC-TMAs_MIBI_leidencelltype5_NeighborhoodCounts_r{radius}.csv"

        if radius in results_coarse:
            neighbor_df = results_coarse[radius]['neighbor_df']
            correlations = validate_against_authors(
                neighbor_df,
                author_csv,
                radius,
                COARSE_TYPES
            )
            if correlations:
                validation_results[radius] = correlations

    # ========================================================================
    # 5. SAVE H5AD WITH SPATIAL GRAPHS
    # ========================================================================
    print_section("SAVING H5AD WITH SPATIAL GRAPHS")

    # Add all connectivity matrices to adata.obsp
    for radius in RADII_COARSE:
        if radius in results_coarse:
            key = f'spatial_connectivities_r{radius}'
            adata.obsp[key] = results_coarse[radius]['connectivity']
            print(f"  ✓ Added: adata.obsp['{key}']")

    # Fine granularity uses same connectivity as coarse for shared radii
    # (connectivity is independent of cell type annotation)

    # Save
    print(f"\n  Saving: {OUTPUT_H5AD}")
    adata.write_h5ad(OUTPUT_H5AD)
    print(f"  ✓ Saved H5AD with spatial graphs")

    # ========================================================================
    # 6. GENERATE VISUALIZATIONS
    # ========================================================================
    generate_visualizations(results_coarse, results_fine, validation_results, fine_type_map=fine_type_map)

    # ========================================================================
    # SUMMARY
    # ========================================================================
    elapsed_total = time.time() - t0_total
    print_section("STEP 3 COMPLETE")
    print(f"\n⏱ Total execution time: {elapsed_total/60:.1f} minutes")
    print("\n✅ Generated outputs:")
    print(f"   • {len(RADII_COARSE)} coarse CSV files (r=25, 40, 100)")
    print(f"   • {len(RADII_FINE)} fine CSV files (r=25, 40)")
    print(f"   • 1 H5AD with spatial graphs: {OUTPUT_H5AD}")
    print(f"   • 3-4 PNG figures (300 DPI): {FIGURES_DIR}/")
    print("\n📊 Key metrics:")
    print(f"   • Cells processed: {adata.n_obs:,}")
    print(f"   • Patients: {adata.obs['Patient'].nunique()}")
    print(f"   • ROIs: {adata.obs['slide_scene'].nunique()}")

    if validation_results:
        avg_corr = np.mean([
            np.mean([m['correlation'] for m in corrs.values()])
            for corrs in validation_results.values()
        ])
        print(f"\n✓ Validation: Avg correlation with authors = {avg_corr:.4f}")

    print("\n🎯 Next step: Step 4 (Cell-cell interactions)")
    print("=" * 80)


# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

def generate_visualizations(results_coarse, results_fine, validation_results, fine_type_map=None):
    """Generate all visualization plots"""
    print_section("GENERATING VISUALIZATIONS (3 core plots)")

    # Plot 1: Neighbor composition heatmap (coarse, r=25)
    plot_neighbor_composition_heatmap(results_coarse, granularity='coarse', radius=25)

    # Plot 2: Neighbor composition heatmap (fine, r=25)
    plot_neighbor_composition_heatmap(
        results_fine,
        granularity='fine',
        radius=25,
        celltype_map=fine_type_map,
        drop_labels=UNCLASSIFIED_LABELS,
    )

    # Plot 3: Radius comparison (coarse only)
    plot_radius_comparison_coarse(results_coarse)

    # Plot 4: Validation plot (if validation was performed)
    if validation_results:
        plot_validation(validation_results)

    print("\n✅ All core visualizations complete (3/3)")


def plot_neighbor_composition_heatmap(results, granularity, radius, celltype_map=None, drop_labels=None):
    """
    Plot heatmap of mean neighbors by source cell type.
    Shows spatial architecture: which cell types are near each other.
    """
    print(f"\n  [1/{granularity}] Neighbor Composition Heatmap (r={radius})...")

    if radius not in results:
        print(f"     ⚠ No data for radius {radius}")
        return

    agg_patient_celltype = results[radius]['agg_patient_celltype']

    # Compute mean across patients
    mean_neighbors = agg_patient_celltype.groupby('celltype').mean()

    # Apply label mapping (fine-grained heatmap uses full names)
    if celltype_map:
        mean_neighbors = mean_neighbors.rename(index=celltype_map, columns=celltype_map)

    # Drop unclassified labels from fine heatmap
    if drop_labels:
        drop_set = {label.strip() for label in drop_labels}
        mean_neighbors = mean_neighbors.loc[
            [idx for idx in mean_neighbors.index
             if str(idx).strip() not in drop_set and 'unclassified' not in str(idx).lower()],
            [col for col in mean_neighbors.columns
             if str(col).strip() not in drop_set and 'unclassified' not in str(col).lower()]
        ]

    # Create heatmap
    fig, ax = plt.subplots(figsize=(14, 12) if granularity == 'fine' else (10, 8), dpi=300)

    sns.heatmap(mean_neighbors, annot=True if granularity == 'coarse' else False,
                fmt='.1f' if granularity == 'coarse' else '',
                cmap='YlOrRd',
                cbar_kws={'label': 'Mean Neighbor Count', 'shrink': 0.5},
                linewidths=0.5, linecolor='white', ax=ax)

    ax.set_title(f'Spatial Neighbor Composition ({granularity.upper()}, r={radius}μm)',
                 fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel('Neighbor Cell Type', fontsize=12, fontweight='bold')
    ax.set_ylabel('Source Cell Type', fontsize=12, fontweight='bold')

    if granularity == 'fine':
        plt.xticks(rotation=45, ha='right', fontsize=7)
        plt.yticks(rotation=0, fontsize=7)

    plt.tight_layout()
    filename = f'01_neighbor_composition_heatmap_{granularity}_r{radius}.png'
    plt.savefig(FIGURES_DIR / filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"     ✓ Saved: {filename}")


def plot_radius_comparison_coarse(results_coarse):
    """
    Plot box plots comparing neighbor counts across radii (coarse only).
    Teaches concept of spatial scale.
    """
    print("\n  [2] Radius Comparison Box Plots (coarse)...")

    # Prepare data
    data_list = []
    for radius in [25, 40, 100]:
        if radius not in results_coarse:
            continue

        df = results_coarse[radius]['agg_patient']
        for col in COARSE_TYPES:
            for val in df[col]:
                data_list.append({
                    'Radius': f'r={radius}μm',
                    'Neighbor Type': col,
                    'Count': val
                })

    if not data_list:
        print("     ⚠ No data for radius comparison")
        return

    plot_df = pd.DataFrame(data_list)

    # Create faceted box plots
    fig, axes = plt.subplots(1, 5, figsize=(20, 5), dpi=300)

    for i, neighbor_type in enumerate(COARSE_TYPES):
        subset = plot_df[plot_df['Neighbor Type'] == neighbor_type]

        sns.boxplot(data=subset, x='Radius', y='Count', ax=axes[i],
                    palette='Set2', width=0.6)

        axes[i].set_title(f'{neighbor_type.capitalize()} Neighbors',
                          fontsize=12, fontweight='bold')
        axes[i].set_xlabel('Spatial Scale', fontsize=10, fontweight='bold')
        axes[i].set_ylabel('Mean Neighbor Count', fontsize=10, fontweight='bold')
        axes[i].grid(axis='y', alpha=0.3)
        axes[i].tick_params(axis='x', rotation=45)

    plt.suptitle('Neighborhood Scale Comparison: r=25, 40, 100 μm (COARSE)',
                 fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / '02_radius_comparison_boxplots.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("     ✓ Saved: 02_radius_comparison_boxplots.png")


def plot_tumor_immune_scatter(results_coarse):
    """
    Plot scatter of tumor fraction vs immune neighbors (r=25).
    Key question: Do tumor-rich patients have poor immune infiltration?
    """
    print("\n  [3] Tumor-Immune Relationship Scatter...")

    if 25 not in results_coarse:
        print("     ⚠ No data for r=25")
        return

    # Load composition fractions
    agg_patient = results_coarse[25]['agg_patient']
    agg_patient_celltype = results_coarse[25]['agg_patient_celltype']

    # Compute epithelial fraction
    neighbor_df = results_coarse[25]['neighbor_df']
    epithelial_frac = neighbor_df.groupby('Patient')['celltype'].apply(
        lambda x: (x == 'epithelial').sum() / len(x)
    )

    # Get immune neighbors for epithelial cells
    epithelial_cells = agg_patient_celltype.loc[
        agg_patient_celltype.index.get_level_values('celltype') == 'epithelial'
    ]
    tumor_immune = epithelial_cells[['immune']].reset_index()
    tumor_immune.columns = ['Patient', 'celltype', 'immune_neighbors']

    # Merge
    merged = pd.DataFrame({
        'Patient': epithelial_frac.index,
        'epithelial_frac': epithelial_frac.values
    }).merge(tumor_immune[['Patient', 'immune_neighbors']], on='Patient', how='inner')

    # Assign cohort
    merged['Cohort'] = merged['Patient'].apply(
        lambda x: 'Zürich' if str(x).startswith('Z') else 'Basel'
    )

    # Create scatter plot
    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)

    for cohort, color in [('Basel', 'steelblue'), ('Zürich', 'coral')]:
        subset = merged[merged['Cohort'] == cohort]
        ax.scatter(subset['epithelial_frac'] * 100, subset['immune_neighbors'],
                   c=color, label=cohort, s=80, alpha=0.7, edgecolor='black', linewidth=0.5)

    # Add trend line
    from scipy.stats import pearsonr
    x = merged['epithelial_frac'] * 100
    y = merged['immune_neighbors']
    z = np.polyfit(x, y, 1)
    p = np.poly1d(z)
    ax.plot(x, p(x), "k--", alpha=0.5, linewidth=2)

    # Compute correlation
    corr, pval = pearsonr(x, y)

    ax.set_xlabel('Tumor Cell Fraction (%)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Mean Immune Neighbors (r=25μm)', fontsize=13, fontweight='bold')
    ax.set_title(f'Tumor Burden vs Immune Infiltration\n(r={corr:.3f}, p={pval:.2e})',
                 fontsize=14, fontweight='bold', pad=15)
    ax.legend(title='Cohort', fontsize=11, title_fontsize=12, loc='upper right')
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / '03_tumor_immune_scatter.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("     ✓ Saved: 03_tumor_immune_scatter.png")


def plot_neighbor_violins(results_coarse):
    """
    Plot violin plots showing distribution of neighbors for top 20 patients.
    Shows intra-patient heterogeneity.
    """
    print("\n  [4] Neighbor Distribution Violins...")

    if 25 not in results_coarse:
        print("     ⚠ No data for r=25")
        return

    neighbor_df = results_coarse[25]['neighbor_df']

    # Get top 20 patients by cell count
    top_patients = neighbor_df['Patient'].value_counts().head(20).index.tolist()
    subset = neighbor_df[neighbor_df['Patient'].isin(top_patients)].copy()

    # Create violin plots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12), dpi=300)
    axes = axes.flatten()

    for i, neighbor_type in enumerate(COARSE_TYPES):
        ax = axes[i]

        # Filter to reasonable range for visualization
        plot_data = subset[subset[neighbor_type] <= subset[neighbor_type].quantile(0.95)]

        sns.violinplot(data=plot_data, x='Patient', y=neighbor_type, ax=ax,
                       palette='Set3', inner='quartile', cut=0)

        ax.set_title(f'{neighbor_type.capitalize()} Neighbors',
                     fontsize=12, fontweight='bold')
        ax.set_xlabel('Patient (Top 20)', fontsize=10, fontweight='bold')
        ax.set_ylabel('Neighbor Count (r=25μm)', fontsize=10, fontweight='bold')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
        ax.grid(axis='y', alpha=0.3)

    # Remove extra subplot
    fig.delaxes(axes[5])

    plt.suptitle('Neighbor Count Distributions by Patient (r=25μm, COARSE)',
                 fontsize=15, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / '04_neighbor_distribution_violins.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("     ✓ Saved: 04_neighbor_distribution_violins.png")


def plot_validation(validation_results):
    """
    Plot validation comparison: squidpy vs authors.
    Shows correlation for each cell type and radius.
    """
    print("\n  [5] Validation: Squidpy vs Authors...")

    # Prepare data
    data = []
    for radius, corrs in validation_results.items():
        for celltype, metrics in corrs.items():
            data.append({
                'Radius': f'r={radius}μm',
                'Cell Type': celltype,
                'Correlation': metrics['correlation'],
                'MAE': metrics['MAE']
            })

    if not data:
        print("     ⚠ No validation data")
        return

    df = pd.DataFrame(data)

    # Create figure with 2 subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), dpi=300)

    # Plot 1: Correlation barplot
    df_pivot = df.pivot(index='Cell Type', columns='Radius', values='Correlation')
    df_pivot.plot(kind='bar', ax=ax1, color=['steelblue', 'coral'], width=0.7)
    ax1.axhline(y=0.95, color='red', linestyle='--', linewidth=2, alpha=0.5, label='Target (0.95)')
    ax1.set_xlabel('Cell Type', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Pearson Correlation', fontsize=12, fontweight='bold')
    ax1.set_title('Squidpy vs Authors: Correlation', fontsize=13, fontweight='bold')
    ax1.legend(title='Radius', fontsize=10)
    ax1.grid(axis='y', alpha=0.3)
    ax1.set_ylim(0.9, 1.0)
    plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')

    # Plot 2: MAE barplot
    df_pivot_mae = df.pivot(index='Cell Type', columns='Radius', values='MAE')
    df_pivot_mae.plot(kind='bar', ax=ax2, color=['steelblue', 'coral'], width=0.7)
    ax2.set_xlabel('Cell Type', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Mean Absolute Error', fontsize=12, fontweight='bold')
    ax2.set_title('Squidpy vs Authors: MAE', fontsize=13, fontweight='bold')
    ax2.legend(title='Radius', fontsize=10)
    ax2.grid(axis='y', alpha=0.3)
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')

    plt.suptitle('Validation: Squidpy vs Authors (Internal Check)',
                 fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / '05_validation_squidpy_vs_authors.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("     ✓ Saved: 05_validation_squidpy_vs_authors.png")


if __name__ == "__main__":
    main()
