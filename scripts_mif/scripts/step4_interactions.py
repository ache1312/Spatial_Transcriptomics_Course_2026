#!/usr/bin/env python3
"""
STEP 4: Cell-Cell Spatial Interactions (Unified)

Purpose:
    Compute DIRECTED cell-cell interaction features from spatial graphs.

Key Concept: ASYMMETRIC/DIRECTED neighbors
    - epithelial → immune ≠ immune → epithelial
    - "How many immune neighbors do epithelial cells have?" (directed)
    - Captures critical biology: infiltration, exclusion, co-localization

What it does:
    1. Loads spatial graphs from Step 3 (pre-computed with squidpy)
    2. Computes directed neighbor features for COARSE and FINE granularities
    3. Aggregates by patient
    4. Generates comprehensive visualizations

Inputs:
    - data/OUTPUT_imc_ERpos_with_graphs.h5ad (with spatial graphs from Step 3)

Outputs (CSV):
    - analysis_imc/ERpos_directed_neighbors_r{25,40,100}_coarse.csv (25 features each)
    - analysis_imc/ERpos_directed_neighbors_r{25,40}_fine.csv (400 features each)

Outputs (Figures):
    - figures/step4_interactions/*.png (8 visualization plots)

Radii by Granularity:
    - COARSE (5 types): r=25, 40, 100 μm
    - FINE (20 types): r=25, 40 μm

Time: ~5-10 minutes
"""

import os
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from pathlib import Path
import warnings
import textwrap
warnings.filterwarnings('ignore')

import anndata as ad

# ============================================================================
# CONFIGURATION
# ============================================================================

# Cell type categories
COARSE_TYPES = ["endothelial", "epithelial", "fibroblast", "immune", "stromal"]

# Radii by granularity
RADII_COARSE = [25, 40, 100]
RADII_FINE = [25, 40]

# Input/Output paths
INPUT_H5AD = "data/OUTPUT_imc_ERpos_with_graphs.h5ad"
OUTPUT_DIR = Path("analysis_imc")
FIGURES_DIR = Path("figures/step4_interactions")

# Matplotlib settings
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def print_section(title, level=1):
    """Print formatted section header"""
    if level == 1:
        print(f"\n{'=' * 80}")
        print(f"{title}")
        print(f"{'=' * 80}")
    else:
        print(f"\n{'─' * 70}")
        print(f"  {title}")
        print(f"{'─' * 70}")


# ============================================================================
# CORE FUNCTIONS: DIRECTED FEATURES
# ============================================================================

def compute_directed_features(adata, conn_key, celltype_col, celltypes, use_full_names=False):
    """
    Compute directed neighbor features: origin_type -> target_type

    Args:
        adata: AnnData with spatial graph in obsp[conn_key]
        conn_key: Key for connectivity matrix (e.g., 'spatial_connectivities_r25')
        celltype_col: Column for cell types
        celltypes: List of cell types to analyze
        use_full_names: If True, use full cell type names for features

    Returns:
        DataFrame with patients × directed_features
    """
    if conn_key not in adata.obsp:
        print(f"  ⚠ Warning: {conn_key} not found in adata.obsp")
        return None

    conn_matrix = adata.obsp[conn_key]

    # Get full names if available
    if use_full_names and 'celltype_full' in adata.obs.columns:
        celltype_full_map = adata.obs.groupby(celltype_col)['celltype_full'].first().to_dict()
    else:
        celltype_full_map = {ct: ct for ct in celltypes}

    # Count neighbors by type for each cell
    neighbor_counts = {}
    for target_ct in celltypes:
        target_mask = (adata.obs[celltype_col] == target_ct).values
        counts = np.array(conn_matrix[:, target_mask].sum(axis=1)).flatten()
        target_name = celltype_full_map.get(target_ct, target_ct)
        neighbor_counts[target_name] = counts

    # Build cell-level DataFrame
    df = pd.DataFrame(neighbor_counts, index=adata.obs.index)
    df['Patient'] = adata.obs['Patient'].values
    df['celltype'] = adata.obs[celltype_col].values

    # Compute directed features per patient
    feature_dfs = []
    target_names = list(neighbor_counts.keys())

    for origin_ct in celltypes:
        origin_cells = df[df['celltype'] == origin_ct]
        if len(origin_cells) == 0:
            continue

        # Average neighbors per target type
        agg = origin_cells.groupby('Patient')[target_names].mean()

        # Add prefix: origin_to_target
        origin_name = celltype_full_map.get(origin_ct, origin_ct)
        agg = agg.add_prefix(f"{origin_name}_to_")

        feature_dfs.append(agg)

    if not feature_dfs:
        return None

    result = pd.concat(feature_dfs, axis=1)
    return result


# ============================================================================
# MAIN PROCESSING PIPELINE
# ============================================================================

def process_granularity(adata, granularity, radii, celltype_col, celltypes):
    """
    Process a specific granularity (coarse or fine).

    Args:
        adata: AnnData object with spatial graphs
        granularity: 'coarse' or 'fine'
        radii: List of radii to process
        celltype_col: Column name for cell type
        celltypes: List of cell types

    Returns:
        Dictionary mapping radius -> directed_features_df
    """
    print_section(f"PROCESSING {granularity.upper()} GRANULARITY", level=2)
    print(f"  Cell types: {len(celltypes)}")
    print(f"  Radii: {radii} μm\n")

    results = {}

    for radius in radii:
        t0 = time.time()
        conn_key = f'spatial_connectivities_r{radius}'

        print(f"  [{granularity}] Radius = {radius} μm")

        # Compute directed features
        use_full_names = (granularity == 'fine')
        df_directed = compute_directed_features(
            adata, conn_key, celltype_col, celltypes, use_full_names=use_full_names
        )

        if df_directed is not None:
            # Save CSV
            output_file = OUTPUT_DIR / f"ERpos_directed_neighbors_r{radius}_{granularity}.csv"
            df_directed.to_csv(output_file)
            print(f"    ✓ Saved: {output_file}")
            print(f"      Shape: {df_directed.shape[0]} patients × {df_directed.shape[1]} features")

            results[radius] = df_directed
        else:
            print(f"    ⚠ No data for radius {radius}")

        elapsed = time.time() - t0
        print(f"    Time: {elapsed:.1f}s\n")

    return results


def main():
    print_section("STEP 4: CELL-CELL SPATIAL INTERACTIONS (UNIFIED)")
    print("\n📊 Computing directed interaction features:")
    print("   • COARSE (5 types): r=25, 40, 100 μm")
    print("   • FINE (20 types): r=25, 40 μm")
    print("\nMethod: Uses pre-computed spatial graphs from Step 3")

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

    # Verify spatial graphs exist
    graph_keys = [k for k in adata.obsp.keys() if k.startswith('spatial_connectivities')]
    print(f"\n  Spatial graphs found: {len(graph_keys)}")
    for key in sorted(graph_keys):
        print(f"    - {key}")

    if len(graph_keys) == 0:
        print("\n  ⚠ ERROR: No spatial graphs found!")
        print("     Please run Step 3 first to generate spatial graphs.")
        return

    # Verify celltype_full exists (should come from Step 3's input h5ad)
    if 'celltype_full' not in adata.obs.columns:
        print("\n  ⚠ ERROR: celltype_full not found in h5ad!")
        print("     Please ensure you're using imc_ERpos_with_embeddings_and_full_names.h5ad")
        return

    print(f"  ✓ Using celltype_full from h5ad (no CSV mapping needed)")

    # Get cell type lists
    fine_types = sorted(adata.obs['leiden'].cat.categories.tolist())

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
    # 4. GENERATE VISUALIZATIONS
    # ========================================================================
    generate_visualizations(adata, results_coarse, results_fine)

    # ========================================================================
    # SUMMARY
    # ========================================================================
    elapsed_total = time.time() - t0_total
    print_section("STEP 4 COMPLETE")
    print(f"\n⏱ Total execution time: {elapsed_total/60:.1f} minutes")
    print("\n✅ Generated outputs:")
    print(f"   • {len(RADII_COARSE)} coarse CSV files (r=25, 40, 100)")
    print(f"   • {len(RADII_FINE)} fine CSV files (r=25, 40)")
    print(f"   • 8 PNG figures (300 DPI): {FIGURES_DIR}/")
    print("\n📊 Key metrics:")
    print(f"   • Cells analyzed: {adata.n_obs:,}")
    print(f"   • Patients: {adata.obs['Patient'].nunique()}")
    print(f"   • COARSE features per patient: 25 (5×5 directed pairs)")
    print(f"   • FINE features per patient: {len(fine_types)**2} ({len(fine_types)}×{len(fine_types)} directed pairs)")
    print("\n🎯 Next step: Step 5 (Neighborhood Enrichment)")
    print("=" * 80)


# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

def generate_visualizations(adata, results_coarse, results_fine):
    """Generate all visualization plots"""
    print_section("GENERATING VISUALIZATIONS (8 plots)")

    # Plot 1: Asymmetry heatmap (coarse, r=25)
    plot_asymmetry_heatmap(results_coarse, granularity='coarse', radius=25)

    # Plot 2: Asymmetry heatmap (fine, r=25) - simplified version
    plot_asymmetry_heatmap(results_fine, granularity='fine', radius=25)

    # Plot 3: Top 20 interactions (coarse, r=25)
    plot_top_interactions(results_coarse, granularity='coarse', radius=25, top_n=20)

    # Plot 4: Top 20 interactions (fine, r=25)
    plot_top_interactions(results_fine, granularity='fine', radius=25, top_n=20)

    # Plot 5: Tumor-immune infiltration scatter (coarse)
    plot_tumor_immune_infiltration(results_coarse)

    # Plot 6: Patient variability (coarse)
    plot_patient_variability(results_coarse, granularity='coarse', radius=25)

    # Plot 7: Interaction network (fine, r=25)
    plot_interaction_network(adata, results_fine, radius=25)

    # Plot 8: Radius comparison (coarse only)
    plot_radius_comparison(results_coarse)

    print("\n✅ All visualizations complete (8/8)")


def plot_asymmetry_heatmap(results, granularity, radius):
    """
    Plot directional asymmetry heatmap.
    Shows origin→target matrix and asymmetry (A - A.T)
    """
    print(f"\n  [1/{granularity}] Asymmetry Heatmap (r={radius})...")

    if radius not in results:
        print(f"     ⚠ No data for radius {radius}")
        return

    df = results[radius]

    if granularity == 'coarse':
        # Build 5×5 matrix
        matrix = np.zeros((5, 5))
        labels = COARSE_TYPES

        for i, origin in enumerate(labels):
            for j, target in enumerate(labels):
                col = f"{origin}_to_{target}"
                if col in df.columns:
                    matrix[i, j] = df[col].mean()

        # Create figure with 2 subplots
        fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=300)

        # Left: Raw directed matrix
        ax1 = axes[0]
        sns.heatmap(matrix, annot=True, fmt='.1f', cmap='YlOrRd',
                    xticklabels=[l.capitalize() for l in labels],
                    yticklabels=[l.capitalize() for l in labels],
                    cbar_kws={'label': 'Mean Neighbor Count', 'shrink': 0.5},
                    linewidths=0.5, linecolor='white', ax=ax1)
        ax1.set_xlabel('Neighbor Type (Target)', fontsize=11, fontweight='bold')
        ax1.set_ylabel('Cell Type (Origin)', fontsize=11, fontweight='bold')
        ax1.set_title(f'Directed Neighborhood Matrix (r={radius}μm)\nRows=Origin, Columns=Target',
                      fontsize=12, fontweight='bold')

        # Right: Asymmetry (A - A.T)
        ax2 = axes[1]
        asymmetry = matrix - matrix.T
        mask = np.eye(5, dtype=bool)

        sns.heatmap(asymmetry, annot=True, fmt='.1f', cmap='RdBu_r', center=0,
                    xticklabels=[l.capitalize() for l in labels],
                    yticklabels=[l.capitalize() for l in labels],
                    cbar_kws={'label': 'Asymmetry (Origin→Target - Target→Origin)', 'shrink': 0.5},
                    linewidths=0.5, linecolor='white', mask=mask, ax=ax2)
        ax2.set_xlabel('Target', fontsize=11, fontweight='bold')
        ax2.set_ylabel('Origin', fontsize=11, fontweight='bold')
        ax2.set_title('Directional Asymmetry\nPositive = Origin has more Target neighbors',
                      fontsize=12, fontweight='bold')

        plt.tight_layout()
        filename = f'01_asymmetry_heatmap_{granularity}_r{radius}.png'
    else:
        # Fine: Just show the raw matrix (no asymmetry, too large)
        # Get top 15 cell types for visualization
        all_features = df.mean().sort_values(ascending=False)
        top_pairs = all_features.head(15).index.tolist()

        # Extract origin and target from feature names
        origins = set()
        targets = set()
        for pair in top_pairs:
            if '_to_' in pair:
                origin, target = pair.split('_to_', 1)
                origins.add(origin)
                targets.add(target)

        # Build matrix for top types
        origins = sorted(origins)
        targets = sorted(targets)

        matrix = np.zeros((len(origins), len(targets)))
        for i, origin in enumerate(origins):
            for j, target in enumerate(targets):
                col = f"{origin}_to_{target}"
                if col in df.columns:
                    matrix[i, j] = df[col].mean()

        fig, ax = plt.subplots(figsize=(12, 10), dpi=300)

        sns.heatmap(matrix, annot=False, cmap='YlOrRd',
                    xticklabels=targets,
                    yticklabels=origins,
                    cbar_kws={'label': 'Mean Neighbor Count', 'shrink': 0.5},
                    linewidths=0.5, linecolor='white', ax=ax)
        ax.set_xlabel('Neighbor Type (Target)', fontsize=11, fontweight='bold')
        ax.set_ylabel('Cell Type (Origin)', fontsize=11, fontweight='bold')
        ax.set_title(f'Top Directed Interactions ({granularity.upper()}, r={radius}μm)',
                     fontsize=12, fontweight='bold')
        plt.xticks(rotation=45, ha='right', fontsize=8)
        plt.yticks(rotation=0, fontsize=8)
        plt.tight_layout()

        filename = f'01_asymmetry_heatmap_{granularity}_r{radius}.png'

    plt.savefig(FIGURES_DIR / filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"     ✓ Saved: {filename}")


def plot_top_interactions(results, granularity, radius, top_n=20):
    """
    Plot top N strongest directed interactions.
    """
    print(f"\n  [2/{granularity}] Top {top_n} Interactions (r={radius})...")

    if radius not in results:
        print(f"     ⚠ No data for radius {radius}")
        return

    df = results[radius]

    # Get mean of each feature across patients
    means = df.mean().sort_values(ascending=False).head(top_n)

    # Parse origin → target from column names
    labels = []
    wrap_width = 24 if granularity == 'fine' else 16
    for col in means.index:
        if '_to_' in col:
            origin, target = col.split('_to_', 1)
            if granularity == 'fine':
                origin_wrapped = "\n".join(
                    textwrap.wrap(origin, width=wrap_width, break_long_words=False, break_on_hyphens=False)
                )
                target_wrapped = "\n".join(
                    textwrap.wrap(target, width=wrap_width, break_long_words=False, break_on_hyphens=False)
                )
                labels.append(f"{origin_wrapped}\n→\n{target_wrapped}")
            else:
                labels.append(f"{origin} → {target}")
        else:
            labels.append(col)

    fig_height = max(8, 0.45 * len(means) + 4)
    fig_width = 12 if granularity == 'fine' else 10
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=300)

    y_pos = np.arange(len(means))
    ax.barh(y_pos, means.values, color='steelblue', edgecolor='white', linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8 if granularity == 'fine' else 9)
    ax.tick_params(axis='y', pad=6)
    ax.set_xlabel('Mean Neighbor Count', fontsize=11, fontweight='bold')
    ax.set_title(f'Top {top_n} Directed Interactions ({granularity.upper()}, r={radius}μm)\n'
                 f'Averaged across {df.shape[0]} patients',
                 fontsize=12, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    ax.invert_yaxis()

    # Add value labels
    max_val = means.values.max() if len(means) else 1.0
    x_offset = max_val * 0.02
    ax.set_xlim(0, max_val * 1.15)
    for i, val in enumerate(means.values):
        ax.text(val + x_offset, i, f'{val:.1f}', va='center', fontsize=8)

    plt.tight_layout()
    filename = f'02_top{top_n}_interactions_{granularity}_r{radius}.png'
    plt.savefig(FIGURES_DIR / filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"     ✓ Saved: {filename}")


def plot_tumor_immune_infiltration(results_coarse):
    """
    Scatter plot: epithelial→immune vs immune→epithelial.
    Critical for understanding tumor-immune spatial relationships.
    """
    print("\n  [3] Tumor-Immune Infiltration Scatter...")

    if 25 not in results_coarse:
        print("     ⚠ No data for r=25")
        return

    df = results_coarse[25]

    x_col = 'epithelial_to_immune'
    y_col = 'immune_to_epithelial'

    if x_col not in df.columns or y_col not in df.columns:
        print(f"     ⚠ Required columns not found: {x_col}, {y_col}")
        return

    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)

    x = df[x_col]
    y = df[y_col]

    # Color by cohort
    df_plot = df.copy()
    df_plot['Cohort'] = df_plot.index.astype(str).map(
        lambda p: 'Zürich' if str(p).startswith('Z') else 'Basel'
    )

    colors = {'Basel': '#3498db', 'Zürich': '#e74c3c'}

    for cohort in ['Basel', 'Zürich']:
        mask = df_plot['Cohort'] == cohort
        ax.scatter(df_plot.loc[mask, x_col], df_plot.loc[mask, y_col],
                   c=colors[cohort], label=cohort, alpha=0.6, s=50, edgecolors='white')

    # Diagonal (symmetry reference)
    max_val = max(x.max(), y.max()) * 1.1
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, label='Symmetry')

    # Correlation
    corr = x.corr(y)
    ax.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax.transAxes,
            fontsize=12, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.set_xlabel('Epithelial → Immune\n(Immune neighbors of tumor cells)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Immune → Epithelial\n(Tumor neighbors of immune cells)', fontsize=11, fontweight='bold')
    ax.set_title(f'Tumor-Immune Infiltration Directionality (r=25μm)\nEach point = 1 patient',
                 fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / '03_tumor_immune_infiltration.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("     ✓ Saved: 03_tumor_immune_infiltration.png")


def plot_patient_variability(results, granularity, radius):
    """
    Box plots showing patient-level variability for key features.
    """
    print(f"\n  [4] Patient Variability ({granularity})...")

    if radius not in results:
        print("     ⚠ No data for radius {radius}")
        return

    df = results[radius]

    if granularity == 'coarse':
        # Select key pairs
        key_pairs = [
            ('epithelial', 'epithelial'),
            ('epithelial', 'immune'),
            ('epithelial', 'fibroblast'),
            ('immune', 'epithelial'),
            ('immune', 'immune'),
            ('stromal', 'epithelial'),
        ]

        data_list = []
        for origin, target in key_pairs:
            col = f"{origin}_to_{target}"
            if col in df.columns:
                for val in df[col]:
                    data_list.append({
                        'Pair': f"{origin[:3].upper()}→{target[:3].upper()}",
                        'Count': val
                    })
    else:
        # For fine, just take top 8 features by mean
        top_features = df.mean().nlargest(8)

        data_list = []
        for col in top_features.index:
            for val in df[col]:
                # Shorten name for display
                short_name = col.replace('_to_', '→')[:25]
                data_list.append({
                    'Pair': short_name,
                    'Count': val
                })

    if not data_list:
        print("     ⚠ No data for variability plot")
        return

    plot_df = pd.DataFrame(data_list)

    fig, ax = plt.subplots(figsize=(12, 6), dpi=300)

    sns.boxplot(data=plot_df, x='Pair', y='Count', palette='Set2', width=0.6, ax=ax)
    sns.stripplot(data=plot_df, x='Pair', y='Count', color='black', alpha=0.3, size=3, ax=ax)

    ax.set_xlabel('Directed Pair (Origin → Target)', fontsize=11, fontweight='bold')
    ax.set_ylabel(f'Mean Neighbor Count (r={radius}μm)', fontsize=11, fontweight='bold')
    ax.set_title(f'Patient Variability in Directed Neighbors ({granularity.upper()})\nEach point = 1 patient',
                 fontsize=12, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)

    # Rotate labels if fine-grained
    if granularity == 'fine':
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)

    plt.tight_layout()
    filename = f'04_patient_variability_{granularity}_r{radius}.png'
    plt.savefig(FIGURES_DIR / filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"     ✓ Saved: {filename}")


def plot_interaction_network(adata, results_fine, radius):
    """
    Network graph of cell-cell interactions (fine granularity).
    Shows directed edges for interactions above threshold.
    """
    print(f"\n  [5] Interaction Network (fine, r={radius})...")

    if radius not in results_fine:
        print("     ⚠ No data for radius {radius}")
        return

    df = results_fine[radius]

    # Compute mean interaction strength for each pair
    interaction_means = df.mean()

    # Create network from significant interactions (threshold at mean > 1.0)
    G = nx.DiGraph()

    threshold = 1.0
    for col in interaction_means.index:
        if interaction_means[col] > threshold and '_to_' in col:
            origin, target = col.split('_to_', 1)
            G.add_edge(origin, target, weight=interaction_means[col])

    if len(G.nodes()) == 0:
        print(f"     ⚠ No interactions above threshold {threshold}")
        return

    fig, ax = plt.subplots(figsize=(14, 14), dpi=300)

    # Position nodes using spring layout
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

    # Node sizes proportional to cell abundance
    if 'celltype_full' in adata.obs.columns:
        celltype_counts = adata.obs['celltype_full'].value_counts()
    else:
        celltype_counts = adata.obs['leiden'].value_counts()

    node_sizes = [celltype_counts.get(node, 100) / 10 for node in G.nodes()]

    nx.draw_networkx_nodes(
        G, pos,
        node_size=node_sizes,
        node_color='lightblue',
        alpha=0.7,
        ax=ax
    )

    # Edge widths proportional to interaction strength
    edges = G.edges()
    weights = [G[u][v]['weight'] for u, v in edges]
    max_weight = max(weights) if weights else 1.0
    edge_widths = [3 * w / max_weight for w in weights]

    nx.draw_networkx_edges(
        G, pos,
        width=edge_widths,
        alpha=0.5,
        edge_color='gray',
        arrows=True,
        arrowsize=10,
        ax=ax
    )

    nx.draw_networkx_labels(
        G, pos,
        font_size=7,
        font_weight='bold',
        ax=ax
    )

    ax.set_title(
        f'Cell-Cell Interaction Network (r={radius}μm)\nEdges: mean neighbors > {threshold}',
        fontsize=13,
        fontweight='bold'
    )
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / f'05_interaction_network_r{radius}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"     ✓ Saved: 05_interaction_network_r{radius}.png")


def plot_radius_comparison(results_coarse):
    """
    Box plots comparing directed interactions across radii (coarse only).
    Shows how spatial scale affects interaction measurements.
    """
    print("\n  [6] Radius Comparison (coarse)...")

    # Prepare data
    data_list = []
    for radius in RADII_COARSE:
        if radius not in results_coarse:
            continue

        df = results_coarse[radius]

        # Select key interaction pairs
        key_pairs = ['epithelial_to_immune', 'immune_to_epithelial',
                     'epithelial_to_epithelial', 'immune_to_immune']

        for pair in key_pairs:
            if pair in df.columns:
                for val in df[pair]:
                    origin, target = pair.replace('_to_', ' → ').split(' → ')
                    data_list.append({
                        'Radius': f'r={radius}μm',
                        'Pair': f'{origin[:3].upper()}→{target[:3].upper()}',
                        'Count': val
                    })

    if not data_list:
        print("     ⚠ No data for radius comparison")
        return

    plot_df = pd.DataFrame(data_list)

    fig, ax = plt.subplots(figsize=(12, 6), dpi=300)

    sns.boxplot(data=plot_df, x='Pair', y='Count', hue='Radius',
                palette='Set2', ax=ax)

    ax.set_xlabel('Interaction Pair', fontsize=11, fontweight='bold')
    ax.set_ylabel('Mean Neighbor Count', fontsize=11, fontweight='bold')
    ax.set_title('Spatial Scale Comparison: Directed Interactions (COARSE)\nr=25, 40, 100 μm',
                 fontsize=12, fontweight='bold')
    ax.legend(title='Radius', fontsize=10)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / '06_radius_comparison_interactions.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("     ✓ Saved: 06_radius_comparison_interactions.png")


if __name__ == "__main__":
    main()
