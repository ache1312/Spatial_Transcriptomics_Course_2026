#!/usr/bin/env python3
"""
STEP 2: Embeddings & Quality Control Visualizations

Purpose:
    Generate comprehensive visualizations of PCA, UMAP, and QC metrics for workshop

What it does:
    1. Loads pre-computed embeddings (PCA, UMAP) from h5ad
    2. Generates 11 high-quality figures (300 DPI):
       - 3 PCA plots (scree, loadings, biplot)
       - 2 UMAP plots (cell types, cohort)
       - 2 Crossbar plots (composition by cohort, by patient)
       - 4 QC plots (cells per patient, cell area, cohort comparison, cell type counts)

Input:
    - data/imc_ERpos_with_embeddings_and_full_names.h5ad

Output:
    - figures/step2_embeddings_qc/*.png (11 figures at 300 DPI)

IMPORTANT - Data Preprocessing:
    All visualizations are based on preprocessed data:
    1. arcsinh transformation: Raw intensities transformed using arcsinh(x)
    2. PCA: Applied to arcsinh-transformed data (20 PCs, 99.9% variance)
    3. UMAP: Applied to first 20 PCs (2D embedding)

"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import time


def print_section(title):
    """Print formatted section header"""
    print(f"\n{'=' * 80}")
    print(f"{title}")
    print(f"{'=' * 80}")


def generate_pca_plots(adata, output_dir):
    """Generate 3 PCA visualization plots"""
    print_section("GENERATING PCA PLOTS (3 figures)")

    if 'X_pca' not in adata.obsm or 'pca' not in adata.uns:
        print("⚠️  PCA data not found, skipping PCA plots")
        return

    # -------------------------------------------------------------------------
    # Figure 1: PCA Scree Plot - Simple line plot
    # -------------------------------------------------------------------------
    print("\n[1/11] PCA Scree Plot...")

    variance_ratio = adata.uns['pca']['variance_ratio']
    n_components = len(variance_ratio)

    fig, ax = plt.subplots(figsize=(12, 6), dpi=300)

    # Simple line plot with dots
    ax.plot(range(1, n_components+1), variance_ratio * 100,
            marker='o', linewidth=2, markersize=6, color='steelblue')

    ax.set_xlabel('Principal Component', fontsize=13, fontweight='bold')
    ax.set_ylabel('Variance Explained (%)', fontsize=13, fontweight='bold')
    ax.set_title('PCA Scree Plot', fontsize=15, fontweight='bold', pad=15)
    ax.set_xlim(0, n_components + 1)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / '01_pca_scree_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("   ✓ Saved: 01_pca_scree_plot.png")

    # -------------------------------------------------------------------------
    # Figure 2: PCA Loadings Heatmap - Only markers with non-zero loadings
    # -------------------------------------------------------------------------
    print("[2/11] PCA Loadings Heatmap...")

    if 'PCs' in adata.varm:
        # Get loadings for first 5 PCs
        pca_loadings_full = pd.DataFrame(
            adata.varm['PCs'][:, :5],
            index=adata.var_names,
            columns=[f'PC{i+1}' for i in range(5)]
        )

        # Filter to only markers with non-zero loadings (used in PCA)
        marker_variance = pca_loadings_full.abs().sum(axis=1)
        valid_markers = marker_variance[marker_variance > 0].index
        pca_loadings = pca_loadings_full.loc[valid_markers]

        print(f"      Using {len(valid_markers)} markers (out of {len(adata.var_names)})")

        fig, ax = plt.subplots(figsize=(10, 8), dpi=300)

        # Smaller colorbar
        sns.heatmap(pca_loadings, cmap='RdBu_r', center=0, ax=ax,
                    cbar_kws={'label': 'Loading', 'shrink': 0.5},
                    yticklabels=True)

        ax.set_title('PCA Loadings', fontsize=14, fontweight='bold')
        ax.set_xlabel('Principal Component', fontsize=12)
        ax.set_ylabel('Marker', fontsize=12)
        plt.tight_layout()
        plt.savefig(output_dir / '02_pca_loadings_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("   ✓ Saved: 02_pca_loadings_heatmap.png")
    else:
        print("   ⚠ WARNING: PCA loadings not found in adata.varm['PCs']")
        print("   Skipping PCA loadings heatmap...")

    # -------------------------------------------------------------------------
    # Figure 3: PCA Biplot
    # -------------------------------------------------------------------------
    print("[3/11] PCA Biplot...")

    # Filter Unclassified Cells only (no outlier filtering)
    adata_plot = adata[adata.obs['celltype_full'] != 'Unclassified cells'].copy()

    fig, ax = plt.subplots(figsize=(14, 11), dpi=300)
    sc.pl.pca(adata_plot, color='celltype_full', ax=ax, show=False,
              legend_loc='right margin', legend_fontsize=11, s=15,
              title='PCA Biplot', frameon=False)
    plt.tight_layout()
    plt.savefig(output_dir / '03_pca_biplot.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("   ✓ Saved: 03_pca_biplot.png")

    print("\n✅ PCA plots complete (3/11)")


def generate_umap_plots(adata, output_dir):
    """Generate 2 UMAP visualization plots (improved)"""
    print_section("GENERATING UMAP PLOTS (2 figures)")

    if 'X_umap' not in adata.obsm:
        print("⚠️  UMAP data not found, skipping UMAP plots")
        return

    # -------------------------------------------------------------------------
    # Figure 4: UMAP by Cell Type
    # -------------------------------------------------------------------------
    print("\n[4/11] UMAP - Cell Types...")

    # Filter Unclassified Cells
    adata_plot = adata[adata.obs['celltype_full'] != 'Unclassified cells'].copy()

    fig, ax = plt.subplots(figsize=(16, 13), dpi=300)
    sc.pl.umap(adata_plot, color='celltype_full', ax=ax, show=False,
               title='UMAP - Cell Types',
               legend_loc='right margin', legend_fontsize=9,
               s=20,
               frameon=False)
    plt.tight_layout()
    plt.savefig(output_dir / '04_umap_celltype.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("   ✓ Saved: 04_umap_celltype.png")

    # -------------------------------------------------------------------------
    # Figure 4b: UMAP by Cell Type - SCALED version (arcsinh + z-score)
    # -------------------------------------------------------------------------
    print("[4b/12] UMAP - Cell Types (Scaled)...")

    if 'X_umap_scaled' in adata.obsm:
        # Temporarily replace X_umap with X_umap_scaled for plotting
        adata_scaled = adata_plot.copy()
        adata_scaled.obsm['X_umap'] = adata[adata.obs['celltype_full'] != 'Unclassified cells'].obsm['X_umap_scaled']

        fig, ax = plt.subplots(figsize=(16, 13), dpi=300)
        sc.pl.umap(adata_scaled, color='celltype_full', ax=ax, show=False,
                   title='UMAP - Cell Types (Scaled)',
                   legend_loc='right margin', legend_fontsize=9,
                   s=20,
                   frameon=False)
        plt.tight_layout()
        plt.savefig(output_dir / '04_umap_celltype_scaled.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("   ✓ Saved: 04_umap_celltype_scaled.png")
    else:
        print("   ⚠ Scaled UMAP not found, skipping...")

    # -------------------------------------------------------------------------
    # Figure 5: UMAP by Cohort
    # -------------------------------------------------------------------------
    print("[5/12] UMAP - Cohorts...")

    fig, ax = plt.subplots(figsize=(12, 10), dpi=300)
    sc.pl.umap(adata, color='cohort', ax=ax, show=False,
               title='UMAP - Cohort',
               palette={'Basel': 'steelblue', 'Zürich': 'coral'},
               s=20,
               alpha=0.75,
               frameon=False)
    plt.tight_layout()
    plt.savefig(output_dir / '05_umap_cohort.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("   ✓ Saved: 05_umap_cohort.png")

    print("\n✅ UMAP plots complete (6/12)")


def generate_crossbar_plots(adata, output_dir):
    """Generate 2 crossbar plots showing cell type composition"""
    print_section("GENERATING CROSSBAR PLOTS (2 figures)")

    # Filter Unclassified Cells for cleaner visualization
    adata_plot = adata[adata.obs['celltype_full'] != 'Unclassified cells'].copy()

    # -------------------------------------------------------------------------
    # Figure 6: Cell Type Composition by Cohort
    # -------------------------------------------------------------------------
    print("\n[7/12] Crossbar - Cell Types by Cohort...")

    # Calculate composition (normalized to 100%)
    comp_cohort = adata_plot.obs.groupby(['cohort', 'celltype_full']).size().unstack(fill_value=0)
    comp_cohort_pct = comp_cohort.div(comp_cohort.sum(axis=1), axis=0) * 100

    # Sort cell types by overall abundance
    celltype_order = adata_plot.obs['celltype_full'].value_counts().index.tolist()
    comp_cohort_pct = comp_cohort_pct[celltype_order]

    # Get consistent colors
    n_celltypes = len(celltype_order)
    colors = sns.color_palette('tab20', n_colors=20)[:n_celltypes]

    fig, ax = plt.subplots(figsize=(8, 7), dpi=300)
    comp_cohort_pct.plot(kind='bar', stacked=True, ax=ax,
                         color=colors, width=0.6, edgecolor='white', linewidth=1.5)

    ax.set_xlabel('Cohort', fontsize=13, fontweight='bold')
    ax.set_ylabel('Percentage (%)', fontsize=13, fontweight='bold')
    ax.set_title('Cell Type Composition by Cohort', fontsize=15, fontweight='bold', pad=15)
    ax.legend(title='Cell Type', fontsize=7, title_fontsize=9,
              loc='center left', bbox_to_anchor=(1, 0.5), frameon=True)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    ax.grid(axis='y', alpha=0.3)
    ax.set_ylim(0, 100)

    plt.tight_layout()
    plt.savefig(output_dir / '06_celltype_by_cohort_crossbar.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("   ✓ Saved: 06_celltype_by_cohort_crossbar.png")

    # -------------------------------------------------------------------------
    # Figure 7: Cell Type Composition by Patient
    # -------------------------------------------------------------------------
    print("[8/12] Crossbar - Cell Types by Patient...")

    # Get top 20 patients by cell count
    top_patients = adata_plot.obs['Patient'].value_counts().head(20).index.tolist()
    adata_top = adata_plot[adata_plot.obs['Patient'].isin(top_patients)].copy()

    # Calculate composition (normalized to 100%)
    comp_patient = adata_top.obs.groupby(['Patient', 'celltype_full']).size().unstack(fill_value=0)
    comp_patient_pct = comp_patient.div(comp_patient.sum(axis=1), axis=0) * 100

    # Sort patients by total cell count
    patient_order = adata_top.obs['Patient'].value_counts().index.tolist()
    comp_patient_pct = comp_patient_pct.loc[patient_order]

    # Only use cell types that exist in top 20 patients
    celltype_order_subset = [ct for ct in celltype_order if ct in comp_patient_pct.columns]
    comp_patient_pct = comp_patient_pct[celltype_order_subset]

    fig, ax = plt.subplots(figsize=(16, 7), dpi=300)
    comp_patient_pct.plot(kind='bar', stacked=True, ax=ax,
                          color=colors, width=0.85, edgecolor='white', linewidth=0.8)

    ax.set_xlabel('Patient', fontsize=13, fontweight='bold')
    ax.set_ylabel('Percentage (%)', fontsize=13, fontweight='bold')
    ax.set_title('Cell Type Composition by Patient', fontsize=15, fontweight='bold', pad=15)

    # Remove legend (too many cell types)
    ax.get_legend().remove()

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=9)
    ax.grid(axis='y', alpha=0.3)
    ax.set_ylim(0, 100)

    plt.tight_layout()
    plt.savefig(output_dir / '07_celltype_by_patient_crossbar.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("   ✓ Saved: 07_celltype_by_patient_crossbar.png")

    print("\n✅ Crossbar plots complete (8/12)")


def generate_qc_plots(adata, output_dir):
    """Generate 4 QC plots"""
    print_section("GENERATING QC PLOTS (4 figures)")

    # -------------------------------------------------------------------------
    # Figure 8: Cells per Patient 
    # -------------------------------------------------------------------------
    print("\n[9/12] QC - Cells per Patient...")

    patient_counts = adata.obs.groupby('Patient').size().sort_values(ascending=False)

    fig, ax = plt.subplots(figsize=(16, 7), dpi=300)

    # Better colors by cohort
    colors_cohort = ['#4682B4' if not str(p).startswith('Z') else '#FF6347'  # More saturated
                     for p in patient_counts.index]

    bars = ax.bar(range(len(patient_counts)), patient_counts.values,
                   color=colors_cohort, edgecolor='black', linewidth=0.5, alpha=0.85)

    ax.set_xlabel('Patient', fontsize=13, fontweight='bold')
    ax.set_ylabel('Number of Cells', fontsize=13, fontweight='bold')
    ax.set_title('Cell Counts per Patient',
                 fontsize=15, fontweight='bold', pad=15)

    # Add median line
    median_val = patient_counts.median()
    ax.axhline(median_val, color='red', linestyle='--', linewidth=2,
               label=f'Median: {median_val:.0f} cells', zorder=10)

    # Grid for easier reading
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.legend(fontsize=11, loc='upper right', frameon=True, shadow=True)

    # Add cohort legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#4682B4', label='Basel'),
                       Patch(facecolor='#FF6347', label='Zürich')]
    ax2 = ax.twinx()
    ax2.set_yticks([])
    ax2.legend(handles=legend_elements, title='Cohort', loc='upper left',
               fontsize=10, title_fontsize=11, frameon=True, shadow=True)

    plt.tight_layout()
    plt.savefig(output_dir / '08_qc_cells_per_patient.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("   ✓ Saved: 08_qc_cells_per_patient.png")

    # -------------------------------------------------------------------------
    # Figure 9: Cell Area Distribution
    # -------------------------------------------------------------------------
    print("[10/12] QC - Cell Area Distribution...")

    if 'area' in adata.obs.columns:
        area_data = adata.obs['area']

        fig, ax = plt.subplots(figsize=(11, 7), dpi=300)

        # Plot only cells <= 8 μm²
        ax.hist(area_data[area_data <= 8], bins=80,
                color='steelblue', alpha=0.75, edgecolor='black', linewidth=0.5)

        ax.set_xlabel('Cell Area (μm²)', fontsize=13, fontweight='bold')
        ax.set_ylabel('Frequency', fontsize=13, fontweight='bold')
        ax.set_title('Cell Area Distribution', fontsize=15, fontweight='bold', pad=15)
        ax.set_xlim(0, 8)

        # Add median line
        median_val = area_data.median()
        ax.axvline(median_val, color='red', linestyle='--', linewidth=2,
                   label=f'Median: {median_val:.2f} μm²')

        ax.grid(axis='y', alpha=0.3)
        ax.legend(fontsize=11, loc='upper right', frameon=True)

        plt.tight_layout()
        plt.savefig(output_dir / '09_qc_cell_area.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("   ✓ Saved: 09_qc_cell_area.png")

    # -------------------------------------------------------------------------
    # Figure 10: Cohort Comparison 
    # -------------------------------------------------------------------------
    print("[11/12] QC - Cohort Comparison...")

    cohort_summary = adata.obs.groupby('cohort', observed=True).agg({
        'Patient': 'nunique'
    }).rename(columns={'Patient': 'n_patients'})

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=300)

    # Cells per cohort
    cohort_cell_counts = adata.obs['cohort'].value_counts()

    bars1 = ax1.bar(cohort_cell_counts.index, cohort_cell_counts.values,
                     color=['#4682B4', '#FF6347'], edgecolor='black',
                     linewidth=1.5, width=0.6, alpha=0.85)

    ax1.set_ylabel('Number of Cells', fontsize=13, fontweight='bold')
    ax1.set_title('Cells per Cohort', fontsize=13, fontweight='bold', pad=10)
    ax1.grid(axis='y', alpha=0.3, linestyle='--')

    for i, (cohort, count) in enumerate(cohort_cell_counts.items()):
        ax1.text(i, count * 1.02, f'{count:,}\n({100*count/len(adata):.1f}%)',
                 ha='center', va='bottom', fontsize=11, fontweight='bold')

    # Patients per cohort
    bars2 = ax2.bar(cohort_summary.index, cohort_summary['n_patients'],
                     color=['#4682B4', '#FF6347'], edgecolor='black',
                     linewidth=1.5, width=0.6, alpha=0.85)

    ax2.set_ylabel('Number of Patients', fontsize=13, fontweight='bold')
    ax2.set_title('Patients per Cohort', fontsize=13, fontweight='bold', pad=10)
    ax2.grid(axis='y', alpha=0.3, linestyle='--')

    for i, (cohort, n) in enumerate(cohort_summary['n_patients'].items()):
        ax2.text(i, n * 1.02, str(int(n)), ha='center', va='bottom',
                 fontsize=11, fontweight='bold')

    plt.suptitle('Cohort Comparison: Basel vs Zürich',
                 fontsize=16, fontweight='bold', y=1.00)
    plt.tight_layout()
    plt.savefig(output_dir / '10_qc_cohort_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("   ✓ Saved: 10_qc_cohort_comparison.png")

    # -------------------------------------------------------------------------
    # Figure 11: Cell Type Counts 
    # -------------------------------------------------------------------------
    print("[12/12] QC - Cell Type Distribution...")

    # Filter Unclassified Cells
    adata_plot = adata[adata.obs['celltype_full'] != 'Unclassified cells'].copy()
    celltype_counts = adata_plot.obs['celltype_full'].value_counts()

    fig, ax = plt.subplots(figsize=(13, 9), dpi=300)

    # Gradient colors by frequency
    colors_gradient = plt.cm.Blues(np.linspace(0.4, 0.9, len(celltype_counts)))

    bars = ax.barh(range(len(celltype_counts)), celltype_counts.values,
                    color=colors_gradient, edgecolor='black', linewidth=0.8)

    ax.set_yticks(range(len(celltype_counts)))
    ax.set_yticklabels(celltype_counts.index, fontsize=10, fontweight='bold')
    ax.set_xlabel('Number of Cells', fontsize=13, fontweight='bold')
    ax.set_title('Cell Type Distribution',
                 fontsize=15, fontweight='bold', pad=15)

    # Add values on bars
    for i, (count, bar) in enumerate(zip(celltype_counts.values, bars)):
        pct = 100 * count / len(adata_plot)
        ax.text(count * 1.01, i, f'{count:,} ({pct:.1f}%)',
                va='center', fontsize=9, fontweight='bold')

    ax.grid(axis='x', alpha=0.3, linestyle='--')
    ax.invert_yaxis()  # Highest at top

    plt.tight_layout()
    plt.savefig(output_dir / '11_qc_celltype_counts.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("   ✓ Saved: 11_qc_celltype_counts.png")

    print("\n✅ QC plots complete (12/12)")


def main():
    """Main function to generate all embeddings and QC visualizations"""

    print_section("STEP 2: EMBEDDINGS & QC VISUALIZATIONS (IMPROVED)")
    print("\n🎨 Generating 12 high-quality figures for workshop")
    print("   • 3 PCA plots (individual variance, loadings, biplot)")
    print("   • 3 UMAP plots (cell types, cell types scaled, cohort)")
    print("   • 2 Crossbar plots (composition by cohort & patient)")
    print("   • 4 QC plots (cells/patient, area, cohort, cell types)")

    start_time = time.time()

    # Load dataset
    print("\n🔄 Loading dataset...")
    adata_path = Path("data/imc_ERpos_with_embeddings_and_full_names.h5ad")

    if not adata_path.exists():
        print(f"\n❌ ERROR: Dataset not found at {adata_path}")
        print("   Please run Step 0 first to ensure dataset exists")
        return 1

    adata = sc.read_h5ad(adata_path)
    print(f"✅ Loaded: {adata.n_obs:,} cells × {adata.n_vars} markers")

    # Create output directory
    output_dir = Path("figures/step2_embeddings_qc")
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"📁 Output directory: {output_dir}/")

    # Set matplotlib parameters for high quality
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['font.size'] = 10
    plt.rcParams['font.family'] = 'sans-serif'

    # Generate all plots
    generate_pca_plots(adata, output_dir)
    generate_umap_plots(adata, output_dir)
    generate_crossbar_plots(adata, output_dir)
    generate_qc_plots(adata, output_dir)

    # Summary
    elapsed = time.time() - start_time

    print_section("STEP 2 COMPLETE")
    print(f"\n✅ All visualizations generated in {elapsed:.1f} seconds")
    print(f"\n📊 Output Summary:")
    print(f"   • Directory: {output_dir}/")
    print(f"   • Total figures: 12 PNG files (300 DPI)")
    print(f"   • PCA plots: 3")
    print(f"   • UMAP plots: 3 (includes scaled version)")
    print(f"   • Crossbar plots: 2")
    print(f"   • QC plots: 4")

    print(f"\n📝 Data Preprocessing Info:")
    print(f"   • arcsinh transformation applied to raw intensities")
    print(f"   • PCA computed on arcsinh-transformed data")
    print(f"   • 20 PCs used for UMAP (99.9% variance explained)")

    print(f"\n🎯 Next step: Step 3 (Spatial neighbor graphs)")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    exit(main())
