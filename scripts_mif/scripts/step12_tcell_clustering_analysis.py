#!/usr/bin/env python3
"""
STEP 12: T Cell Clustering Analysis

Purpose:
    Analyze whether CD3+ T lymphocytes cluster with each other (homotypic clustering),
    and test if clustering differs between high and low proliferation tumors.

Key Finding (Keren et al. 2018, Eng et al. 2025):
    - High-proliferation tumors have MORE T cell clustering than low-proliferation
    - Quantified by mean number of T cell neighbors per T cell
    - Tukey's HSD P=0.01

Biological Rationale:
    - High-prolif tumors are immunogenic → T cell recruitment
    - T cells cluster in "hot spots" near antigen-rich areas
    - Low-prolif tumors → sparse, isolated T cells
    - Clustering = coordinated immune response

Input:
    - data/imc_ERpos_with_embeddings_and_full_names.h5ad

Outputs:
    - analysis_imc/ERpos_tcell_clustering_by_prolif.csv
    - analysis_imc/ERpos_tcell_clustering_stats.csv
    - figures/step12_tcell_clustering/*.png
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd

warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================
INPUT_H5AD = "data/imc_ERpos_with_embeddings_and_full_names.h5ad"
OUTPUT_DIR = Path("analysis_imc")
FIGURES_DIR = Path("figures/step12_tcell_clustering")

OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
FIGURES_DIR.mkdir(exist_ok=True, parents=True)

# Set matplotlib parameters
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

# Clustering parameters
RADIUS = 25  # μm radius for T-T cell neighbors
PROLIF_QUANTILE = 0.50  # Median split for proliferation

print("=" * 80)
print("STEP 12: T CELL CLUSTERING ANALYSIS")
print("=" * 80)
print(f"Spatial radius: {RADIUS} μm")
print(f"Proliferation cutoff: {PROLIF_QUANTILE} (median)")
print()

# ============================================================================
# STEP 1: LOAD DATA AND CALCULATE PROLIFERATION STATUS
# ============================================================================
print("[1/5] Loading dataset and calculating proliferation fractions...")
t0_total = __import__('time').time()

adata = sc.read_h5ad(INPUT_H5AD)

# Use full cell type names
if 'celltype_full' in adata.obs.columns:
    adata.obs['celltype'] = adata.obs['celltype_full'].astype('category')
else:
    raise ValueError("Column 'celltype_full' not found")

print(f"Loaded: {adata.n_obs:,} cells × {adata.n_vars} markers")
print(f"Patients: {adata.obs['Patient'].nunique()}")

# Calculate proliferation fraction per patient
print("\nCalculating proliferation fractions per patient...")

patient_prolif = []

for patient in adata.obs['Patient'].unique():
    patient_cells = adata.obs[adata.obs['Patient'] == patient]
    n_total = len(patient_cells)

    # Proliferating tumor cells (both Ki67+ types)
    n_prolif = (
        (patient_cells['celltype'] == 'Ki67+ proliferative tumor cells') |
        (patient_cells['celltype'] == 'HER2+ Ki67+ proliferative tumor cells')
    ).sum()
    prolif_frac = n_prolif / n_total

    patient_prolif.append({
        'Patient': patient,
        'prolif_frac': prolif_frac
    })

df_prolif = pd.DataFrame(patient_prolif)

# Classify into high/low proliferation
prolif_cutoff = df_prolif['prolif_frac'].quantile(PROLIF_QUANTILE)
df_prolif['prolif_status'] = df_prolif['prolif_frac'].apply(
    lambda x: 'high' if x >= prolif_cutoff else 'low'
)

print(f"Proliferation cutoff (median): {prolif_cutoff:.6f}")
print(f"Proliferation groups: {df_prolif['prolif_status'].value_counts().to_dict()}")

# ============================================================================
# STEP 2: FILTER FOR CD3+ T LYMPHOCYTES
# ============================================================================
print("\n[2/5] Filtering for CD3+ T lymphocytes...")

tcells = adata[adata.obs['celltype'] == 'CD3+ T lymphocytes'].copy()
print(f"Found {tcells.n_obs:,} CD3+ T lymphocytes")
print(f"Patients with T cells: {tcells.obs['Patient'].nunique()}")

if tcells.n_obs == 0:
    raise ValueError("No CD3+ T lymphocytes found in dataset!")

# ============================================================================
# STEP 3: BUILD SPATIAL GRAPH AND CALCULATE T-T CLUSTERING
# ============================================================================
print(f"\n[3/5] Building spatial graph at radius={RADIUS}μm...")

# Build spatial graph for T cells only
sq.gr.spatial_neighbors(
    tcells,
    radius=RADIUS,
    coord_type="generic",
    spatial_key='spatial'
)

print("Calculating T cell neighbors per T cell...")

# Get adjacency matrix (T cells × T cells)
adj_matrix = tcells.obsp['spatial_connectivities']

# Count neighbors per T cell (sum across rows)
tcell_neighbor_counts = np.array(adj_matrix.sum(axis=1)).flatten()

tcells.obs['tcell_neighbors'] = tcell_neighbor_counts

print(f"T cell neighbors: mean={tcell_neighbor_counts.mean():.2f}, "
      f"median={np.median(tcell_neighbor_counts):.2f}")

# ============================================================================
# STEP 4: AGGREGATE BY PATIENT AND MERGE WITH PROLIFERATION STATUS
# ============================================================================
print("\n[4/5] Aggregating by patient and merging with proliferation status...")

# Calculate mean T-T clustering per patient
df_clustering = tcells.obs.groupby('Patient').agg({
    'tcell_neighbors': ['mean', 'median', 'count']
}).reset_index()

df_clustering.columns = ['Patient', 'mean_tcell_clustering', 'median_tcell_clustering', 'n_tcells']

print(f"Calculated clustering for {len(df_clustering)} patients")

# Merge with proliferation status
df_merged = df_clustering.merge(df_prolif, on='Patient', how='inner')

print(f"Merged dataset: {len(df_merged)} patients")

# ============================================================================
# STEP 5: STATISTICAL COMPARISON BY PROLIFERATION STATUS
# ============================================================================
print("\n[5/5] Comparing T cell clustering between high and low proliferation tumors...")

# Split by proliferation status
high = df_merged[df_merged['prolif_status'] == 'high']['mean_tcell_clustering']
low = df_merged[df_merged['prolif_status'] == 'low']['mean_tcell_clustering']

print(f"\nHigh proliferation (N={len(high)}):")
print(f"  Mean: {high.mean():.3f} ± {high.std():.3f}")
print(f"  Median: {high.median():.3f}")
print(f"  Range: [{high.min():.3f}, {high.max():.3f}]")

print(f"\nLow proliferation (N={len(low)}):")
print(f"  Mean: {low.mean():.3f} ± {low.std():.3f}")
print(f"  Median: {low.median():.3f}")
print(f"  Range: [{low.min():.3f}, {low.max():.3f}]")

# Statistical tests
# 1. Mann-Whitney U test (non-parametric)
statistic_mw, pval_mw = stats.mannwhitneyu(high, low, alternative='two-sided')

# 2. Kruskal-Wallis (equivalent to Mann-Whitney for 2 groups)
statistic_kw, pval_kw = stats.kruskal(high, low)

# 3. Tukey's HSD (as used in paper)
tukey_result = pairwise_tukeyhsd(
    df_merged['mean_tcell_clustering'],
    df_merged['prolif_status'],
    alpha=0.05
)

# 4. Effect size (Cohen's d)
cohens_d = (high.mean() - low.mean()) / np.sqrt(((len(high)-1)*high.var() + (len(low)-1)*low.var()) / (len(high) + len(low) - 2))

print("\n" + "=" * 60)
print("STATISTICAL TESTS")
print("=" * 60)
print(f"Mann-Whitney U test: P={pval_mw:.4f}")
print(f"Kruskal-Wallis test: P={pval_kw:.4f}")
print(f"\nTukey's HSD:")
print(tukey_result)
print(f"\nEffect size (Cohen's d): {cohens_d:.3f}")

if cohens_d > 0.8:
    effect_interpretation = "large"
elif cohens_d > 0.5:
    effect_interpretation = "medium"
elif cohens_d > 0.2:
    effect_interpretation = "small"
else:
    effect_interpretation = "negligible"

print(f"Effect size interpretation: {effect_interpretation}")

# Save results
results_summary = {
    'mann_whitney_p': pval_mw,
    'kruskal_wallis_p': pval_kw,
    'tukey_reject': tukey_result.reject[0],
    'tukey_pval': tukey_result.pvalues[0],
    'cohens_d': cohens_d,
    'high_mean': high.mean(),
    'high_std': high.std(),
    'high_median': high.median(),
    'high_n': len(high),
    'low_mean': low.mean(),
    'low_std': low.std(),
    'low_median': low.median(),
    'low_n': len(low)
}

df_stats = pd.DataFrame([results_summary])
df_stats.to_csv(OUTPUT_DIR / "ERpos_tcell_clustering_stats.csv", index=False)
print(f"\nSaved: {OUTPUT_DIR}/ERpos_tcell_clustering_stats.csv")

df_merged.to_csv(OUTPUT_DIR / "ERpos_tcell_clustering_by_prolif.csv", index=False)
print(f"Saved: {OUTPUT_DIR}/ERpos_tcell_clustering_by_prolif.csv")

# ============================================================================
# VISUALIZATION
# ============================================================================
print("\nGenerating visualizations...")

# Figure 1: Boxplot + Stripplot
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left panel: Boxplot with individual points
ax1 = axes[0]
sns.boxplot(
    data=df_merged,
    x='prolif_status',
    y='mean_tcell_clustering',
    ax=ax1,
    palette={'high': '#d62728', 'low': '#1f77b4'},
    order=['high', 'low']
)
sns.stripplot(
    data=df_merged,
    x='prolif_status',
    y='mean_tcell_clustering',
    ax=ax1,
    color='black',
    alpha=0.5,
    size=4,
    order=['high', 'low']
)

ax1.set_title(f'T Cell Clustering by Proliferation\nMann-Whitney P={pval_mw:.4f}', fontsize=12)
ax1.set_ylabel('Mean T cell neighbors per T cell', fontsize=11)
ax1.set_xlabel('Proliferation status', fontsize=11)
ax1.set_xticklabels(['High', 'Low'])

# Add sample sizes
ax1.text(0, ax1.get_ylim()[1]*0.95, f'n={len(high)}', ha='center', fontsize=9)
ax1.text(1, ax1.get_ylim()[1]*0.95, f'n={len(low)}', ha='center', fontsize=9)

# Right panel: Violin plot
ax2 = axes[1]
sns.violinplot(
    data=df_merged,
    x='prolif_status',
    y='mean_tcell_clustering',
    ax=ax2,
    palette={'high': '#d62728', 'low': '#1f77b4'},
    order=['high', 'low'],
    inner='box'
)

ax2.set_title(f'Distribution of T Cell Clustering\nCohen\'s d={cohens_d:.3f} ({effect_interpretation})', fontsize=12)
ax2.set_ylabel('Mean T cell neighbors per T cell', fontsize=11)
ax2.set_xlabel('Proliferation status', fontsize=11)
ax2.set_xticklabels(['High', 'Low'])

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'tcell_clustering_by_prolif.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: tcell_clustering_by_prolif.png")

# Figure 2: Histogram comparison
fig, ax = plt.subplots(figsize=(10, 6))

ax.hist(high, bins=20, alpha=0.6, color='#d62728', label=f'High prolif (n={len(high)})', edgecolor='black')
ax.hist(low, bins=20, alpha=0.6, color='#1f77b4', label=f'Low prolif (n={len(low)})', edgecolor='black')

ax.axvline(high.mean(), color='#d62728', linestyle='--', linewidth=2, label=f'High mean: {high.mean():.2f}')
ax.axvline(low.mean(), color='#1f77b4', linestyle='--', linewidth=2, label=f'Low mean: {low.mean():.2f}')

ax.set_xlabel('Mean T cell neighbors per T cell', fontsize=11)
ax.set_ylabel('Number of patients', fontsize=11)
ax.set_title(f'Distribution of T Cell Clustering\nMann-Whitney P={pval_mw:.4f}', fontsize=12)
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'tcell_clustering_histogram.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: tcell_clustering_histogram.png")

# Figure 3: Scatter plot (prolif_frac vs tcell_clustering)
fig, ax = plt.subplots(figsize=(10, 6))

# Color by proliferation status
colors = df_merged['prolif_status'].map({'high': '#d62728', 'low': '#1f77b4'})

ax.scatter(
    df_merged['prolif_frac'],
    df_merged['mean_tcell_clustering'],
    c=colors,
    s=50,
    alpha=0.6,
    edgecolors='black',
    linewidth=0.5
)

# Add correlation
from scipy.stats import spearmanr
corr, pval_corr = spearmanr(df_merged['prolif_frac'], df_merged['mean_tcell_clustering'])

ax.set_xlabel('Proliferation fraction', fontsize=11)
ax.set_ylabel('Mean T cell neighbors per T cell', fontsize=11)
ax.set_title(f'T Cell Clustering vs Proliferation\nSpearman r={corr:.3f}, P={pval_corr:.4f}', fontsize=12)
ax.grid(alpha=0.3)

# Add legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#d62728', edgecolor='black', label=f'High prolif (n={len(high)})'),
    Patch(facecolor='#1f77b4', edgecolor='black', label=f'Low prolif (n={len(low)})')
]
ax.legend(handles=legend_elements, fontsize=10)

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'tcell_clustering_vs_proliferation.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: tcell_clustering_vs_proliferation.png")

# ============================================================================
# SUMMARY
# ============================================================================
t_total = __import__('time').time() - t0_total

print("\n" + "=" * 80)
print("STEP 12 COMPLETE")
print("=" * 80)

print("\nKey Findings:")
print(f"  High proliferation: mean T-T clustering = {high.mean():.2f} neighbors")
print(f"  Low proliferation: mean T-T clustering = {low.mean():.2f} neighbors")
print(f"  Difference: {high.mean() - low.mean():.2f} neighbors")

if high.mean() > low.mean():
    print("  → High-prolif tumors have MORE T cell clustering ✓")
else:
    print("  → Low-prolif tumors have MORE T cell clustering (opposite to paper)")

if pval_mw < 0.05:
    print(f"  Statistical significance: P={pval_mw:.4f} (SIGNIFICANT)")
else:
    print(f"  Statistical significance: P={pval_mw:.4f} (not significant)")

print(f"\nBiological Interpretation:")
if high.mean() > low.mean():
    print("  High-prolif tumors are more immunogenic")
    print("  → Active T cell recruitment and clustering")
    print("  → Coordinated immune response")
    print("")
    print("  Low-prolif tumors are less immunogenic")
    print("  → Sparse, isolated T cells")
    print("  → Minimal immune coordination")
else:
    print("  Unexpected pattern - requires further investigation")

print(f"\nComparison to Keren et al. 2018:")
print(f"  Paper: Tukey's HSD P=0.01 (high > low)")
print(f"  Our dataset: Mann-Whitney P={pval_mw:.4f}")

print(f"\nOutputs saved to:")
print(f"  - {OUTPUT_DIR}/ERpos_tcell_clustering_by_prolif.csv")
print(f"  - {OUTPUT_DIR}/ERpos_tcell_clustering_stats.csv")
print(f"  - {FIGURES_DIR}/tcell_clustering_by_prolif.png")
print(f"  - {FIGURES_DIR}/tcell_clustering_histogram.png")
print(f"  - {FIGURES_DIR}/tcell_clustering_vs_proliferation.png")

print(f"\nTotal time: {t_total:.1f}s")
print("=" * 80)
