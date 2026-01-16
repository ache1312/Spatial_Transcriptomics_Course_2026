#!/usr/bin/env python3
"""
STEP 13: Lineage-Level Neighbor Analysis

Purpose:
    Analyze spatial neighbors at the lineage level (epithelial, stromal, immune, endothelial)
    rather than individual cell types. Test whether epithelial-stromal proximity predicts survival.

Key Finding (Keren et al. 2018, Eng et al. 2025):
    - Epithelial cells with MORE stromal neighbors → better RFS
    - Discovery cohort: log-rank P=0.018
    - Validation cohort: FDR=0.028
    - Multivariable Cox: P=0.02 (independent of clinical variables)

Biological Rationale:
    - Stromal fibroblasts provide ECM scaffolding
    - Cancer-associated fibroblasts (CAFs) secrete growth factors
    - Epithelial-stromal crosstalk regulates invasion
    - "Reactive stroma" may constrain tumor spread

Input:
    - data/imc_ERpos_with_embeddings_and_full_names.h5ad

Outputs:
    - analysis_imc/ERpos_lineage_neighbors_40um.csv
    - analysis_imc/ERpos_lineage_neighbors_survival.csv
    - figures/step13_lineage_neighbors/*.png (KM curves)
"""

import sys
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.spatial import cKDTree

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import multivariate_logrank_test

warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================
INPUT_H5AD = "data/imc_ERpos_with_embeddings_and_full_names.h5ad"
OUTPUT_DIR = Path("analysis_imc")
FIGURES_DIR = Path("figures/step13_lineage_neighbors")

OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
FIGURES_DIR.mkdir(exist_ok=True, parents=True)

# Set matplotlib parameters
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

# Spatial parameters
RADIUS = 40  # μm radius for lineage neighbors (standard in Keren et al.)

# Lineage mapping (full cell type names)
LINEAGE_MAP = {
    # Epithelial (all tumor cells)
    'Luminal ER+ tumor cells': 'epithelial',
    'Luminal tumor cells': 'epithelial',
    'HER2+ ER+ tumor cells': 'epithelial',
    'HER2+ ER+ PR+ tumor cells': 'epithelial',
    'Ki67+ proliferative tumor cells': 'epithelial',
    'HER2+ Ki67+ proliferative tumor cells': 'epithelial',
    'Basal tumor cells': 'epithelial',
    'Cytokeratin-low tumor cells': 'epithelial',
    'CD44+ tumor cells': 'epithelial',
    'EGFR+ tumor cells': 'epithelial',
    'Myoepithelial cells': 'epithelial',

    # Immune
    'CD3+ T lymphocytes': 'immune',
    'CD20+ B lymphocytes': 'immune',
    'CD68+ macrophages': 'immune',

    # Stromal
    'Quiescent stroma': 'stromal',
    'Vimentin+ fibroblasts': 'stromal',
    'Fibronectin+ fibroblasts': 'stromal',
    'SMA+ pericytes and fibroblasts': 'stromal',

    # Endothelial
    'CD31+ endothelial cells': 'endothelial',

    # Exclude unclassified
    'Unclassified cells': 'unclassified'
}

print("=" * 80)
print("STEP 13: LINEAGE-LEVEL NEIGHBOR ANALYSIS")
print("=" * 80)
print(f"Spatial radius: {RADIUS} μm")
print(f"Lineages: epithelial, stromal, immune, endothelial")
print()

# ============================================================================
# STEP 1: LOAD DATA AND MAP TO LINEAGES
# ============================================================================
print("[1/5] Loading dataset and mapping cell types to lineages...")
t0_total = __import__('time').time()

adata = sc.read_h5ad(INPUT_H5AD)

# Use full cell type names
if 'celltype_full' in adata.obs.columns:
    adata.obs['celltype'] = adata.obs['celltype_full'].astype('category')
else:
    raise ValueError("Column 'celltype_full' not found")

print(f"Loaded: {adata.n_obs:,} cells × {adata.n_vars} markers")
print(f"Patients: {adata.obs['Patient'].nunique()}")

# Map cell types to lineages
adata.obs['lineage'] = adata.obs['celltype'].map(LINEAGE_MAP)

# Check for unmapped cell types
unmapped = adata.obs[adata.obs['lineage'].isna()]['celltype'].unique()
if len(unmapped) > 0:
    print(f"\nWarning: Unmapped cell types (will be excluded): {unmapped}")

# Exclude unclassified cells
adata = adata[adata.obs['lineage'] != 'unclassified'].copy()

# Count cells per lineage
lineage_counts = adata.obs['lineage'].value_counts()
print(f"\nLineage counts:")
for lineage, count in lineage_counts.items():
    pct = count / adata.n_obs * 100
    print(f"  {lineage}: {count:,} cells ({pct:.1f}%)")

# ============================================================================
# STEP 2: CALCULATE LINEAGE NEIGHBORS (PATIENT-BY-PATIENT)
# ============================================================================
def calculate_lineage_neighbors(adata, source_lineage, target_lineage, radius=40, roi_col='slide_scene'):
    """
    Calculate mean number of target_lineage neighbors for each source_lineage cell.
    Uses KDTree for memory-efficient neighbor search, processing patient-by-patient.
    """
    print(f"\n[Calculating {target_lineage} neighbors of {source_lineage} at radius={radius}μm...]")

    # Get spatial coordinates
    if 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
    else:
        raise ValueError("No 'spatial' coordinates found in adata.obsm")

    if roi_col not in adata.obs.columns:
        raise ValueError(f"Expected ROI column '{roi_col}' not found in adata.obs (needed to avoid cross-ROI neighbors)")

    # Process patient-by-patient to reduce memory usage
    patient_results = []
    patients = adata.obs['Patient'].unique()

    for i, patient in enumerate(patients):
        if (i + 1) % 20 == 0:
            print(f"  Progress: {i+1}/{len(patients)} patients")

        patient_mask = (adata.obs['Patient'] == patient).to_numpy()
        if not patient_mask.any():
            continue

        roi_values = adata.obs.loc[patient_mask, roi_col].unique()

        all_neighbor_counts = []
        total_source = 0
        total_target = 0

        for roi in roi_values:
            roi_mask = patient_mask & (adata.obs[roi_col].to_numpy() == roi)
            if not roi_mask.any():
                continue

            roi_coords = coords[roi_mask]
            roi_lineages = adata.obs.loc[roi_mask, 'lineage'].values

            source_mask = roi_lineages == source_lineage
            target_mask = roi_lineages == target_lineage

            n_source = int(source_mask.sum())
            n_target = int(target_mask.sum())

            if n_source == 0:
                continue

            total_source += n_source
            total_target += n_target

            if n_target == 0:
                all_neighbor_counts.append(np.zeros(n_source, dtype=int))
                continue

            target_coords = roi_coords[target_mask]
            tree = cKDTree(target_coords)

            source_coords = roi_coords[source_mask]
            neighbor_counts = tree.query_ball_point(source_coords, r=radius, return_length=True)
            all_neighbor_counts.append(np.asarray(neighbor_counts, dtype=int))

        if total_source == 0:
            continue

        neighbor_counts_all = np.concatenate(all_neighbor_counts) if all_neighbor_counts else np.zeros(total_source, dtype=int)

        patient_results.append({
            'Patient': patient,
            'mean_neighbors': float(np.mean(neighbor_counts_all)),
            'median_neighbors': float(np.median(neighbor_counts_all)),
            'n_source_cells': int(total_source),
            'n_target_cells': int(total_target)
        })

    df_result = pd.DataFrame(patient_results)

    print(f"  Calculated for {len(df_result)} patients")
    print(f"  Mean neighbors: {df_result['mean_neighbors'].mean():.2f}")
    print(f"  Median neighbors: {df_result['mean_neighbors'].median():.2f}")

    return df_result


print("\n[2/5] Calculating lineage-level neighbors at 40μm radius...")

# Stromal neighbors of epithelial (main finding from paper)
df_stromal_epithelial = calculate_lineage_neighbors(adata, 'epithelial', 'stromal', radius=RADIUS)
df_stromal_epithelial = df_stromal_epithelial.set_index('Patient')
df_stromal_epithelial.columns = [
    'stromal_neighbors_of_epithelial_mean',
    'stromal_neighbors_of_epithelial_median',
    'n_epithelial_cells',
    'n_stromal_cells'
]

# Immune neighbors of epithelial
df_immune_epithelial = calculate_lineage_neighbors(adata, 'epithelial', 'immune', radius=RADIUS)
df_immune_epithelial = df_immune_epithelial.set_index('Patient')
df_immune_epithelial.columns = [
    'immune_neighbors_of_epithelial_mean',
    'immune_neighbors_of_epithelial_median',
    'n_epithelial_cells_immune',
    'n_immune_cells'
]

# Endothelial neighbors of epithelial
df_endothelial_epithelial = calculate_lineage_neighbors(adata, 'epithelial', 'endothelial', radius=RADIUS)
df_endothelial_epithelial = df_endothelial_epithelial.set_index('Patient')
df_endothelial_epithelial.columns = [
    'endothelial_neighbors_of_epithelial_mean',
    'endothelial_neighbors_of_epithelial_median',
    'n_epithelial_cells_endo',
    'n_endothelial_cells'
]

# Merge all
df_features = df_stromal_epithelial.join(df_immune_epithelial, how='outer')
df_features = df_features.join(df_endothelial_epithelial, how='outer')

# Drop duplicate epithelial count columns (keep first)
df_features = df_features.drop(columns=['n_epithelial_cells_immune', 'n_epithelial_cells_endo'])

print(f"\nMerged features for {len(df_features)} patients")

# Save features
df_features.to_csv(OUTPUT_DIR / "ERpos_lineage_neighbors_40um.csv")
print(f"Saved: {OUTPUT_DIR}/ERpos_lineage_neighbors_40um.csv")

# ============================================================================
# STEP 3: EXTRACT CLINICAL DATA FROM H5AD
# ============================================================================
print("\n[3/5] Extracting clinical data from h5ad...")

# Extract clinical data per patient
clinical_data = []

for patient in adata.obs['Patient'].unique():
    patient_cells = adata.obs[adata.obs['Patient'] == patient]

    # Get first row (all clinical vars are same within patient)
    row = patient_cells.iloc[0]

    clinical_data.append({
        'Patient': patient,
        'OS': row.get('Survival_time', np.nan),
        'OSstate': row.get('Survival', np.nan),
        'RFS': row.get('Recurrence_time', np.nan),
        'RFSstate': row.get('Recurrence', np.nan),
        # Optional covariates for adjusted Cox (if available)
        'age': row.get('age', np.nan),
        'grade': row.get('grade', np.nan),
        'tumor_size': row.get('tumor_size', np.nan),
        'Stage': row.get('Stage', np.nan),
        'N_Stage': row.get('N_Stage', np.nan),
        'cohort': row.get('cohort', np.nan),
    })

df_clinical = pd.DataFrame(clinical_data).set_index('Patient')

# Filter to patients with survival data
df_clinical_os = df_clinical[df_clinical['OS'].notna()].copy()
df_clinical_rfs = df_clinical[df_clinical['RFS'].notna()].copy()

print(f"Patients with OS data: {len(df_clinical_os)}")
print(f"  OS events: {int(df_clinical_os['OSstate'].sum())}")
print(f"Patients with RFS data: {len(df_clinical_rfs)}")
print(f"  RFS events: {int(df_clinical_rfs['RFSstate'].sum())}")

# ============================================================================
# STEP 4: SURVIVAL ANALYSIS
# ============================================================================
def survival_analysis(df_features, df_clinical, metric_name, outcome='RFS'):
    """
    Perform survival analysis for a given metric.
    """
    time_col = outcome
    event_col = f'{outcome}state'

    # Merge
    df = df_features[[metric_name]].merge(df_clinical, left_index=True, right_index=True, how='inner')
    df = df.dropna(subset=[metric_name, time_col, event_col])

    if len(df) < 10:
        print(f"Insufficient data for {metric_name} ({len(df)} patients)")
        return {}

    # Split by median
    median_val = df[metric_name].median()
    df['group'] = np.where(df[metric_name] >= median_val, 'high', 'low')

    # Log-rank test
    try:
        lr_result = multivariate_logrank_test(
            df[time_col],
            df['group'],
            df[event_col]
        )
        lr_pval = lr_result.p_value
    except Exception as e:
        print(f"Log-rank failed: {e}")
        lr_pval = np.nan

    # Univariate Cox
    try:
        df_cox = df[[time_col, event_col, metric_name]].copy()
        cph = CoxPHFitter()
        cph.fit(df_cox, duration_col=time_col, event_col=event_col)
        hr = cph.hazard_ratios_[metric_name]
        ci_low, ci_high = cph.confidence_intervals_.loc[metric_name].values
        ci_low, ci_high = float(np.exp(ci_low)), float(np.exp(ci_high))
        p_cox = cph.summary.loc[metric_name, 'p']
    except Exception as e:
        print(f"Cox failed: {e}")
        hr, ci_low, ci_high, p_cox = np.nan, np.nan, np.nan, np.nan

    results = {
        'metric': metric_name,
        'outcome': outcome,
        'n_patients': len(df),
        'n_events': int(df[event_col].sum()),
        'median_value': median_val,
        'logrank_p': lr_pval,
        'cox_hr': hr,
        'cox_ci_low': ci_low,
        'cox_ci_high': ci_high,
        'cox_p': p_cox
    }

    return results, df


def _zscore(series):
    s = pd.to_numeric(series, errors="coerce")
    mu = s.mean()
    sd = s.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return s * np.nan
    return (s - mu) / sd


def fit_multivariable_cox(df_features, df_clinical, metric_name, outcome='RFS', covariates=None):
    """
    Fit adjusted Cox model: metric + clinical covariates.
    Uses z-scored numeric covariates for stability; keeps `metric_name` on original scale.
    """
    time_col = outcome
    event_col = f'{outcome}state'

    if covariates is None:
        covariates = ['age', 'grade', 'tumor_size', 'Stage', 'N_Stage', 'cohort']

    keep_covs = [c for c in covariates if c in df_clinical.columns]
    if not keep_covs:
        return pd.DataFrame()

    df = df_features[[metric_name]].merge(
        df_clinical[[time_col, event_col, *keep_covs]],
        left_index=True,
        right_index=True,
        how='inner',
    )
    df = df.dropna(subset=[metric_name, time_col, event_col])

    model_df = df[[time_col, event_col, metric_name]].copy()
    for c in keep_covs:
        if c == 'cohort':
            continue
        model_df[f'{c}_z'] = _zscore(df[c])

    if 'cohort' in keep_covs:
        cohort = df['cohort'].astype(str)
        dummies = pd.get_dummies(cohort, prefix='cohort', drop_first=True)
        model_df = pd.concat([model_df, dummies], axis=1)

    model_df = model_df.dropna()
    if len(model_df) < 30:
        return pd.DataFrame()

    cph = CoxPHFitter(penalizer=0.1)
    cph.fit(model_df, duration_col=time_col, event_col=event_col)

    summ = cph.summary.reset_index()
    if 'index' in summ.columns:
        summ = summ.rename(columns={'index': 'covariate'})
    elif 'covariate' not in summ.columns:
        summ = summ.rename(columns={summ.columns[0]: 'covariate'})

    summ['HR'] = np.exp(summ['coef'])
    summ['HR_lower_95'] = np.exp(summ['coef lower 95%'])
    summ['HR_upper_95'] = np.exp(summ['coef upper 95%'])
    summ['outcome'] = outcome
    summ['n'] = len(model_df)

    return summ[['outcome', 'covariate', 'HR', 'HR_lower_95', 'HR_upper_95', 'p', 'n']]


def forest_plot_univariate(df_results, outpath):
    """
    Step-10-style forest plot with two panels (OS/RFS) for univariate Cox results.
    """
    if df_results.empty:
        return

    label_map = {
        'stromal_neighbors_of_epithelial_mean': 'Epithelial→Stromal (mean)',
        'immune_neighbors_of_epithelial_mean': 'Epithelial→Immune (mean)',
        'endothelial_neighbors_of_epithelial_mean': 'Epithelial→Endothelial (mean)',
    }
    order = [k for k in label_map.keys() if k in set(df_results['metric'].astype(str))]
    if not order:
        order = list(pd.unique(df_results['metric'].astype(str)))

    fig, axes = plt.subplots(2, 1, figsize=(10.5, 7.0), dpi=250, sharex=True)
    for ax, outcome in zip(axes, ['OS', 'RFS']):
        sub = df_results[df_results['outcome'] == outcome].copy()
        if sub.empty:
            ax.axis('off')
            continue
        sub['metric'] = pd.Categorical(sub['metric'].astype(str), categories=order, ordered=True)
        sub = sub.sort_values('metric')

        hr = sub['cox_hr'].to_numpy(dtype=float)
        lo = sub['cox_ci_low'].to_numpy(dtype=float)
        hi = sub['cox_ci_high'].to_numpy(dtype=float)
        pvals = sub['cox_p'].to_numpy(dtype=float)

        y_pos = np.arange(len(sub))
        colors = ['red' if h > 1 else 'blue' for h in hr]
        ax.scatter(hr, y_pos, s=90, c=colors, zorder=3, alpha=0.8)
        for idx, (l, u) in enumerate(zip(lo, hi)):
            ax.plot([l, u], [idx, idx], c=colors[idx], linewidth=2, alpha=0.6, zorder=2)

        ax.axvline(x=1.0, color='black', linestyle='--', linewidth=1, alpha=0.5, zorder=1)
        ax.set_yticks(y_pos)
        ax.set_yticklabels([label_map.get(m, m) for m in sub['metric'].astype(str)], fontsize=11)
        ax.set_title(outcome, fontsize=12, fontweight='bold', pad=10)
        ax.grid(axis='x', alpha=0.3, zorder=0)

        for idx, (h, u, p) in enumerate(zip(hr, hi, pvals)):
            sig_marker = " *" if np.isfinite(p) and p < 0.05 else ""
            hr_text = f"HR={h:.3f}\np={p:.4f}{sig_marker}"
            x_pos = max(u, ax.get_xlim()[1] * 0.7)
            bg_color = 'lightyellow' if np.isfinite(p) and p < 0.05 else 'white'
            edge_color = 'orange' if np.isfinite(p) and p < 0.05 else 'none'
            ax.text(
                x_pos,
                idx,
                hr_text,
                va='center',
                fontsize=9,
                bbox=dict(boxstyle='round,pad=0.3', facecolor=bg_color, alpha=0.8, edgecolor=edge_color, linewidth=1.5),
            )

    axes[-1].set_xlabel('Hazard Ratio (95% CI)', fontsize=12, fontweight='bold')
    xmax = float(np.nanmax(df_results['cox_ci_high'].to_numpy(dtype=float)))
    xmin = float(np.nanmin(df_results['cox_ci_low'].to_numpy(dtype=float)))
    axes[-1].set_xlim(max(0.05, xmin * 0.9), max(1.2, xmax * 1.1))
    fig.suptitle('Univariate Cox: Lineage Neighbors (40μm)', fontsize=14, fontweight='bold', y=0.98)

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


def forest_plot_multivariable(os_summ, rfs_summ, outpath):
    if (os_summ is None or os_summ.empty) and (rfs_summ is None or rfs_summ.empty):
        return

    def _display_name(cov):
        base = {
            'stromal_neighbors_of_epithelial_mean': 'Epithelial→Stromal (per +1 neighbor)',
            'age_z': 'Age (z)',
            'grade_z': 'Grade (z)',
            'tumor_size_z': 'Tumor size (z)',
            'Stage_z': 'Stage (z)',
            'N_Stage_z': 'N stage (z)',
        }
        if cov in base:
            return base[cov]
        if cov.startswith('cohort_'):
            return f"Cohort: {cov.replace('cohort_', '')}"
        return cov

    fig, axes = plt.subplots(2, 1, figsize=(10.5, 7.0), dpi=250, sharex=True)
    for ax, outcome, summ in zip(axes, ['OS', 'RFS'], [os_summ, rfs_summ]):
        if summ is None or summ.empty:
            ax.axis('off')
            continue
        sub = summ.copy()
        sub['covariate'] = sub['covariate'].astype(str)
        cov_order = []
        if 'stromal_neighbors_of_epithelial_mean' in set(sub['covariate']):
            cov_order.append('stromal_neighbors_of_epithelial_mean')
        cov_order += [c for c in ['age_z', 'grade_z', 'tumor_size_z', 'Stage_z', 'N_Stage_z'] if c in set(sub['covariate'])]
        cov_order += [c for c in sub['covariate'] if c.startswith('cohort_')]
        cov_order += [c for c in sub['covariate'] if c not in set(cov_order)]
        sub['covariate'] = pd.Categorical(sub['covariate'], categories=cov_order, ordered=True)
        sub = sub.sort_values('covariate')

        y_pos = np.arange(len(sub))
        hr = sub['HR'].to_numpy(dtype=float)
        lo = sub['HR_lower_95'].to_numpy(dtype=float)
        hi = sub['HR_upper_95'].to_numpy(dtype=float)
        pvals = sub['p'].to_numpy(dtype=float)
        colors = ['red' if h > 1 else 'blue' for h in hr]

        ax.scatter(hr, y_pos, s=90, c=colors, zorder=3, alpha=0.8)
        for idx, (l, u) in enumerate(zip(lo, hi)):
            ax.plot([l, u], [idx, idx], c=colors[idx], linewidth=2, alpha=0.6, zorder=2)

        ax.axvline(x=1.0, color='black', linestyle='--', linewidth=1, alpha=0.5, zorder=1)
        ax.set_yticks(y_pos)
        ax.set_yticklabels([_display_name(c) for c in sub['covariate'].astype(str)], fontsize=11)
        ax.set_title(outcome, fontsize=12, fontweight='bold', pad=10)
        ax.grid(axis='x', alpha=0.3, zorder=0)

        for idx, (h, u, p) in enumerate(zip(hr, hi, pvals)):
            sig_marker = " *" if np.isfinite(p) and p < 0.05 else ""
            hr_text = f"HR={h:.3f}\np={p:.4f}{sig_marker}"
            x_pos = max(u, ax.get_xlim()[1] * 0.7)
            bg_color = 'lightyellow' if np.isfinite(p) and p < 0.05 else 'white'
            edge_color = 'orange' if np.isfinite(p) and p < 0.05 else 'none'
            ax.text(
                x_pos,
                idx,
                hr_text,
                va='center',
                fontsize=9,
                bbox=dict(boxstyle='round,pad=0.3', facecolor=bg_color, alpha=0.8, edgecolor=edge_color, linewidth=1.5),
            )

    axes[-1].set_xlabel('Hazard Ratio (95% CI)', fontsize=12, fontweight='bold')
    all_hi = []
    all_lo = []
    for s in [os_summ, rfs_summ]:
        if s is None or s.empty:
            continue
        all_hi.append(s['HR_upper_95'].to_numpy(dtype=float))
        all_lo.append(s['HR_lower_95'].to_numpy(dtype=float))
    xmax = float(np.nanmax(np.concatenate(all_hi))) if all_hi else 2.0
    xmin = float(np.nanmin(np.concatenate(all_lo))) if all_lo else 0.5
    axes[-1].set_xlim(max(0.05, xmin * 0.9), max(1.2, xmax * 1.1))
    fig.suptitle('Adjusted Cox: Stromal Neighbors + Clinical Covariates', fontsize=14, fontweight='bold', y=0.98)

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


def plot_survival_curves(df, metric_name, outcome='RFS', figname=None):
    """Plot KM curves"""
    time_col = outcome
    event_col = f'{outcome}state'

    fig, ax = plt.subplots(figsize=(10, 7))

    kmf = KaplanMeierFitter()

    for group in ['high', 'low']:
        df_group = df[df['group'] == group]

        kmf.fit(
            df_group[time_col],
            df_group[event_col],
            label=f'{group.capitalize()} (n={len(df_group)})'
        )
        kmf.plot_survival_function(ax=ax, ci_show=True)

    # Log-rank test
    lr_result = multivariate_logrank_test(
        df[time_col],
        df['group'],
        df[event_col]
    )

    ax.set_xlabel('Time (months)', fontsize=12)
    ax.set_ylabel('Survival Probability', fontsize=12)
    ax.set_title(f'{metric_name}\n{outcome} - Log-rank P={lr_result.p_value:.4f}', fontsize=13)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=11)

    plt.tight_layout()

    if figname is None:
        figname = f'KM_{metric_name}_{outcome}.png'

    plt.savefig(FIGURES_DIR / figname, dpi=300, bbox_inches='tight')
    plt.close()

    return lr_result.p_value


print("\n[4/5] Performing survival analysis...")

# Primary metric for KM curves (kept as the headline finding)
metric = 'stromal_neighbors_of_epithelial_mean'

results_list = []

metrics_to_test = [
    'stromal_neighbors_of_epithelial_mean',
    'immune_neighbors_of_epithelial_mean',
    'endothelial_neighbors_of_epithelial_mean',
]
metrics_to_test = [m for m in metrics_to_test if m in df_features.columns]

res_rfs = None
res_os = None
df_rfs = None
df_os = None

for m in metrics_to_test:
    print(f"\n{m} → RFS:")
    res_rfs_m, df_rfs_m = survival_analysis(df_features, df_clinical_rfs, m, outcome='RFS')
    results_list.append(res_rfs_m)
    if m == metric:
        res_rfs, df_rfs = res_rfs_m, df_rfs_m
    print(f"  N={res_rfs_m['n_patients']}, Events={res_rfs_m['n_events']}")
    print(f"  Median cutoff: {res_rfs_m['median_value']:.2f}")
    print(f"  Log-rank P={res_rfs_m['logrank_p']:.4f}")
    print(f"  Cox HR={res_rfs_m['cox_hr']:.3f} (95% CI: {res_rfs_m['cox_ci_low']:.3f}-{res_rfs_m['cox_ci_high']:.3f}), P={res_rfs_m['cox_p']:.4f}")

    print(f"\n{m} → OS:")
    res_os_m, df_os_m = survival_analysis(df_features, df_clinical_os, m, outcome='OS')
    results_list.append(res_os_m)
    if m == metric:
        res_os, df_os = res_os_m, df_os_m
    print(f"  N={res_os_m['n_patients']}, Events={res_os_m['n_events']}")
    print(f"  Median cutoff: {res_os_m['median_value']:.2f}")
    print(f"  Log-rank P={res_os_m['logrank_p']:.4f}")
    print(f"  Cox HR={res_os_m['cox_hr']:.3f} (95% CI: {res_os_m['cox_ci_low']:.3f}-{res_os_m['cox_ci_high']:.3f}), P={res_os_m['cox_p']:.4f}")

# Save results
df_results = pd.DataFrame(results_list)
df_results.to_csv(OUTPUT_DIR / "ERpos_lineage_neighbors_survival.csv", index=False)
print(f"\nSaved: {OUTPUT_DIR}/ERpos_lineage_neighbors_survival.csv")

# Forest plot (univariate)
forest_plot_univariate(df_results, FIGURES_DIR / "03_cox_forest_lineage_neighbors_univariate.png")
print("  Saved: 03_cox_forest_lineage_neighbors_univariate.png")

# Optional adjusted Cox model (metric + clinical covariates)
try:
    os_adj = fit_multivariable_cox(df_features, df_clinical_os, metric, outcome='OS')
    rfs_adj = fit_multivariable_cox(df_features, df_clinical_rfs, metric, outcome='RFS')
    if not os_adj.empty or not rfs_adj.empty:
        df_adj = pd.concat([os_adj, rfs_adj], ignore_index=True)
        df_adj.to_csv(OUTPUT_DIR / "ERpos_lineage_neighbors_multivariable_cox.csv", index=False)
        print("  Saved: analysis_imc/ERpos_lineage_neighbors_multivariable_cox.csv")

        forest_plot_multivariable(os_adj, rfs_adj, FIGURES_DIR / "04_cox_forest_lineage_neighbors_multivariable.png")
        print("  Saved: 04_cox_forest_lineage_neighbors_multivariable.png")
except Exception as e:  # noqa: BLE001
    print(f"Multivariable Cox failed (skipping): {type(e).__name__}: {e}")

# ============================================================================
# STEP 5: GENERATE KAPLAN-MEIER CURVES
# ============================================================================
print("\n[5/5] Generating Kaplan-Meier curves...")

# RFS
p_rfs = plot_survival_curves(df_rfs, metric, outcome='RFS', figname='KM_stromal_epithelial_RFS.png')
print(f"  Saved: KM_stromal_epithelial_RFS.png (P={p_rfs:.4f})")

# OS
p_os = plot_survival_curves(df_os, metric, outcome='OS', figname='KM_stromal_epithelial_OS.png')
print(f"  Saved: KM_stromal_epithelial_OS.png (P={p_os:.4f})")

# ============================================================================
# SUMMARY
# ============================================================================
t_total = __import__('time').time() - t0_total

print("\n" + "=" * 80)
print("STEP 13 COMPLETE")
print("=" * 80)

print("\nKey Findings:")
print(f"  Stromal neighbors of epithelial (mean): {df_features[metric].mean():.2f}")
print(f"  Range: [{df_features[metric].min():.2f}, {df_features[metric].max():.2f}]")

print(f"\nSurvival associations:")
print(f"  RFS: HR={res_rfs['cox_hr']:.3f}, P={res_rfs['cox_p']:.4f}")
print(f"  OS: HR={res_os['cox_hr']:.3f}, P={res_os['cox_p']:.4f}")

if res_rfs['cox_hr'] < 1:
    print("  → High stromal neighbors = PROTECTIVE (HR<1) ✓")
else:
    print("  → High stromal neighbors = HARMFUL (HR>1)")

print(f"\nComparison to Keren et al. 2018:")
print(f"  Paper: RFS log-rank P=0.018 (significant)")
print(f"  Our dataset: RFS log-rank P={res_rfs['logrank_p']:.4f}")

if res_rfs['logrank_p'] < 0.05:
    print("  ✓ REPLICATED in our ER+ cohort")
else:
    print("  Not significant in our ER+ cohort (may need larger N)")

print(f"\nBiological Interpretation:")
if res_rfs['cox_hr'] < 1:
    print("  Epithelial cells with more stromal neighbors:")
    print("  → Better RFS (slower recurrence)")
    print("  → Stromal cells may constrain tumor invasion")
    print("  → CAFs provide structural scaffolding")
    print("  → 'Desmoplastic' stroma as physical barrier")
else:
    print("  Unexpected pattern - requires further investigation")

print(f"\nOutputs saved to:")
print(f"  - {OUTPUT_DIR}/ERpos_lineage_neighbors_40um.csv")
print(f"  - {OUTPUT_DIR}/ERpos_lineage_neighbors_survival.csv")
print(f"  - {FIGURES_DIR}/KM_stromal_epithelial_RFS.png")
print(f"  - {FIGURES_DIR}/KM_stromal_epithelial_OS.png")

print(f"\nTotal time: {t_total:.1f}s")
print("=" * 80)
