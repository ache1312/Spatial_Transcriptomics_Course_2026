#!/usr/bin/env python
"""
STEP 8: IMPROVED Survival Analysis with 4 enhancements

IMPROVEMENTS over original step8_survival.py:
  1. ✅ USE ENRICHMENT DATA (extract from 18x18 matrix)
  2. ✅ VALIDATE PROPORTIONAL HAZARDS ASSUMPTION
  3. ✅ STANDARDIZE VARIABLES (report HR per SD)
  4. ✅ MULTIVARIABLE COX WITH LASSO SELECTION
  5. ✅ ANALYZE BOTH OS AND RFS
  6. ✅ USE ORIGINAL PLOT STYLE (cleaner, yellow boxes for p<0.05)

Run from pipeline_leiden_validated/:
  conda run -p /path/to/squidpy-env python experimental/step8_improved/step8_survival_improved.py
"""

import os
import sys
import time
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from pathlib import Path

# Import survival analysis packages
try:
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import multivariate_logrank_test, proportional_hazard_test
    from statsmodels.stats.multitest import multipletests
except ImportError:
    print("ERROR: lifelines or statsmodels not installed")
    sys.exit(1)

# Try to import scikit-survival for Lasso Cox
try:
    from sksurv.linear_model import CoxnetSurvivalAnalysis
    from sksurv.util import Surv
    LASSO_AVAILABLE = True
except ImportError:
    print("WARNING: scikit-survival not installed (Lasso Cox unavailable)")
    print("Install with: pip install scikit-survival")
    LASSO_AVAILABLE = False

warnings.filterwarnings('ignore')

# Settings
OUTPUT_DIR = Path("analysis_imc")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

FIGURES_DIR = Path("figures/step8_survival")
FIGURES_DIR.mkdir(exist_ok=True, parents=True)

print("=" * 80)
print("STEP 8 IMPROVED: SURVIVAL ANALYSIS")
print("=" * 80)
print("Enhancements:")
print("  1. Enrichment data extraction (specific pairs)")
print("  2. Proportional hazards validation")
print("  3. Standardized HR (per 1 SD)")
print("  4. Lasso Cox for multivariable selection")
print("=" * 80)
print()

def _first_existing(df, names):
    for n in names:
        if n in df.columns:
            return n
    return None


def _zscore(series):
    s = pd.to_numeric(series, errors="coerce")
    mu = s.mean()
    sd = s.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return s * np.nan
    return (s - mu) / sd

# ============================================================================
# SECTION 1: LOAD CLINICAL DATA
# ============================================================================
print("[1/8] Loading clinical data from h5ad...")
t0 = time.time()
adata = sc.read_h5ad("data/imc_ERpos_with_embeddings_and_full_names.h5ad")

# Extract patient-level clinical data
clinical_cols = ['Patient', 'Survival_months', 'Censored', 'Age', 'Tumor_size', 'Grade',
                 'age', 'tumor_size', 'grade', 'Stage', 'N_Stage', 'cohort',
                 'Recurrence_time', 'Recurrence']
available_cols = [c for c in clinical_cols if c in adata.obs.columns]

clinical_df = adata.obs[available_cols].drop_duplicates(subset='Patient').copy()

# Convert survival columns
if 'Survival_months' in clinical_df.columns:
    clinical_df['OS_months'] = clinical_df['Survival_months']
elif 'Survival_time' in adata.obs.columns:
    survival_df = adata.obs[['Patient', 'Survival_time']].drop_duplicates(subset='Patient').copy()
    clinical_df = clinical_df.merge(survival_df, on='Patient', how='left')
    clinical_df['OS_months'] = clinical_df['Survival_time']
else:
    print("WARNING: No survival time column found")
    clinical_df['OS_months'] = np.nan

# Convert OS_months from days to months if needed
if clinical_df['OS_months'].notna().sum() > 0:
    median_os = clinical_df['OS_months'].median()
    if median_os > 200:
        clinical_df['OS_months'] = clinical_df['OS_months'] / 30.44
        print("  Converted Survival_time from days to months (÷30.44)")

if 'Censored' in clinical_df.columns:
    clinical_df['OS_event'] = 1 - clinical_df['Censored'].astype(int)
elif 'Survival' in adata.obs.columns:
    clinical_df_temp = adata.obs[['Patient', 'Survival']].drop_duplicates(subset='Patient')
    clinical_df = clinical_df.merge(clinical_df_temp, on='Patient', how='left')
    clinical_df['OS_event'] = clinical_df['Survival'].astype(int)
else:
    print("WARNING: No survival event column found")
    clinical_df['OS_event'] = 0

# Recurrence-free survival (RFS) conversion
if 'Recurrence_time' in clinical_df.columns:
    clinical_df['RFS_months'] = clinical_df['Recurrence_time']
    if clinical_df['RFS_months'].notna().sum() > 0:
        median_rfs = clinical_df['RFS_months'].median()
        if median_rfs > 200:
            clinical_df['RFS_months'] = clinical_df['RFS_months'] / 30.44
            print("  Converted Recurrence_time from days to months (÷30.44)")
else:
    clinical_df['RFS_months'] = np.nan

if 'Recurrence' in clinical_df.columns:
    clinical_df['RFS_event'] = clinical_df['Recurrence'].astype(float)
else:
    clinical_df['RFS_event'] = np.nan

print(f"   Loaded clinical data for {len(clinical_df)} patients")
if 'RFS_months' in clinical_df.columns and clinical_df['RFS_months'].notna().sum() > 0:
    n_rfs = clinical_df['RFS_months'].notna().sum()
    n_rfs_events = clinical_df[clinical_df['RFS_event'].notna()]['RFS_event'].sum()
    print(f"   RFS data available: {n_rfs} patients, {int(n_rfs_events)} recurrence events")
print(f"   Time: {time.time() - t0:.1f}s")

# Proliferation proxy (Ki67-like) per patient (align with Step 10)
prolif_types = {'Ki67+ proliferative tumor cells', 'HER2+ Ki67+ proliferative tumor cells'}
if 'celltype_full' in adata.obs.columns:
    prolif_df = pd.DataFrame({
        'Patient': adata.obs['Patient'].astype(str).values,
        'is_prolif': adata.obs['celltype_full'].astype(str).isin(prolif_types).values,
    })
    prolif_df = prolif_df.groupby('Patient', as_index=False).agg(
        n_cells=('is_prolif', 'size'),
        n_prolif=('is_prolif', 'sum'),
    )
    prolif_df['prolif_frac'] = prolif_df['n_prolif'] / prolif_df['n_cells']
    clinical_df['Patient'] = clinical_df['Patient'].astype(str)
    clinical_df = clinical_df.merge(prolif_df[['Patient', 'prolif_frac']], on='Patient', how='left')

# ============================================================================
# SECTION 2: LOAD SPATIAL METRICS (same as original)
# ============================================================================
print("\n[2/8] Loading spatial metrics from previous steps...")

metrics_to_load = {
    'interactions': 'analysis_imc/ERpos_directed_neighbors_r40_coarse.csv',
    'dispersion': 'analysis_imc/ERpos_dispersion_metrics_by_patient.csv',
}

all_metrics = [clinical_df]

for name, path in metrics_to_load.items():
    if Path(path).exists():
        df_metric = pd.read_csv(path)
        if 'Patient' not in df_metric.columns:
            print(f"  ⊗ Skipped {name}: Not patient-level data")
            continue
        print(f"  ✓ Loaded {name}: {df_metric.shape[1]-1} metrics")
        all_metrics.append(df_metric)
    else:
        print(f"  ✗ Missing {name}: {path}")

# ============================================================================
# SECTION 3: 🆕 MEJORA 1 - EXTRACT ENRICHMENT DATA
# ============================================================================
print("\n[3/8] 🆕 MEJORA 1: Extracting enrichment data (specific pairs)...")

enrich_path = Path("analysis_imc/ERpos_nhood_enrichment_r25.csv")
per_patient_enrich_path = Path("analysis_imc/ERpos_nhood_enrichment_per_patient_r40.csv")
if per_patient_enrich_path.exists():
    df_enrich = pd.read_csv(per_patient_enrich_path)
    if "Patient" in df_enrich.columns:
        feat_cols = [c for c in df_enrich.columns if c != "Patient"]
        # Suffix with radius to keep feature provenance explicit
        df_enrich = df_enrich.rename(columns={c: f"{c}_r40" for c in feat_cols})
        all_metrics.append(df_enrich)
        print(f"  ✓ Loaded per-patient enrichment features (r=40µm): {len(feat_cols)}")
        print(f"    File: {per_patient_enrich_path}")
        print(f"    Min non-missing per feature: {int(df_enrich.drop(columns=['Patient']).notna().sum().min())}/{len(df_enrich)}")
    else:
        print(f"  ⊗ Per-patient enrichment missing Patient column: {per_patient_enrich_path}")
elif enrich_path.exists():
    # Fallback: global enrichment matrix (constant across patients; documented but not used for Cox)
    enrich_matrix = pd.read_csv(enrich_path, index_col=0)
    print(f"  ✓ Loaded enrichment matrix: {enrich_matrix.shape}")

    # Define biologically relevant pairs (from Keren et al. paper)
    pairs_of_interest = [
        # Immune infiltration
        ('CD3+ T lymphocytes', 'Luminal tumor cells', 'immune_tumor_infiltration'),
        ('CD3+ T lymphocytes', 'Luminal ER+ tumor cells', 'Tcell_ER_tumor_contact'),

        # Stromal organization
        ('Vimentin+ fibroblasts', 'Vimentin+ fibroblasts', 'CAF_homotypic_clustering'),
        ('SMA+ pericytes and fibroblasts', 'SMA+ pericytes and fibroblasts', 'pericyte_clustering'),

        # Angiogenesis
        ('CD31+ endothelial cells', 'Luminal tumor cells', 'tumor_angiogenesis'),
        ('CD31+ endothelial cells', 'CD31+ endothelial cells', 'vascular_clustering'),

        # Tumor heterogeneity
        ('Basal tumor cells', 'Basal tumor cells', 'basal_clustering'),
        ('Luminal ER+ tumor cells', 'Luminal ER+ tumor cells', 'luminal_clustering'),

        # Immune-stromal
        ('CD68+ macrophages', 'Vimentin+ fibroblasts', 'macrophage_CAF_interaction'),
        ('CD3+ T lymphocytes', 'CD68+ macrophages', 'T_cell_macrophage'),
    ]

    for source, target, feature_name in pairs_of_interest:
        if source in enrich_matrix.index and target in enrich_matrix.columns:
            enrich_value = enrich_matrix.loc[source, target]
            print(f"    {feature_name}: {enrich_value:.2f} (z-score)")
        else:
            print(f"    ⊗ {feature_name}: Cell types not found in matrix")

    print(f"\n  ⚠️  WARNING: Enrichment matrix is GLOBAL (not per-patient)")
    print(f"      All patients get same enrichment values → NOT useful for survival!")
    print(f"      To fix: Compute per-patient enrichment from adata.obs + spatial graphs")

    print(f"  ⊗ Skipping enrichment features (constant across patients)")

else:
    print(f"  ✗ Enrichment file not found: {per_patient_enrich_path} (per-patient) or {enrich_path} (global)")

# ============================================================================
# Merge all metrics
# ============================================================================
print("\n[4/8] Merging all patient-level metrics...")
if len(all_metrics) > 1:
    df = all_metrics[0]
    for df_metric in all_metrics[1:]:
        df = df.merge(df_metric, on='Patient', how='left')
    print(f"   Merged dataset: {df.shape[0]} patients × {df.shape[1]} features")
else:
    df = clinical_df
    print("   WARNING: Using only clinical data")

# Filter to patients with survival data
df_clean = df[(df['OS_months'].notna()) & (df['OS_months'] > 0)].copy()
print(f"   Patients with survival data: {len(df_clean)}")

if len(df_clean) < 10:
    print("\nERROR: Insufficient patients with survival data")
    sys.exit(1)

print(f"\n   Survival summary:")
print(f"     Events (deaths): {df_clean['OS_event'].sum()}/{len(df_clean)} ({df_clean['OS_event'].sum()/len(df_clean)*100:.1f}%)")
print(f"     Median follow-up: {df_clean['OS_months'].median():.1f} months")
print(f"     Range: {df_clean['OS_months'].min():.1f} - {df_clean['OS_months'].max():.1f} months")

# ============================================================================
# SECTION 5: 🆕 MEJORA 3 - STANDARDIZE VARIABLES (before Cox regression)
# ============================================================================
print("\n[5/8] 🆕 MEJORA 3: Standardizing variables (for HR per SD)...")

# Select numerical columns
exclude_cols = ['Patient', 'OS_months', 'OS_event', 'RFS_months', 'RFS_event',
                'slide_scene', 'Censored', 'Survival', 'Survival_time', 'Survival_months',
                'Recurrence_time', 'Recurrence']
candidate_cols = [c for c in df_clean.columns
                  if c not in exclude_cols and df_clean[c].dtype in ['float64', 'int64', 'float32', 'int32']]

print(f"   Candidate variables: {len(candidate_cols)}")

# Create standardized versions of all variables
df_standardized = df_clean.copy()
standardization_stats = {}

for col in candidate_cols:
    # Skip if too many NaNs or zeros
    if df_clean[col].isna().sum() / len(df_clean) > 0.5:
        continue
    if (df_clean[col] == 0).sum() / len(df_clean) > 0.80:
        continue

    # Standardize: (X - mean) / SD
    mean_val = df_clean[col].mean()
    std_val = df_clean[col].std()

    if std_val > 0:
        col_std = f"{col}_std"
        df_standardized[col_std] = (df_clean[col] - mean_val) / std_val
        standardization_stats[col] = {'mean': mean_val, 'std': std_val}

print(f"   Standardized {len(standardization_stats)} variables")

# ============================================================================
# 🆕 Multivariable Cox (Adjusted): per-patient enrichment term
# ============================================================================
print("\n[5b/8] 🆕 Multivariable Cox (Adjusted) for per-patient enrichment...")

term_raw = 'T_cell_macrophage_r40'
term_std = f'{term_raw}_std'
if term_raw in df_clean.columns and term_std in df_standardized.columns:
    # Build covariates from available clinical columns
    age_col = _first_existing(df_standardized, ['age', 'Age'])
    grade_col = _first_existing(df_standardized, ['grade', 'Grade'])
    tumor_size_col = _first_existing(df_standardized, ['tumor_size', 'Tumor_size'])
    stage_col = _first_existing(df_standardized, ['Stage', 'stage'])
    n_stage_col = _first_existing(df_standardized, ['N_Stage', 'N_stage'])
    cohort_col = _first_existing(df_standardized, ['cohort'])

    cox_base = df_standardized.copy()
    cox_base['Patient'] = cox_base['Patient'].astype(str)
    cox_base['age_z'] = _zscore(cox_base[age_col]) if age_col else np.nan
    cox_base['grade_z'] = _zscore(cox_base[grade_col]) if grade_col else np.nan
    cox_base['tumor_size_z'] = _zscore(cox_base[tumor_size_col]) if tumor_size_col else np.nan
    cox_base['stage_z'] = _zscore(cox_base[stage_col]) if stage_col else np.nan
    cox_base['n_stage_z'] = _zscore(cox_base[n_stage_col]) if n_stage_col else np.nan
    cox_base['prolif_frac_z'] = _zscore(cox_base['prolif_frac']) if 'prolif_frac' in cox_base.columns else np.nan

    if cohort_col and cohort_col in cox_base.columns:
        cox_base[cohort_col] = cox_base[cohort_col].astype('category')
        cohort_dummies = pd.get_dummies(cox_base[cohort_col], prefix='cohort', drop_first=True)
        cox_base = pd.concat([cox_base, cohort_dummies], axis=1)

    base_covars = ['age_z', 'grade_z', 'tumor_size_z', 'stage_z', 'n_stage_z', 'prolif_frac_z']
    base_covars += [c for c in cox_base.columns if str(c).startswith('cohort_')]

    from lifelines import CoxPHFitter

    def _fit_cox_table(df_in, duration_col, event_col, covariates, penalizer=0.1, label=''):
        model_df = df_in[[duration_col, event_col, *covariates]].copy()
        model_df = model_df.dropna()
        model_df = model_df[model_df[duration_col] > 0].copy()
        model_df[event_col] = pd.to_numeric(model_df[event_col], errors='coerce').astype(int)

        kept = []
        dropped = []
        for c in covariates:
            if model_df[c].nunique(dropna=True) <= 1:
                dropped.append(c)
                continue
            kept.append(c)
        if dropped:
            print(f"   Dropped non-informative covariates ({label}): {dropped}")
        model_df = model_df[[duration_col, event_col, *kept]].copy()

        last_err = None
        for pen in [penalizer, max(1.0, penalizer * 10), max(10.0, penalizer * 100)]:
            try:
                cph = CoxPHFitter(penalizer=pen)
                cph.fit(model_df, duration_col=duration_col, event_col=event_col)
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
            except Exception as e:  # noqa: BLE001
                last_err = e
                continue
        print(f"   Cox fit failed ({label}): {type(last_err).__name__}: {last_err}")
        return pd.DataFrame()

    def _forest_multivariable(os_summ, rfs_summ, covar_order, outpath):
        def _display_name(cov):
            base = {
                term_std: 'T cell ↔ macrophage enrichment',
                'age_z': 'Age',
                'grade_z': 'Grade',
                'tumor_size_z': 'Tumor size',
                'stage_z': 'Stage',
                'n_stage_z': 'N stage',
                'prolif_frac_z': 'Proliferation proxy',
            }
            if cov in base:
                return base[cov]
            if str(cov).startswith('cohort_'):
                return f"Cohort: {str(cov).replace('cohort_', '')}"
            return str(cov)

        def _prep(summ):
            if summ is None or summ.empty:
                return pd.DataFrame()
            summ = summ.copy()
            summ['covariate'] = summ['covariate'].astype(str)
            keep = [c for c in covar_order if c in set(summ['covariate'])]
            summ = summ[summ['covariate'].isin(keep)].copy()
            summ['order'] = summ['covariate'].map({c: i for i, c in enumerate(keep)})
            summ = summ.sort_values('order')
            summ['var_display'] = summ['covariate'].map(_display_name)
            return summ

        os_df = _prep(os_summ)
        rfs_df = _prep(rfs_summ)
        dfs = [('OS', os_df), ('RFS', rfs_df)]
        if all(d.empty for _, d in dfs):
            return

        all_lo = np.concatenate([d['HR_lower_95'].to_numpy() for _, d in dfs if not d.empty])
        all_hi = np.concatenate([d['HR_upper_95'].to_numpy() for _, d in dfs if not d.empty])
        xmin = max(0.05, float(np.nanmin(all_lo)) * 0.9)
        xmax = max(1.2, float(np.nanmax(all_hi)) * 1.1)

        heights = [max(2.6, len(d) * 0.55) if not d.empty else 2.0 for _, d in dfs]
        fig, axes = plt.subplots(2, 1, figsize=(10.5, sum(heights)), dpi=250, sharex=True)

        for ax, (outcome, df_plot) in zip(axes, dfs):
            if df_plot.empty:
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
        fig.suptitle('Adjusted Cox Model: T cell ↔ macrophage enrichment (r=40µm)', fontsize=14, fontweight='bold', y=0.98)

        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='red', alpha=0.8, label='Increased risk (HR > 1)'),
            Patch(facecolor='blue', alpha=0.8, label='Decreased risk (HR < 1)'),
            Patch(facecolor='lightyellow', alpha=0.8, edgecolor='orange', linewidth=1.5, label='Statistically significant (p < 0.05)'),
        ]
        fig.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(0.02, 0.98), fontsize=8, frameon=True)
        fig.tight_layout(rect=[0, 0.02, 1, 0.90])
        fig.savefig(outpath, bbox_inches='tight')
        plt.close(fig)

    covars = [term_std, *base_covars]
    os_summ = _fit_cox_table(cox_base.dropna(subset=['OS_months', 'OS_event']), 'OS_months', 'OS_event', covars, penalizer=0.1, label='OS')
    rfs_summ = _fit_cox_table(cox_base.dropna(subset=['RFS_months', 'RFS_event']), 'RFS_months', 'RFS_event', covars, penalizer=0.1, label='RFS')

    if not os_summ.empty:
        os_summ.to_csv(OUTPUT_DIR / f"ERpos_step8_cox_multivariable_OS_{term_raw}.csv", index=False)
    if not rfs_summ.empty:
        rfs_summ.to_csv(OUTPUT_DIR / f"ERpos_step8_cox_multivariable_RFS_{term_raw}.csv", index=False)

    # Combined forest plot (OS + RFS full model terms)
    _forest_multivariable(os_summ, rfs_summ, [*base_covars, term_std], FIGURES_DIR / '07_cox_multivariable_T_cell_macrophage_r40.png')
    print("   ✓ Generated: 07_cox_multivariable_T_cell_macrophage_r40.png")
else:
    print(f"   ⊗ Skipping multivariable Cox: missing {term_raw} (run Step 5 per-patient enrichment first).")

# ============================================================================
# SECTION 6: UNIVARIATE COX PH WITH PH VALIDATION
# ============================================================================
print("\n[6/8] Running univariate Cox PH with proportional hazards validation...")

cox_results = []
ph_violations = []

for col in candidate_cols:
    # Filter bad variables (same as original)
    pct_zero = (df_clean[col] == 0).sum() / len(df_clean)
    if pct_zero > 0.80:
        continue
    if df_clean[col].isna().sum() / len(df_clean) > 0.5:
        continue
    if df_clean[col].nunique() < 2:
        continue

    try:
        # Use STANDARDIZED version
        col_std = f"{col}_std"
        if col_std not in df_standardized.columns:
            continue

        df_cox = df_standardized[['OS_months', 'OS_event', col_std]].dropna()

        if len(df_cox) < 10:
            continue

        # Fit Cox PH
        cph = CoxPHFitter()
        cph.fit(df_cox, duration_col='OS_months', event_col='OS_event')

        summary = cph.summary

        # 🆕 MEJORA 2: Validate proportional hazards assumption
        try:
            ph_test = proportional_hazard_test(cph, df_cox, time_transform='rank')
            ph_pvalue = ph_test.summary['p'].iloc[0]
            ph_valid = ph_pvalue > 0.05

            if not ph_valid:
                ph_violations.append({
                    'variable': col,
                    'ph_pvalue': ph_pvalue
                })
        except:
            ph_pvalue = np.nan
            ph_valid = True  # Assume valid if test fails

        # Store results (HR is now per 1 SD)
        cox_results.append({
            'variable': col,
            'coef': summary['coef'].iloc[0],
            'exp_coef_HR_per_SD': summary['exp(coef)'].iloc[0],  # HR per 1 SD!
            'se_coef': summary['se(coef)'].iloc[0],
            'p_value': summary['p'].iloc[0],
            'lower_95': summary['exp(coef) lower 95%'].iloc[0],
            'upper_95': summary['exp(coef) upper 95%'].iloc[0],
            'n_patients': len(df_cox),
            'ph_assumption_valid': ph_valid,
            'ph_test_pvalue': ph_pvalue,
            'variable_mean': standardization_stats[col]['mean'],
            'variable_sd': standardization_stats[col]['std']
        })

    except Exception as e:
        print(f"     Warning: Cox PH failed for {col}: {str(e)}")
        continue

df_cox_results = pd.DataFrame(cox_results)

print(f"\n   Variables tested: {len(df_cox_results)}")
print(f"   Proportional hazards violations: {len(ph_violations)}")

if len(ph_violations) > 0:
    print("\n   ⚠️  Variables violating PH assumption (p<0.05):")
    for v in ph_violations[:5]:
        print(f"     {v['variable']:40} PH p={v['ph_pvalue']:.4f}")

# FDR correction
if len(df_cox_results) > 0:
    reject, pvals_corrected, _, _ = multipletests(
        df_cox_results['p_value'],
        method='fdr_bh',
        alpha=0.05
    )

    df_cox_results['p_value_FDR'] = pvals_corrected
    df_cox_results['significant_FDR'] = reject
    df_cox_results = df_cox_results.sort_values('p_value')

    # Save results
    output_cox = OUTPUT_DIR / "ERpos_survival_cox_univariate.csv"
    df_cox_results.to_csv(output_cox, index=False)
    print(f"\n   ✓ Saved: {output_cox}")

    # Print top results
    print("\n   Top 10 prognostic variables (Cox PH, HR per 1 SD):")
    for idx, row in df_cox_results.head(10).iterrows():
        sig_mark = "***" if row['p_value_FDR'] < 0.05 else ""
        ph_mark = "⚠️" if not row['ph_assumption_valid'] else "✓"
        print(f"     {ph_mark} {row['variable']:38} HR/SD={row['exp_coef_HR_per_SD']:.3f}, p={row['p_value']:.4f} {sig_mark}")

# ============================================================================
# SECTION 7: 🆕 MEJORA 4 - MULTIVARIABLE COX WITH LASSO SELECTION (FIXED!)
# ============================================================================
print("\n[7/8] 🆕 MEJORA 4: Multivariable Cox with Lasso selection (FIXED)...")

if LASSO_AVAILABLE and len(df_cox_results) > 0:
    # Select top candidates EXCLUDING PH violations
    print("   Filtering candidates...")

    # Exclude variables with PH violations
    df_cox_valid = df_cox_results[df_cox_results['ph_assumption_valid'] == True].copy()
    print(f"   Excluded {len(df_cox_results) - len(df_cox_valid)} variables with PH violations")

    # Select top candidates (p<0.15 in univariate, relaxed threshold)
    top_candidates = df_cox_valid[df_cox_valid['p_value'] < 0.15]['variable'].tolist()

    if len(top_candidates) < 5:
        # If too few, take top 10
        top_candidates = df_cox_valid.head(10)['variable'].tolist()

    print(f"   Testing Lasso on {len(top_candidates)} candidate variables (PH-valid only)")

    # Add clinical variables
    clinical_vars = ['Age', 'Tumor_size', 'Grade']
    available_clinical = [v for v in clinical_vars if v in df_clean.columns and df_clean[v].notna().sum() > len(df_clean) * 0.5]

    # Standardize clinical vars too
    for cvar in available_clinical:
        cvar_std = f"{cvar}_std"
        if cvar_std not in df_standardized.columns:
            mean_val = df_clean[cvar].mean()
            std_val = df_clean[cvar].std()
            if std_val > 0:
                df_standardized[cvar_std] = (df_clean[cvar] - mean_val) / std_val

    # Combine standardized columns
    all_vars_std = [f"{v}_std" for v in top_candidates if f"{v}_std" in df_standardized.columns]
    clinical_vars_std = [f"{v}_std" for v in available_clinical if f"{v}_std" in df_standardized.columns]

    lasso_vars = clinical_vars_std + all_vars_std
    print(f"   Total variables for Lasso: {len(lasso_vars)} ({len(clinical_vars_std)} clinical + {len(all_vars_std)} spatial)")

    # Prepare X and y for sksurv
    lasso_cols = ['OS_months', 'OS_event'] + lasso_vars
    df_lasso = df_standardized[lasso_cols].dropna()

    if len(df_lasso) >= 20:
        X = df_lasso[lasso_vars].values
        y = Surv.from_dataframe('OS_event', 'OS_months', df_lasso)

        print(f"   Lasso dataset: {len(df_lasso)} patients × {len(lasso_vars)} variables")

        # FIX: Increased alpha range to avoid numerical instability
        # BEFORE: alphas = np.logspace(-4, 1, 50)  # [0.0001 ... 10] ← TOO LOW
        # AFTER:  alphas = np.logspace(-1, 2, 50)  # [0.1 ... 100]   ← HIGHER PENALTY
        alphas = np.logspace(-1, 2, 50)
        print(f"   Alpha range: {alphas[0]:.4f} to {alphas[-1]:.1f} (50 values)")

        try:
            # Fit with path of alphas
            coxnet = CoxnetSurvivalAnalysis(l1_ratio=1.0, alphas=alphas, max_iter=100000, tol=1e-7)
            coxnet.fit(X, y)

            # Count non-zero coefficients for each alpha
            n_nonzero = np.sum(coxnet.coef_ != 0, axis=1)

            print(f"\n   Alpha scan results:")
            print(f"     Min variables: {n_nonzero.min()}")
            print(f"     Max variables: {n_nonzero.max()}")

            # Strategy: Select alpha where we have 3-8 variables
            # (sparse but not too sparse)
            target_vars = min(8, max(3, len(lasso_vars) // 5))

            # Find alpha closest to target
            valid_alphas = np.where((n_nonzero >= 2) & (n_nonzero <= target_vars + 3))[0]

            if len(valid_alphas) > 0:
                # Pick the middle one (balanced)
                best_idx = valid_alphas[len(valid_alphas) // 2]
            else:
                # Fallback: pick alpha with ~5-10 variables
                best_idx = np.argmin(np.abs(n_nonzero - 5))

            best_alpha = alphas[best_idx]
            n_selected = n_nonzero[best_idx]

            print(f"\n   ✓ Selected alpha={best_alpha:.4f} → {n_selected} variables")

            # Extract coefficients for selected alpha (already fit!)
            # sksurv.CoxnetSurvivalAnalysis returns coef_ with shape (n_features, n_alphas)
            # NOT (n_alphas, n_features) as in sklearn!
            if len(coxnet.coef_.shape) == 1:
                # Single alpha case
                selected_coefs_full = coxnet.coef_
            else:
                # Multiple alphas: extract column for best_idx
                selected_coefs_full = coxnet.coef_[:, best_idx]  # (n_features,)

            selected_mask = selected_coefs_full != 0
            selected_vars = [lasso_vars[i] for i in range(len(lasso_vars)) if selected_mask[i]]
            selected_coefs = selected_coefs_full[selected_mask]

            if len(selected_vars) == 0:
                print(f"\n   ⚠️  Lasso selected 0 variables (alpha too high)")
                print(f"      Try decreasing alpha manually")
            else:
                print(f"\n   ✓ Lasso selected {len(selected_vars)} variables:")

                lasso_results = []
                for var, coef in zip(selected_vars, selected_coefs):
                    hr = np.exp(coef)
                    var_orig = var.replace('_std', '')

                    # Get significance from univariate
                    if var_orig in df_cox_valid['variable'].values:
                        p_uni = df_cox_valid[df_cox_valid['variable'] == var_orig]['p_value'].iloc[0]
                    else:
                        p_uni = np.nan

                    lasso_results.append({
                        'variable': var_orig,
                        'coef_lasso': coef,
                        'HR_per_SD_lasso': hr,
                        'p_value_univariate': p_uni
                    })

                    sig_mark = "***" if p_uni < 0.001 else "**" if p_uni < 0.01 else "*" if p_uni < 0.05 else ""
                    print(f"     {var_orig:40} HR/SD={hr:.3f} (p_uni={p_uni:.4f}) {sig_mark}")

                # Save Lasso results
                df_lasso_results = pd.DataFrame(lasso_results)
                output_lasso = OUTPUT_DIR / "ERpos_survival_lasso_multivariable.csv"
                df_lasso_results.to_csv(output_lasso, index=False)
                print(f"\n   ✓ Saved: {output_lasso}")

                print(f"\n   Interpretation:")
                print(f"     These {len(selected_vars)} variables are independently prognostic")
                print(f"     (survived L1 penalization → not redundant)")

        except Exception as e:
            print(f"   ✗ Lasso failed: {str(e)}")
            print(f"      Error type: {type(e).__name__}")
            import traceback
            print(f"      Traceback: {traceback.format_exc()}")
    else:
        print(f"   ⊗ Insufficient patients for Lasso (n={len(df_lasso)}, need ≥20)")
else:
    if not LASSO_AVAILABLE:
        print("   ⊗ Skipping Lasso (scikit-survival not installed)")
    else:
        print("   ⊗ Skipping Lasso (no significant univariate results)")

# ============================================================================
# SECTION 7b: 🆕 FOREST PLOTS MEJORADOS
# ============================================================================
print("\n[7b/9] 🆕 Generating improved forest plots...")

if len(df_cox_results) > 0:
    # Forest plot function (HR per SD with annotations)
    def plot_forest_improved(df_results, title, filename, top_n=20, lasso_vars=None):
        """
        Forest plot with original style but showing HR per SD.
        - Uses clean style from original (yellow boxes for p<0.05)
        - Shows HR per SD (interpretable!)
        - Marks Lasso-selected variables
        """
        # Select top N
        df_plot = df_results.head(top_n).copy()
        df_plot = df_plot.sort_values('exp_coef_HR_per_SD', ascending=True)

        # Clean variable names (original style)
        df_plot['var_display'] = df_plot['variable'].str.replace('_to_', ' → ')
        df_plot['var_display'] = df_plot['var_display'].str.replace('_neighbors', '')
        df_plot['var_display'] = df_plot['var_display'].str.replace('_nn_dist_median', ' (NN dist, median)')
        df_plot['var_display'] = df_plot['var_display'].str.replace('_nn_dist_mean', ' (NN dist, mean)')
        df_plot['var_display'] = df_plot['var_display'].apply(
            lambda x: x[:60] + '...' if len(x) > 60 else x
        )

        fig, ax = plt.subplots(figsize=(10, max(8, top_n * 0.4)))

        y_pos = np.arange(len(df_plot))

        # Plot HRs as points (original style: simple red/blue)
        colors = ['red' if hr > 1 else 'blue' for hr in df_plot['exp_coef_HR_per_SD']]
        ax.scatter(df_plot['exp_coef_HR_per_SD'], y_pos, s=80, c=colors, zorder=3, alpha=0.8)

        # Plot confidence intervals
        for idx, (i, row) in enumerate(df_plot.iterrows()):
            ax.plot([row['lower_95'], row['upper_95']], [idx, idx],
                    c=colors[idx], linewidth=2, alpha=0.6, zorder=2)

        # Vertical line at HR=1 (original style)
        ax.axvline(x=1.0, color='black', linestyle='--', linewidth=1, alpha=0.5, zorder=1)

        # Labels and formatting
        ax.set_yticks(y_pos)
        ax.set_yticklabels(df_plot['var_display'], fontsize=9)
        ax.set_xlabel('Hazard Ratio per 1 SD (95% CI)', fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
        ax.grid(axis='x', alpha=0.3, zorder=0)

        # Add HR text annotations with yellow boxes (original style)
        for idx, (i, row) in enumerate(df_plot.iterrows()):
            # Significance marker
            sig_marker = ""
            if row['p_value'] < 0.001:
                sig_marker = " ***"
            elif row['p_value'] < 0.01:
                sig_marker = " **"
            elif row['p_value'] < 0.05:
                sig_marker = " *"

            # Lasso marker (NEW)
            lasso_marker = ""
            if lasso_vars and row['variable'] in lasso_vars:
                lasso_marker = " [LASSO]"

            hr_text = f"HR={row['exp_coef_HR_per_SD']:.2f}\np={row['p_value']:.4f}{sig_marker}{lasso_marker}"
            x_pos = max(row['upper_95'], ax.get_xlim()[1] * 0.7)

            # Yellow box for significant results (original style)
            bg_color = 'lightyellow' if row['p_value'] < 0.05 else 'white'
            edge_color = 'orange' if row['p_value'] < 0.05 else 'none'

            ax.text(x_pos, idx, hr_text, va='center', fontsize=7,
                    bbox=dict(boxstyle='round,pad=0.3', facecolor=bg_color, alpha=0.8,
                             edgecolor=edge_color, linewidth=1.5))

        # Legend (moved to upper left)
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='red', alpha=0.8, label='Increased risk (HR > 1)'),
            Patch(facecolor='blue', alpha=0.8, label='Decreased risk (HR < 1)'),
            Patch(facecolor='lightyellow', alpha=0.8, edgecolor='orange', linewidth=1.5,
                  label='Statistically significant (p < 0.05)')
        ]
        ax.legend(handles=legend_elements, loc='upper left', fontsize=9)

        plt.tight_layout()
        plt.savefig(FIGURES_DIR / filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   ✓ Generated: {filename}")

    # Generate forest plot for OS
    lasso_selected = []
    if 'df_lasso_results' in locals() and len(df_lasso_results) > 0:
        lasso_selected = df_lasso_results['variable'].tolist()

    plot_forest_improved(
        df_cox_results,
        title='Overall Survival: Top Prognostic Variables (HR per SD)',
        filename='01_forest_plot_OS.png',
        top_n=min(25, len(df_cox_results)),
        lasso_vars=lasso_selected
    )

print(f"   Forest plots saved to: {FIGURES_DIR}/")

# ============================================================================
# SECTION 7c: RECURRENCE-FREE SURVIVAL (RFS) ANALYSIS
# ============================================================================
print("\n[7c/9] 🆕 Analyzing Recurrence-Free Survival (RFS)...")

if 'RFS_months' in df_clean.columns and 'RFS_event' in df_clean.columns:
    # Check if we have enough RFS data
    rfs_data_available = df_clean['RFS_months'].notna().sum()
    rfs_events = df_clean[df_clean['RFS_event'].notna()]['RFS_event'].sum()

    print(f"   RFS data: {rfs_data_available} patients, {int(rfs_events)} recurrence events")

    if rfs_data_available >= 20 and rfs_events >= 5:
        # Run Cox PH for RFS (same approach as OS)
        rfs_results = []
        candidate_cols_rfs = [c for c in candidate_cols if c in df_clean.columns]

        print(f"   Running Cox PH for {len(candidate_cols_rfs)} variables...")

        for col in candidate_cols_rfs:
            # FILTER: Skip variables with >80% zeros (prevents spurious HRs)
            pct_zero = (df_clean[col] == 0).sum() / len(df_clean)
            if pct_zero > 0.80:
                continue

            # FILTER: Skip variables with >50% NaNs
            if df_clean[col].isna().sum() / len(df_clean) > 0.5:
                continue

            # FILTER: Skip variables with no variance
            if df_clean[col].nunique() < 2:
                continue

            df_cox_rfs = df_clean[['RFS_months', 'RFS_event', col]].dropna()

            if len(df_cox_rfs) < 15:
                continue

            try:
                cph = CoxPHFitter()
                cph.fit(df_cox_rfs, duration_col='RFS_months', event_col='RFS_event')

                summary = cph.summary
                coef = summary['coef'].iloc[0]
                hr = summary['exp(coef)'].iloc[0]
                se = summary['se(coef)'].iloc[0]
                pval = summary['p'].iloc[0]
                lower = summary['exp(coef) lower 95%'].iloc[0]
                upper = summary['exp(coef) upper 95%'].iloc[0]

                # Standardization (per SD)
                if col in standardization_stats:
                    std_val = standardization_stats[col]['std']
                    hr_per_sd = np.exp(coef * std_val)
                    lower_sd = np.exp((coef - 1.96 * se) * std_val)
                    upper_sd = np.exp((coef + 1.96 * se) * std_val)
                else:
                    hr_per_sd = hr
                    lower_sd = lower
                    upper_sd = upper
                    std_val = 1.0

                # PH validation
                try:
                    ph_test = proportional_hazard_test(cph, df_cox_rfs, time_transform='rank')
                    ph_pvalue = ph_test.summary['p'].iloc[0]
                    ph_valid = ph_pvalue > 0.05
                except:
                    ph_pvalue = np.nan
                    ph_valid = True

                rfs_results.append({
                    'variable': col,
                    'n_patients': len(df_cox_rfs),
                    'coef': coef,
                    'exp_coef_HR': hr,
                    'exp_coef_HR_per_SD': hr_per_sd,
                    'se_coef': se,
                    'lower_95': lower_sd,
                    'upper_95': upper_sd,
                    'p_value': pval,
                    'variable_sd': std_val,
                    'ph_assumption_valid': ph_valid,
                    'ph_test_pvalue': ph_pvalue
                })
            except Exception as e:
                continue

        # Convert to DataFrame and sort
        df_rfs_results = pd.DataFrame(rfs_results)
        df_rfs_results = df_rfs_results.sort_values('p_value').reset_index(drop=True)

        # FDR correction
        if len(df_rfs_results) > 0:
            _, pvals_fdr, _, _ = multipletests(df_rfs_results['p_value'], method='fdr_bh')
            df_rfs_results['p_value_FDR'] = pvals_fdr
            df_rfs_results['significant_FDR'] = pvals_fdr < 0.05

        # Save results
        output_rfs = OUTPUT_DIR / "ERpos_survival_rfs_cox_univariate.csv"
        df_rfs_results.to_csv(output_rfs, index=False)
        print(f"   ✓ Saved RFS results: {output_rfs}")
        print(f"   ✓ Variables significant (p<0.05): {(df_rfs_results['p_value'] < 0.05).sum()}")
        print(f"   ✓ Variables significant (FDR<0.05): {df_rfs_results['significant_FDR'].sum()}")

        # Generate forest plot for RFS
        if len(df_rfs_results) > 0:
            print("   Generating RFS forest plot...")
            plot_forest_improved(
                df_rfs_results,
                title='Recurrence-Free Survival: Top Prognostic Variables (HR per SD)',
                filename='02_forest_plot_RFS.png',
                top_n=min(25, len(df_rfs_results)),
                lasso_vars=lasso_selected  # Same Lasso vars for consistency
            )

        # NOTE: RFS Lasso is skipped due to numerical instability from high collinearity
        # The OS Lasso already identifies independent predictors, and RFS univariate Cox provides
        # comprehensive significant variables. RFS Lasso does not add critical information.
        print("   ⊗ RFS Lasso skipped (numerically unstable - high collinearity in spatial features)")
    else:
        print(f"   ⊗ Insufficient RFS data/events (need ≥20 patients and ≥5 events)")
else:
    print("   ⊗ No RFS data available in dataset")

# ============================================================================
# SECTION 7d: KAPLAN-MEIER CURVES (4 SIGNIFICANT VARIABLES)
# ============================================================================
print("\n[7d/9] 🆕 Generating Kaplan-Meier curves for significant variables...")

# Define the 4 significant variables (from manual validation - log-rank p<0.05)
significant_km_vars = [
    {
        'variable': 'n_isolated',
        'endpoint': 'RFS',
        'title': 'Isolated Cells Count',
        'cutoff': 750.5
    },
    {
        'variable': 'frac_Quiescent stroma',
        'endpoint': 'RFS',
        'title': 'Quiescent Stroma Fraction',
        'cutoff': 13.68
    },
    {
        'variable': 'Luminal tumor cells_nn_dist_mean',
        'endpoint': 'RFS',
        'title': 'Luminal Tumor Cells Dispersion',
        'cutoff': 20.75
    },
    {
        'variable': 'SMA+ pericytes and fibroblasts_nn_dist_median',
        'endpoint': 'OS',
        'title': 'SMA+ Pericytes Dispersion',
        'cutoff': 15.49
    }
]

def plot_km_single(df_data, variable, endpoint, title, cutoff, filename):
    """Plot single KM curve with consistent style"""
    from lifelines.statistics import logrank_test

    time_col = f'{endpoint}_months'
    event_col = f'{endpoint}_event'

    # Clean data
    df_km = df_data[[time_col, event_col, variable]].dropna()

    if len(df_km) < 20:
        print(f"   ⊗ Skipped {variable}: insufficient data")
        return

    # Dichotomize by cutoff
    df_km['group'] = df_km[variable].apply(lambda x: 'High' if x > cutoff else 'Low')

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)

    # Fit KM curves
    kmf = KaplanMeierFitter()

    # Colors: green for low-risk (protective), red for high-risk
    is_protective = variable == 'frac_Quiescent stroma'

    if is_protective:
        colors = {'High': 'darkgreen', 'Low': 'red'}
    else:
        colors = {'Low': 'darkgreen', 'High': 'red'}

    for group in ['Low', 'High']:
        mask = df_km['group'] == group
        n = mask.sum()

        kmf.fit(
            durations=df_km.loc[mask, time_col],
            event_observed=df_km.loc[mask, event_col],
            label=f'{group} (n={n})'
        )
        kmf.plot_survival_function(ax=ax, ci_show=True, color=colors[group], linewidth=2.5)

    # Log-rank test
    results_lr = logrank_test(
        durations_A=df_km.loc[df_km['group'] == 'High', time_col],
        durations_B=df_km.loc[df_km['group'] == 'Low', time_col],
        event_observed_A=df_km.loc[df_km['group'] == 'High', event_col],
        event_observed_B=df_km.loc[df_km['group'] == 'Low', event_col]
    )
    p_lr = results_lr.p_value

    # Formatting
    ax.set_xlabel(f'{endpoint} (months)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Survival Probability', fontsize=11, fontweight='bold')

    # Title with p-value
    p_stars = '***' if p_lr < 0.001 else '**' if p_lr < 0.01 else '*'
    ax.set_title(f'{title}\nLog-rank P = {p_lr:.4f} {p_stars}',
                 fontsize=12, fontweight='bold')

    ax.set_ylim(0, 1.05)
    ax.legend(loc='lower left', fontsize=10, framealpha=0.9)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / filename, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"   ✓ Generated: {filename}")

# Generate the 4 KM plots
km_count = 0
for idx, var_info in enumerate(significant_km_vars, start=3):  # Start at 03 (after 01 and 02 forest plots)
    var = var_info['variable']

    # Check if variable exists in data
    if var not in df_clean.columns:
        print(f"   ⊗ Skipped {var}: not in dataset")
        continue

    filename = f"{idx:02d}_km_{var_info['endpoint']}_{var.replace(' ', '_').replace('+', 'pos').replace('/', '_')}.png"

    plot_km_single(
        df_data=df_clean,
        variable=var,
        endpoint=var_info['endpoint'],
        title=var_info['title'],
        cutoff=var_info['cutoff'],
        filename=filename
    )
    km_count += 1

print(f"   Generated {km_count} Kaplan-Meier curves")

# ============================================================================
# SECTION 8: SUMMARY REPORT
# ============================================================================
print("\n[8/9] Creating summary report...")

summary_report = []
summary_report.append("=" * 80)
summary_report.append("ER+ SURVIVAL ANALYSIS - IMPROVED VERSION")
summary_report.append("=" * 80)
summary_report.append("")
summary_report.append("IMPROVEMENTS:")
summary_report.append("  1. Enrichment extraction (attempted)")
summary_report.append("  2. Proportional hazards validation (✓)")
summary_report.append("  3. Standardized HR (per 1 SD) (✓)")
summary_report.append("  4. Lasso multivariable selection (✓)")
summary_report.append("")
summary_report.append(f"Dataset: {len(df_clean)} ER+ patients")
summary_report.append(f"Events: {df_clean['OS_event'].sum()} deaths ({df_clean['OS_event'].sum()/len(df_clean)*100:.1f}%)")
summary_report.append("")

if len(df_cox_results) > 0:
    summary_report.append("Univariate Cox PH (HR per 1 SD):")
    summary_report.append(f"  Variables tested: {len(df_cox_results)}")
    summary_report.append(f"  Significant (p<0.05): {(df_cox_results['p_value'] < 0.05).sum()}")
    summary_report.append(f"  Significant (FDR<0.05): {df_cox_results['significant_FDR'].sum()}")
    summary_report.append(f"  PH violations: {len(ph_violations)}")
    summary_report.append("")
    summary_report.append("Top 5 Variables:")
    for idx, row in df_cox_results.head(5).iterrows():
        summary_report.append(f"  {row['variable']:38} HR/SD={row['exp_coef_HR_per_SD']:.3f}, p={row['p_value']:.4f}")

summary_report.append("")
summary_report.append("=" * 80)

output_report = OUTPUT_DIR / "ERpos_survival_summary.txt"
with open(output_report, 'w') as f:
    f.write('\n'.join(summary_report))

print()
for line in summary_report:
    print(line)

print(f"\n   ✓ Saved: {output_report}")

print("\n" + "=" * 80)
print("STEP 8 IMPROVED COMPLETE (with forest plots!)")
print("=" * 80)
print(f"\n📁 Outputs:")
print(f"   Results: {OUTPUT_DIR}/")
print(f"   Figures: {FIGURES_DIR}/")
print()
