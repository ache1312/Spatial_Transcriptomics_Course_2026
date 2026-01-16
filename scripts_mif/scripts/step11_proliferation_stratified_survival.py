#!/usr/bin/env python3
"""
STEP 11: Proliferation-Stratified Survival Analysis

Purpose:
    Test for interaction between proliferation status and CD3+ T cell infiltration
    in predicting overall survival and recurrence-free survival.

Key Finding from Eng et al. 2025 (Figure 3E-G):
    - CD3+ T cells are prognostic ONLY in high-proliferation ER+ tumors
    - High-prolif + High CD3 → Best prognosis (HR~0.16, P=0.028)
    - Low-prolif tumors → CD3 NOT prognostic (P~0.99)
    - 4-group analysis → P=0.0071

Biological Rationale:
    - High-proliferation tumors are more immunogenic (neoantigen presentation)
    - Immune surveillance is more effective in rapidly dividing tumors
    - Low-proliferation (hormone-driven) tumors rely less on immune evasion

Input:
    - data/imc_ERpos_with_embeddings_and_full_names.h5ad

Outputs:
    - analysis_imc/ERpos_proliferation_stratified_survival.csv
    - analysis_imc/ERpos_4group_survival.csv
    - figures/step11_proliferation_stratified/*.png (KM curves)
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import multivariate_logrank_test

warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================
INPUT_H5AD = "data/imc_ERpos_with_embeddings_and_full_names.h5ad"
OUTPUT_DIR = Path("analysis_imc")
FIGURES_DIR = Path("figures/step11_proliferation_stratified")

OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
FIGURES_DIR.mkdir(exist_ok=True, parents=True)

# Set matplotlib parameters
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

# Stratification parameters
PROLIF_QUANTILE = 0.50  # Median split for proliferation
CD3_QUANTILE = 0.33     # Lower tertile for CD3 (following paper)

print("=" * 80)
print("STEP 11: PROLIFERATION-STRATIFIED SURVIVAL ANALYSIS")
print("=" * 80)
print(f"Proliferation cutoff: {PROLIF_QUANTILE} (median)")
print(f"CD3 cutoff: {CD3_QUANTILE} (lower tertile)")
print()

# ============================================================================
# STEP 1: LOAD DATA AND CALCULATE CELL TYPE FRACTIONS
# ============================================================================
print("[1/5] Loading dataset and computing cell type fractions...")
t0_total = __import__('time').time()

adata = sc.read_h5ad(INPUT_H5AD)

# Use full cell type names
if 'celltype_full' in adata.obs.columns:
    adata.obs['celltype'] = adata.obs['celltype_full'].astype('category')
else:
    raise ValueError("Column 'celltype_full' not found")

print(f"Loaded: {adata.n_obs:,} cells × {adata.n_vars} markers")
print(f"Patients: {adata.obs['Patient'].nunique()}")

# Calculate cell type fractions per patient
print("\nCalculating cell type fractions per patient...")

patient_fractions = []

for patient in adata.obs['Patient'].unique():
    patient_mask = adata.obs['Patient'] == patient
    patient_cells = adata.obs[patient_mask]

    n_total = len(patient_cells)

    # Proliferating tumor cells (both Ki67+ types)
    n_prolif = ((patient_cells['celltype'] == 'Ki67+ proliferative tumor cells') |
                (patient_cells['celltype'] == 'HER2+ Ki67+ proliferative tumor cells')).sum()
    prolif_frac = n_prolif / n_total

    # CD3+ T lymphocytes
    n_cd3 = (patient_cells['celltype'] == 'CD3+ T lymphocytes').sum()
    cd3_frac = n_cd3 / n_total

    # CD20+ B lymphocytes
    n_cd20 = (patient_cells['celltype'] == 'CD20+ B lymphocytes').sum()
    cd20_frac = n_cd20 / n_total

    # CD68+ macrophages
    n_cd68 = (patient_cells['celltype'] == 'CD68+ macrophages').sum()
    cd68_frac = n_cd68 / n_total

    # All tumor cells
    tumor_types = [
        'Basal tumor cells', 'Luminal tumor cells', 'Luminal ER+ tumor cells',
        'Cytokeratin-low tumor cells', 'HER2+ ER+ tumor cells',
        'Ki67+ proliferative tumor cells', 'HER2+ Ki67+ proliferative tumor cells',
        'HER2+ ER+ PR+ tumor cells', 'Myoepithelial cells',
        'CD44+ tumor cells', 'EGFR+ tumor cells'
    ]
    n_tumor = patient_cells['celltype'].isin(tumor_types).sum()
    tumor_frac = n_tumor / n_total

    patient_fractions.append({
        'Patient': patient,
        'n_cells': n_total,
        'prolif_frac': prolif_frac,
        'cd3_frac': cd3_frac,
        'cd20_frac': cd20_frac,
        'cd68_frac': cd68_frac,
        'tumor_frac': tumor_frac
    })

df_fractions = pd.DataFrame(patient_fractions)

print(f"Computed fractions for {len(df_fractions)} patients")
print(f"\nProliferation fraction: mean={df_fractions['prolif_frac'].mean():.4f}, median={df_fractions['prolif_frac'].median():.4f}")
print(f"CD3+ T cell fraction: mean={df_fractions['cd3_frac'].mean():.4f}, median={df_fractions['cd3_frac'].median():.4f}")

# ============================================================================
# STEP 2: MERGE WITH CLINICAL DATA
# ============================================================================
print("\n[2/5] Merging with clinical data...")

# Extract clinical data from h5ad
clinical_cols = ['Patient', 'Survival_time', 'Survival', 'Recurrence_time', 'Recurrence']
available_cols = [c for c in clinical_cols if c in adata.obs.columns]

clinical_df = adata.obs[available_cols].drop_duplicates(subset='Patient').copy()

# Convert survival time to months
if 'Survival_time' in clinical_df.columns:
    if clinical_df['Survival_time'].max() > 500:  # Likely in days
        clinical_df['OS_months'] = clinical_df['Survival_time'] / 30.44
    else:
        clinical_df['OS_months'] = clinical_df['Survival_time']
    clinical_df['OS_event'] = clinical_df['Survival'].astype(int)
else:
    raise ValueError("Survival_time column not found in h5ad")

# Convert recurrence time to months
if 'Recurrence_time' in clinical_df.columns:
    if clinical_df['Recurrence_time'].max() > 500:
        clinical_df['RFS_months'] = clinical_df['Recurrence_time'] / 30.44
    else:
        clinical_df['RFS_months'] = clinical_df['Recurrence_time']
    clinical_df['RFS_event'] = clinical_df['Recurrence'].astype(int)

# Merge fractions with clinical data
df = df_fractions.merge(clinical_df, on='Patient', how='inner')
df = df[df['OS_months'].notna() & (df['OS_months'] > 0)].copy()

print(f"Patients with OS data: {len(df)}")
print(f"OS events: {df['OS_event'].sum()}")

if 'RFS_months' in df.columns:
    df_rfs = df[df['RFS_months'].notna() & (df['RFS_months'] > 0)].copy()
    print(f"Patients with RFS data: {len(df_rfs)}")
    print(f"RFS events: {df_rfs['RFS_event'].sum()}")

# ============================================================================
# STEP 3: CLASSIFY PATIENTS INTO HIGH/LOW GROUPS
# ============================================================================
print("\n[3/5] Classifying patients into high/low groups...")

# Proliferation: median split
prolif_cutoff = df['prolif_frac'].quantile(PROLIF_QUANTILE)
df['prolif_group'] = np.where(df['prolif_frac'] > prolif_cutoff, 'high', 'low')

print(f"\nProliferation cutoff ({PROLIF_QUANTILE} quantile): {prolif_cutoff:.6f}")
print(f"Proliferation groups: {df['prolif_group'].value_counts().to_dict()}")

# CD3: lower tertile (following paper methodology)
cd3_cutoff = df['cd3_frac'].quantile(CD3_QUANTILE)
df['cd3_group'] = np.where(df['cd3_frac'] > cd3_cutoff, 'high', 'low')

print(f"\nCD3 cutoff ({CD3_QUANTILE} quantile): {cd3_cutoff:.6f}")
print(f"CD3 groups: {df['cd3_group'].value_counts().to_dict()}")

# Create 4-group variable
df['prolif_cd3'] = df['prolif_group'] + '_' + df['cd3_group']
print(f"\n4-group distribution:")
for group in ['high_high', 'high_low', 'low_high', 'low_low']:
    n = (df['prolif_cd3'] == group).sum()
    print(f"  {group}: {n} patients")

# ============================================================================
# STEP 4: STRATIFIED SURVIVAL ANALYSIS
# ============================================================================
print("\n[4/5] Performing stratified survival analysis...")

def analyze_stratified_survival(df_input, prolif_status, outcome='OS'):
    """Analyze survival stratified by proliferation status"""

    time_col = f'{outcome}_months'
    event_col = f'{outcome}_event'

    subset = df_input[df_input['prolif_group'] == prolif_status].copy()
    subset = subset.dropna(subset=[time_col, event_col, 'cd3_group'])

    if len(subset) < 10:
        return None

    n_patients = len(subset)
    n_events = int(subset[event_col].sum())

    results = {
        'prolif_status': prolif_status,
        'outcome': outcome,
        'n_patients': n_patients,
        'n_events': n_events
    }

    # Log-rank test
    try:
        lr_result = multivariate_logrank_test(
            subset[time_col],
            subset['cd3_group'],
            subset[event_col]
        )
        results['logrank_p'] = lr_result.p_value
    except Exception as e:
        print(f"  Log-rank failed for {prolif_status}: {e}")
        results['logrank_p'] = np.nan

    # Univariate Cox
    try:
        subset_cox = subset[[time_col, event_col, 'cd3_group']].copy()
        subset_cox['cd3_high'] = (subset_cox['cd3_group'] == 'high').astype(int)

        cph = CoxPHFitter(penalizer=0.1)
        cph.fit(subset_cox[[time_col, event_col, 'cd3_high']],
                duration_col=time_col, event_col=event_col)

        results['cox_HR'] = cph.hazard_ratios_['cd3_high']
        results['cox_p'] = cph.summary.loc['cd3_high', 'p']
        results['cox_CI_lower'] = np.exp(cph.summary.loc['cd3_high', 'coef lower 95%'])
        results['cox_CI_upper'] = np.exp(cph.summary.loc['cd3_high', 'coef upper 95%'])
    except Exception as e:
        print(f"  Cox PH failed for {prolif_status}: {e}")
        results.update({'cox_HR': np.nan, 'cox_p': np.nan,
                       'cox_CI_lower': np.nan, 'cox_CI_upper': np.nan})

    return results

# Analyze for OS
print("\nOVERALL SURVIVAL (OS):")
print("-" * 60)

os_results = []
for prolif_status in ['high', 'low']:
    result = analyze_stratified_survival(df, prolif_status, outcome='OS')
    if result:
        os_results.append(result)
        print(f"\n{prolif_status.upper()} proliferation:")
        print(f"  N={result['n_patients']}, Events={result['n_events']}")
        print(f"  Log-rank P={result['logrank_p']:.4f}")
        print(f"  Cox HR={result['cox_HR']:.3f} (95% CI: {result['cox_CI_lower']:.3f}-{result['cox_CI_upper']:.3f}), P={result['cox_p']:.4f}")

# Analyze for RFS (if available)
rfs_results = []
if 'RFS_months' in df.columns:
    print("\n" + "=" * 60)
    print("RECURRENCE-FREE SURVIVAL (RFS):")
    print("-" * 60)

    for prolif_status in ['high', 'low']:
        result = analyze_stratified_survival(df, prolif_status, outcome='RFS')
        if result:
            rfs_results.append(result)
            print(f"\n{prolif_status.upper()} proliferation:")
            print(f"  N={result['n_patients']}, Events={result['n_events']}")
            print(f"  Log-rank P={result['logrank_p']:.4f}")
            print(f"  Cox HR={result['cox_HR']:.3f} (95% CI: {result['cox_CI_lower']:.3f}-{result['cox_CI_upper']:.3f}), P={result['cox_p']:.4f}")

# 4-group analysis
print("\n" + "=" * 60)
print("4-GROUP ANALYSIS (Proliferation × CD3):")
print("-" * 60)

four_group_results = []

for outcome, time_col, event_col in [('OS', 'OS_months', 'OS_event'),
                                      ('RFS', 'RFS_months', 'RFS_event')]:
    if time_col not in df.columns:
        continue

    subset = df.dropna(subset=[time_col, event_col, 'prolif_cd3'])

    # Log-rank test
    try:
        lr_result = multivariate_logrank_test(
            subset[time_col],
            subset['prolif_cd3'],
            subset[event_col]
        )
        lr_pval = lr_result.p_value
    except:
        lr_pval = np.nan

    print(f"\n{outcome}:")
    print(f"  4-group log-rank P={lr_pval:.4f}")

    # Group statistics
    for group in ['high_high', 'high_low', 'low_high', 'low_low']:
        grp = subset[subset['prolif_cd3'] == group]
        if len(grp) > 0:
            kmf = KaplanMeierFitter()
            kmf.fit(grp[time_col], grp[event_col])

            four_group_results.append({
                'outcome': outcome,
                'group': group,
                'n': len(grp),
                'events': int(grp[event_col].sum()),
                'median_survival': kmf.median_survival_time_,
                'logrank_p': lr_pval
            })

            print(f"  {group}: n={len(grp)}, events={int(grp[event_col].sum())}, median={kmf.median_survival_time_:.1f}")

# Save results
df_os_results = pd.DataFrame(os_results)
df_os_results.to_csv(OUTPUT_DIR / 'ERpos_proliferation_stratified_survival.csv', index=False)
print(f"\nSaved: {OUTPUT_DIR / 'ERpos_proliferation_stratified_survival.csv'}")

if rfs_results:
    df_rfs_results = pd.DataFrame(rfs_results)
    df_combined = pd.concat([df_os_results, df_rfs_results], ignore_index=True)
    df_combined.to_csv(OUTPUT_DIR / 'ERpos_proliferation_stratified_survival.csv', index=False)

df_4group = pd.DataFrame(four_group_results)
df_4group.to_csv(OUTPUT_DIR / 'ERpos_4group_survival.csv', index=False)
print(f"Saved: {OUTPUT_DIR / 'ERpos_4group_survival.csv'}")

# ============================================================================
# STEP 5: GENERATE KAPLAN-MEIER PLOTS
# ============================================================================
print("\n[5/5] Generating Kaplan-Meier curves...")

def plot_km_stratified(df_input, outcome='OS'):
    """Plot KM curves for stratified analysis"""

    time_col = f'{outcome}_months'
    event_col = f'{outcome}_event'

    fig, axes = plt.subplots(1, 3, figsize=(18, 5), dpi=300)

    # Plot 1: High proliferation
    ax = axes[0]
    subset = df_input[df_input['prolif_group'] == 'high'].dropna(subset=[time_col, event_col, 'cd3_group'])

    if len(subset) > 0:
        kmf = KaplanMeierFitter()
        colors = {'high': 'green', 'low': 'red'}

        for group in ['high', 'low']:
            grp = subset[subset['cd3_group'] == group]
            if len(grp) > 0:
                kmf.fit(grp[time_col], grp[event_col], label=f'CD3 {group} (n={len(grp)})')
                kmf.plot_survival_function(ax=ax, ci_show=True, color=colors[group])

        try:
            lr = multivariate_logrank_test(subset[time_col], subset['cd3_group'], subset[event_col])
            p_val = lr.p_value
        except:
            p_val = np.nan

        ax.set_title(f'High Proliferation\nLog-rank P={p_val:.4f}', fontsize=12, fontweight='bold')
        ax.set_xlabel(f'{outcome} (months)', fontsize=11)
        ax.set_ylabel('Survival Probability', fontsize=11)
        ax.legend(loc='lower left', fontsize=9)
        ax.grid(alpha=0.3)

    # Plot 2: Low proliferation
    ax = axes[1]
    subset = df_input[df_input['prolif_group'] == 'low'].dropna(subset=[time_col, event_col, 'cd3_group'])

    if len(subset) > 0:
        kmf = KaplanMeierFitter()
        colors = {'high': 'green', 'low': 'red'}

        for group in ['high', 'low']:
            grp = subset[subset['cd3_group'] == group]
            if len(grp) > 0:
                kmf.fit(grp[time_col], grp[event_col], label=f'CD3 {group} (n={len(grp)})')
                kmf.plot_survival_function(ax=ax, ci_show=True, color=colors[group])

        try:
            lr = multivariate_logrank_test(subset[time_col], subset['cd3_group'], subset[event_col])
            p_val = lr.p_value
        except:
            p_val = np.nan

        ax.set_title(f'Low Proliferation\nLog-rank P={p_val:.4f}', fontsize=12, fontweight='bold')
        ax.set_xlabel(f'{outcome} (months)', fontsize=11)
        ax.set_ylabel('Survival Probability', fontsize=11)
        ax.legend(loc='lower left', fontsize=9)
        ax.grid(alpha=0.3)

    # Plot 3: 4-group analysis
    ax = axes[2]
    subset = df_input.dropna(subset=[time_col, event_col, 'prolif_cd3'])

    if len(subset) > 0:
        kmf = KaplanMeierFitter()
        colors = {'high_high': 'darkgreen', 'high_low': 'red',
                 'low_high': 'blue', 'low_low': 'orange'}

        for group in ['high_high', 'high_low', 'low_high', 'low_low']:
            grp = subset[subset['prolif_cd3'] == group]
            if len(grp) > 0:
                label = group.replace('_', ' / ')
                kmf.fit(grp[time_col], grp[event_col], label=f'{label} (n={len(grp)})')
                kmf.plot_survival_function(ax=ax, ci_show=False, color=colors.get(group, 'gray'))

        try:
            lr = multivariate_logrank_test(subset[time_col], subset['prolif_cd3'], subset[event_col])
            p_val = lr.p_value
        except:
            p_val = np.nan

        ax.set_title(f'4-Group (Prolif/CD3)\nLog-rank P={p_val:.4f}', fontsize=12, fontweight='bold')
        ax.set_xlabel(f'{outcome} (months)', fontsize=11)
        ax.set_ylabel('Survival Probability', fontsize=11)
        ax.legend(loc='lower left', fontsize=8)
        ax.grid(alpha=0.3)

    plt.suptitle(f'Proliferation-Stratified Survival Analysis ({outcome})',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / f'KM_proliferation_stratified_{outcome}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: KM_proliferation_stratified_{outcome}.png")

# Generate plots
plot_km_stratified(df, outcome='OS')

if 'RFS_months' in df.columns and df['RFS_months'].notna().sum() > 10:
    plot_km_stratified(df, outcome='RFS')

# ============================================================================
# SUMMARY
# ============================================================================
total_time = __import__('time').time() - t0_total

print("\n" + "=" * 80)
print("STEP 11 COMPLETE")
print("=" * 80)

print("\nKey Findings:")
if os_results:
    high_result = [r for r in os_results if r['prolif_status'] == 'high'][0]
    low_result = [r for r in os_results if r['prolif_status'] == 'low'][0]

    print(f"\nHigh proliferation tumors:")
    print(f"  CD3+ T cells: HR={high_result['cox_HR']:.3f}, P={high_result['cox_p']:.4f}")
    if high_result['cox_p'] < 0.05:
        print(f"  → CD3 IS PROGNOSTIC ✓")
    else:
        print(f"  → CD3 not significant")

    print(f"\nLow proliferation tumors:")
    print(f"  CD3+ T cells: HR={low_result['cox_HR']:.3f}, P={low_result['cox_p']:.4f}")
    if low_result['cox_p'] < 0.05:
        print(f"  → CD3 IS PROGNOSTIC")
    else:
        print(f"  → CD3 not significant ✓")

print("\nBiological Interpretation:")
print("  High-prolif tumors: More immunogenic, immune surveillance active")
print("  Low-prolif tumors: Hormone-driven, less dependent on immune evasion")

print(f"\nOutputs saved to:")
print(f"  - {OUTPUT_DIR / 'ERpos_proliferation_stratified_survival.csv'}")
print(f"  - {OUTPUT_DIR / 'ERpos_4group_survival.csv'}")
print(f"  - {FIGURES_DIR / 'KM_proliferation_stratified_*.png'}")

print(f"\nTotal time: {total_time:.1f}s")
print("=" * 80)
