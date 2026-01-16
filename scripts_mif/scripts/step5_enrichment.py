
#!/usr/bin/env python
"""
Neighborhood enrichment (z-score) for ER+ IMC using Squidpy.

ISOLATED PIPELINE - Run from pipeline_leiden_validated/ directory:
  cd /path/to/pipeline_leiden_validated
  conda run -p /path/to/squidpy-env python scripts/step5_enrichment.py

Input:
  data/imc_ERpos_with_embeddings_and_full_names.h5ad
Outputs:
  analysis_imc/ERpos_nhood_enrichment_r25.csv/.png/.top20.csv
  analysis_imc/ERpos_nhood_enrichment_r40.csv/.png/.top20.csv
  analysis_imc/ERpos_nhood_enrichment_r{radius}_n_rois.csv (conteo de ROIs por par)
  analysis_imc/ERpos_nhood_enrichment_per_patient_r{radius}.csv (features por paciente; para Step 8)
Notes:
  - Enrichment por slide_scene; sin aristas entre ROIs.
  - Radios del paper: 25µm, 40µm (solo radius, sin n_neighs cap).
  - Filtra celltypes raros (<200 celdas), alinea categorías y reindexa al set global.
"""
import os
import sys
import time
import argparse
from pathlib import Path

# Ensure Numba stays deterministic + avoids cache permission issues.
# The pipeline runner also sets these, but keeping this here makes standalone runs safer.
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp")

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq

MIN_CELLS = 200
RADII = (25, 40)  # Both radii as per paper methods


def _per_patient_pairs_of_interest():
    # These are computed from each patient's averaged ROI enrichment matrix (z-scores).
    return [
        ("CD3+ T lymphocytes", "Luminal tumor cells", "Tcell_tumor_infiltration"),
        ("CD3+ T lymphocytes", "HER2+ ER+ tumor cells", "Tcell_ERposHER2pos_contact"),
        ("Vimentin+ fibroblasts", "Vimentin+ fibroblasts", "CAF_homotypic_clustering"),
        ("SMA+ pericytes and fibroblasts", "SMA+ pericytes and fibroblasts", "pericyte_clustering"),
        ("Luminal tumor cells", "CD31+ endothelial cells", "tumor_angiogenesis"),
        ("CD31+ endothelial cells", "CD31+ endothelial cells", "vascular_clustering"),
        ("Basal tumor cells", "Basal tumor cells", "basal_clustering"),
        ("Luminal tumor cells", "Luminal tumor cells", "luminal_clustering"),
        ("CD68+ macrophages", "Vimentin+ fibroblasts", "macrophage_CAF_interaction"),
        ("CD3+ T lymphocytes", "CD68+ macrophages", "T_cell_macrophage"),
        ("Cytokeratin-low tumor cells", "Cytokeratin-low tumor cells", "CKlow_clustering"),
        ("Quiescent stroma", "Quiescent stroma", "quiescent_stroma_clustering"),
        ("CD3+ T lymphocytes", "Basal tumor cells", "Tcell_basal_contact"),
        ("CD68+ macrophages", "Luminal tumor cells", "macrophage_tumor_contact"),
    ]


def run_enrichment(radius, compute_global=True, compute_per_patient=True):
    t0 = time.time()
    adata = sc.read_h5ad("data/imc_ERpos_with_embeddings_and_full_names.h5ad")
    # FIX: Use fine-grained leiden annotations (20 phenotypes) instead of coarse celltype (4 types)
    adata.obs["celltype"] = adata.obs["leiden"].astype("category")

    # Create cell type name mapping from h5ad (short → full names)
    name_mapping = adata.obs[['celltype_abbrev', 'celltype_full']].drop_duplicates()
    name_map_dict = dict(zip(name_mapping['celltype_abbrev'], name_mapping['celltype_full']))

    # Ensure output dirs exist
    OUTPUT_DIR = Path("analysis_imc")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Filtrar celltypes raros globalmente y alinear categorías
    counts = adata.obs["celltype"].value_counts()
    keep = counts[counts >= MIN_CELLS].index
    adata = adata[adata.obs["celltype"].isin(keep)].copy()
    adata.obs["celltype"] = adata.obs["celltype"].cat.remove_unused_categories()
    celltypes = adata.obs["celltype"].cat.categories.tolist()
    if len(celltypes) < 2:
        print(f"No celltypes with >= {MIN_CELLS} cells for radius {radius}")
        return

    z_mats = []
    support_mats = []

    # Per-patient accumulators (sum + counts over ROIs)
    patient_sum = {}
    patient_n = {}
    pairs_of_interest = _per_patient_pairs_of_interest()

    scenes = adata.obs["slide_scene"].unique()
    total = len(scenes)
    print(f"\n{'='*60}", flush=True)
    print(f"Processing {total} ROIs with radius={radius}µm", flush=True)
    print(f"{'='*60}", flush=True)

    roi_start_time = time.time()
    for idx, scene in enumerate(scenes, 1):
        sub = adata[adata.obs["slide_scene"] == scene].copy()
        if sub.n_obs == 0:
            continue
        sub.obs["celltype"] = sub.obs["celltype"].cat.remove_unused_categories()
        present = sub.obs["celltype"].cat.categories.tolist()
        if len(present) < 2:
            continue

        sq.gr.spatial_neighbors(
            sub,
            coord_type="generic",
            radius=radius,
            delaunay=False,
        )
        sq.gr.nhood_enrichment(sub, cluster_key="celltype")

        z = sub.uns["celltype_nhood_enrichment"]["zscore"]
        z_df = pd.DataFrame(z, index=present, columns=present)
        z_df = z_df.reindex(index=celltypes, columns=celltypes)

        if compute_global:
            z_mats.append(z_df)
            support_mats.append(~z_df.isna())

        if compute_per_patient:
            patient_id = str(sub.obs["Patient"].iloc[0])
            z_arr = z_df.to_numpy()
            finite = np.isfinite(z_arr)
            if patient_id not in patient_sum:
                patient_sum[patient_id] = np.zeros_like(z_arr, dtype=float)
                patient_n[patient_id] = np.zeros_like(z_arr, dtype=np.int32)
            patient_sum[patient_id] += np.where(finite, z_arr, 0.0)
            patient_n[patient_id] += finite.astype(np.int32)

        # Progress logging every 10 ROIs
        if idx % 10 == 0 or idx == total:
            elapsed = time.time() - roi_start_time
            rate = idx / elapsed if elapsed > 0 else 0
            remaining = (total - idx) / rate if rate > 0 else 0
            pct = 100 * idx / total
            print(f"[r={radius}µm] {idx}/{total} ({pct:.1f}%) | "
                  f"{elapsed/60:.1f}min elapsed | "
                  f"ETA: {remaining/60:.1f}min | "
                  f"Rate: {rate:.2f} ROI/s", flush=True)

    if compute_global and not z_mats:
        print(f"No scenes processed for radius {radius}")
        return

    df = None
    support = None
    df_renamed = None

    if compute_global:
        df = pd.concat(z_mats).groupby(level=0).mean()
        support = pd.concat(support_mats).groupby(level=0).sum()

        # Rename cell types to full names for better readability
        df_renamed = df.rename(index=name_map_dict, columns=name_map_dict)
        support_renamed = support.rename(index=name_map_dict, columns=name_map_dict)

        df_renamed.to_csv(OUTPUT_DIR / f"ERpos_nhood_enrichment_r{radius}.csv")
        support_renamed.to_csv(OUTPUT_DIR / f"ERpos_nhood_enrichment_r{radius}_n_rois.csv")

    # Compute per-patient features (independent of global reporting)
    if compute_per_patient and patient_sum:
        print(f"\n{'='*60}", flush=True)
        print("Computing per-patient enrichment values...", flush=True)
        print(f"{'='*60}", flush=True)

        records = []
        for patient_id in sorted(patient_sum.keys(), key=lambda x: (len(x), x)):
            s = patient_sum[patient_id]
            n = patient_n[patient_id]
            with np.errstate(invalid="ignore", divide="ignore"):
                avg = s / np.where(n > 0, n, np.nan)
            avg_df = pd.DataFrame(avg, index=celltypes, columns=celltypes)
            avg_full = avg_df.rename(index=name_map_dict, columns=name_map_dict)

            rec = {"Patient": patient_id}
            for src, tgt, feat in pairs_of_interest:
                rec[feat] = float(avg_full.loc[src, tgt]) if (src in avg_full.index and tgt in avg_full.columns) else np.nan
            records.append(rec)

        df_patient = pd.DataFrame(records)
        out_path = OUTPUT_DIR / f"ERpos_nhood_enrichment_per_patient_r{radius}.csv"
        df_patient.to_csv(out_path, index=False)

        print(f"\n✓ Per-patient enrichment saved: {out_path}", flush=True)
        print(f"  Shape: {df_patient.shape} (patients × features)", flush=True)
        print(f"  Features: {len(pairs_of_interest)} cell-cell interaction pairs", flush=True)
        print(f"  Missing values per feature:", flush=True)
        for col in df_patient.columns:
            if col == "Patient":
                continue
            n_missing = int(df_patient[col].isna().sum())
            pct_missing = 100.0 * n_missing / len(df_patient)
            print(f"    {col}: {n_missing}/{len(df_patient)} ({pct_missing:.1f}%)", flush=True)

    if not compute_global:
        elapsed = time.time() - t0
        print(f"\n{'='*60}", flush=True)
        print(f"✓ COMPLETED r={radius}µm in {elapsed/60:.1f} min (per-patient={compute_per_patient})", flush=True)
        print(f"{'='*60}\n", flush=True)
        return

    # Create figures directory (global reporting)
    FIGURES_DIR = Path("figures/step5_enrichment")
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # ========================================================================
    # PLOT 1: Enrichment heatmap (using full cell type names)
    # ========================================================================
    plt.figure(figsize=(14, 12), dpi=300)  # Larger to accommodate full names
    sns.heatmap(df_renamed, cmap="coolwarm", center=0, cbar_kws={'label': 'Z-score', 'shrink': 0.5})
    plt.title(f"Neighborhood Enrichment Z-scores (r={radius}μm)", fontsize=13, fontweight='bold')
    plt.xlabel('Neighbor Cell Type', fontsize=11)
    plt.ylabel('Origin Cell Type', fontsize=11)
    plt.xticks(rotation=45, ha='right', fontsize=7)  # Smaller font for full names
    plt.yticks(rotation=0, fontsize=7)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / f'01_enrichment_heatmap_r{radius}.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Ensure index has a name for melting
    df_renamed.index.name = "celltype" if df_renamed.index.name is None else df_renamed.index.name

    # Prepare melted data for additional plots (using full names)
    df_melted = (
        df_renamed.reset_index()
        .melt(id_vars="celltype", var_name="neighbor", value_name="z")
        .dropna()
        .assign(pair=lambda d: d["celltype"].astype(str) + " → " + d["neighbor"].astype(str))
    )

    # Calculate p-values and FDR correction for statistical significance
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    df_melted['p_value'] = 2 * (1 - stats.norm.cdf(np.abs(df_melted['z'])))

    # FDR correction (Benjamini-Hochberg)
    reject_fdr, pvals_fdr, _, _ = multipletests(
        df_melted['p_value'].values,
        method='fdr_bh',
        alpha=0.05
    )
    df_melted['FDR'] = pvals_fdr
    df_melted['significant_FDR'] = reject_fdr

    # Significance levels
    df_melted['significance'] = 'ns'
    df_melted.loc[np.abs(df_melted['z']) > 1.96, 'significance'] = '*'
    df_melted.loc[np.abs(df_melted['z']) > 2.58, 'significance'] = '**'
    df_melted.loc[np.abs(df_melted['z']) > 3.29, 'significance'] = '***'
    df_melted.loc[df_melted['significant_FDR'], 'significance'] = 'FDR<0.05'

    # Save significance table
    df_melted.to_csv(f"analysis_imc/ERpos_nhood_enrichment_r{radius}_significance.csv", index=False)

    # ------------------------------------------------------------------------
    # Postprocessing: curated biological questions + summaries
    #
    # These are derived from the Step 5 enrichment/significance table and help
    # students (and paper validation) find signals that would otherwise be buried
    # in the full 18×18 matrix.
    # ------------------------------------------------------------------------

    def _celltype_group(name: str) -> str:
        s = str(name).lower()
        # Keep immune matching strict to avoid false positives (e.g., CD44 matching CD4).
        if (
            ("cd3+" in s)
            or ("cd20+" in s)
            or ("cd68+" in s)
            or ("t lymph" in s)
            or ("b lymph" in s)
            or ("macroph" in s)
        ):
            return "immune"
        if ("cd31+" in s) or ("endothelial" in s):
            return "endothelial"
        if ("fibroblast" in s) or ("pericyte" in s) or ("stroma" in s):
            return "stromal"
        if ("tumor" in s) or ("luminal" in s) or ("basal" in s) or ("myoepithelial" in s):
            return "epithelial"
        return "other"

    df_annot = df_melted.copy()
    df_annot["origin_group"] = df_annot["celltype"].map(_celltype_group)
    df_annot["neighbor_group"] = df_annot["neighbor"].map(_celltype_group)
    df_annot.to_csv(
        f"analysis_imc/ERpos_nhood_enrichment_r{radius}_significance_annotated.csv",
        index=False,
    )

    def _summarize_subset(sub: pd.DataFrame, label: str) -> dict:
        sig_sub = sub[sub["significant_FDR"] == True]
        return {
            "label": label,
            "pairs_total": int(len(sub)),
            "pairs_sig_fdr": int(len(sig_sub)),
            "pairs_sig_enriched_z_gt_0": int((sig_sub["z"] > 0).sum()),
            "pairs_sig_depleted_z_lt_0": int((sig_sub["z"] < 0).sum()),
            "min_z": float(sub["z"].min()) if len(sub) else np.nan,
            "median_z": float(sub["z"].median()) if len(sub) else np.nan,
            "max_z": float(sub["z"].max()) if len(sub) else np.nan,
        }

    df_group_summary = (
        df_annot.groupby(["origin_group", "neighbor_group"], dropna=False)
        .agg(
            pairs_total=("z", "size"),
            mean_z=("z", "mean"),
            median_z=("z", "median"),
            pairs_sig_fdr=("significant_FDR", "sum"),
            pairs_sig_enriched_z_gt_0=("z", lambda s: int(((df_annot.loc[s.index, "significant_FDR"] == True) & (s > 0)).sum())),
            pairs_sig_depleted_z_lt_0=("z", lambda s: int(((df_annot.loc[s.index, "significant_FDR"] == True) & (s < 0)).sum())),
        )
        .reset_index()
        .sort_values(["pairs_sig_fdr", "pairs_total"], ascending=[False, False])
    )
    df_group_summary.to_csv(
        f"analysis_imc/ERpos_nhood_enrichment_r{radius}_group_pair_summary.csv",
        index=False,
    )

    def _is_fibro_like(name: str) -> bool:
        s = str(name).lower()
        return ("fibroblast" in s) or ("pericyte" in s)

    # Curated paper-aligned questions (derived from Step 5 outputs)
    # - "stromal" includes Quiescent stroma (broad barrier/compartment signal)
    # - "fibro_like" focuses on fibroblast/pericyte programs (closer to CAF narratives)
    stromal_to_epithelial = df_annot[(df_annot["origin_group"] == "stromal") & (df_annot["neighbor_group"] == "epithelial")]
    epithelial_to_stromal = df_annot[(df_annot["origin_group"] == "epithelial") & (df_annot["neighbor_group"] == "stromal")]
    fibro_like_to_epithelial = df_annot[df_annot["celltype"].map(_is_fibro_like) & (df_annot["neighbor_group"] == "epithelial")]
    epithelial_to_fibro_like = df_annot[(df_annot["origin_group"] == "epithelial") & df_annot["neighbor"].map(_is_fibro_like)]
    immune_to_epithelial = df_annot[(df_annot["origin_group"] == "immune") & (df_annot["neighbor_group"] == "epithelial")]
    epithelial_to_immune = df_annot[(df_annot["origin_group"] == "epithelial") & (df_annot["neighbor_group"] == "immune")]
    immune_to_immune = df_annot[(df_annot["origin_group"] == "immune") & (df_annot["neighbor_group"] == "immune")]

    hyp_summary = pd.DataFrame(
        [
            _summarize_subset(stromal_to_epithelial, "stromal->epithelial (CAF barrier / stromal adjacency)"),
            _summarize_subset(epithelial_to_stromal, "epithelial->stromal (CAF barrier / stromal adjacency)"),
            _summarize_subset(fibro_like_to_epithelial, "fibro_like->epithelial (fibro/pericyte adjacency)"),
            _summarize_subset(epithelial_to_fibro_like, "epithelial->fibro_like (fibro/pericyte adjacency)"),
            _summarize_subset(immune_to_epithelial, "immune->epithelial (immune exclusion / infiltration)"),
            _summarize_subset(epithelial_to_immune, "epithelial->immune (immune exclusion / infiltration)"),
            _summarize_subset(immune_to_immune, "immune->immune (immune aggregates / TLS-like)"),
        ]
    )
    hyp_summary.to_csv(
        f"analysis_imc/ERpos_nhood_enrichment_r{radius}_curated_hypotheses_summary.csv",
        index=False,
    )

    top_sig = (
        df_annot[df_annot["significant_FDR"] == True]
        .assign(abs_z=lambda d: d["z"].abs())
        .sort_values("abs_z", ascending=False)
        .head(50)
        .drop(columns=["abs_z"])
    )
    top_sig.to_csv(
        f"analysis_imc/ERpos_nhood_enrichment_top50_significant_r{radius}.csv",
        index=False,
    )

    print(f"\n  Curated Step 5 summaries (r={radius}µm):")
    print(f"    - analysis_imc/ERpos_nhood_enrichment_r{radius}_significance_annotated.csv")
    print(f"    - analysis_imc/ERpos_nhood_enrichment_r{radius}_group_pair_summary.csv")
    print(f"    - analysis_imc/ERpos_nhood_enrichment_r{radius}_curated_hypotheses_summary.csv")
    print(f"    - analysis_imc/ERpos_nhood_enrichment_top50_significant_r{radius}.csv")

    # Save top 20 for backward compatibility
    top = (
        df_melted
        .assign(abs_z=lambda d: d["z"].abs())
        .sort_values("abs_z", ascending=False)
        .head(20)
    )
    top.to_csv(f"analysis_imc/ERpos_nhood_enrichment_top20_r{radius}.csv", index=False)

    # Print significance summary
    print(f"\n  📊 Statistical Significance (r={radius}µm):")
    print(f"     Total pairs: {len(df_melted)}")
    print(f"     Significant (FDR < 0.05): {df_melted['significant_FDR'].sum()} ({100*df_melted['significant_FDR'].sum()/len(df_melted):.1f}%)")
    print(f"     Enriched (Z>0, FDR<0.05): {((df_melted['z'] > 0) & df_melted['significant_FDR']).sum()}")
    print(f"     Depleted (Z<0, FDR<0.05): {((df_melted['z'] < 0) & df_melted['significant_FDR']).sum()}")

    # ========================================================================
    # PLOT 2: Top 10 enriched pairs (barplot)
    # ========================================================================
    print("Generating additional plots...")

    top_enriched = df_melted.nlargest(10, 'z').copy()

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

    # Create horizontal barplot
    bars = ax.barh(range(len(top_enriched)), top_enriched['z'].values, color='firebrick', alpha=0.8)

    ax.set_yticks(range(len(top_enriched)))
    ax.set_yticklabels(top_enriched['pair'].values, fontsize=9)
    ax.set_xlabel('Enrichment Z-score', fontsize=12)
    ax.set_title(f'Top 10 Enriched Cell-Cell Pairs (r={radius}μm)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    ax.axvline(0, color='black', linewidth=0.8)

    # Add value labels
    for i, (idx, row) in enumerate(top_enriched.iterrows()):
        ax.text(row['z'] + 0.1, i, f"{row['z']:.2f}", va='center', fontsize=8)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / f'02_top_enriched_pairs_r{radius}.png', dpi=300, bbox_inches='tight')
    plt.close()

    # ========================================================================
    # PLOT 3: Top 10 depleted pairs (barplot)
    # ========================================================================
    top_depleted = df_melted.nsmallest(10, 'z').copy()

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

    # Create horizontal barplot
    bars = ax.barh(range(len(top_depleted)), top_depleted['z'].values, color='steelblue', alpha=0.8)

    ax.set_yticks(range(len(top_depleted)))
    ax.set_yticklabels(top_depleted['pair'].values, fontsize=9)
    ax.set_xlabel('Depletion Z-score', fontsize=12)
    ax.set_title(f'Top 10 Depleted Cell-Cell Pairs (r={radius}μm)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    ax.axvline(0, color='black', linewidth=0.8)

    # Add value labels
    for i, (idx, row) in enumerate(top_depleted.iterrows()):
        ax.text(row['z'] - 0.1, i, f"{row['z']:.2f}", va='center', ha='right', fontsize=8)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / f'03_top_depleted_pairs_r{radius}.png', dpi=300, bbox_inches='tight')
    plt.close()

    elapsed = time.time() - t0
    print(f"\n{'='*60}", flush=True)
    print(f"✓ COMPLETED r={radius}µm in {elapsed/60:.1f} min ({len(z_mats)} ROIs processed)", flush=True)
    print(f"✅ Step 5 visualizations complete! (figures/step5_enrichment/)", flush=True)
    print(f"{'='*60}\n", flush=True)


def main():
    parser = argparse.ArgumentParser(description="Step 5: Neighborhood enrichment (global + per-patient).")
    parser.add_argument(
        "--mode",
        choices=["both", "global", "per-patient"],
        default="both",
        help="Which outputs to compute: global matrix/plots, per-patient features, or both.",
    )
    parser.add_argument(
        "--radii",
        nargs="*",
        type=int,
        default=list(RADII),
        help="Radii in µm to run (default: 25 40).",
    )
    args = parser.parse_args()

    compute_global = args.mode in {"both", "global"}
    compute_per_patient = args.mode in {"both", "per-patient"}

    for r in args.radii:
        run_enrichment(r, compute_global=compute_global, compute_per_patient=compute_per_patient)


if __name__ == "__main__":
    main()
