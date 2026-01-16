#!/usr/bin/env python3
"""
Step 5 Postprocessing (No Squidpy Required)

Creates curated, teaching-friendly summaries from the Step 5 enrichment significance tables.

Inputs (must already exist):
  - analysis_imc/ERpos_nhood_enrichment_r{25,40}_significance.csv

Outputs (per radius):
  - analysis_imc/ERpos_nhood_enrichment_r{radius}_significance_annotated.csv
  - analysis_imc/ERpos_nhood_enrichment_r{radius}_group_pair_summary.csv
  - analysis_imc/ERpos_nhood_enrichment_r{radius}_curated_hypotheses_summary.csv
  - analysis_imc/ERpos_nhood_enrichment_top50_significant_r{radius}.csv
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def celltype_group(name: str) -> str:
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


def summarize_subset(sub: pd.DataFrame, label: str) -> dict:
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


def process_radius(radius: int) -> None:
    in_file = Path(f"analysis_imc/ERpos_nhood_enrichment_r{radius}_significance.csv")
    if not in_file.exists():
        print(f"Skip r={radius}: missing {in_file}")
        return

    n_rois_file = Path(f"analysis_imc/ERpos_nhood_enrichment_r{radius}_n_rois.csv")
    if not n_rois_file.exists():
        print(f"Skip r={radius}: missing {n_rois_file} (ROI support matrix)")
        return

    df = pd.read_csv(in_file)
    required = {"celltype", "neighbor", "z", "FDR", "significant_FDR"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{in_file} missing columns: {sorted(missing)}")

    df_n_rois = pd.read_csv(n_rois_file, index_col=0)
    # Standardize celltype strings for lookup safety
    df["celltype"] = df["celltype"].astype(str)
    df["neighbor"] = df["neighbor"].astype(str)

    df["origin_group"] = df["celltype"].map(celltype_group)
    df["neighbor_group"] = df["neighbor"].map(celltype_group)

    def _lookup_support(row) -> float:
        ct = row["celltype"]
        nb = row["neighbor"]
        try:
            return float(df_n_rois.loc[ct, nb])
        except Exception:
            return np.nan

    df["n_rois_support"] = df.apply(_lookup_support, axis=1)

    out_annot = Path(f"analysis_imc/ERpos_nhood_enrichment_r{radius}_significance_annotated.csv")
    df.to_csv(out_annot, index=False)

    group_summary = (
        df.groupby(["origin_group", "neighbor_group"], dropna=False)
        .agg(
            pairs_total=("z", "size"),
            mean_z=("z", "mean"),
            median_z=("z", "median"),
            pairs_sig_fdr=("significant_FDR", "sum"),
            pairs_sig_enriched_z_gt_0=("z", lambda s: int(((df.loc[s.index, "significant_FDR"] == True) & (s > 0)).sum())),
            pairs_sig_depleted_z_lt_0=("z", lambda s: int(((df.loc[s.index, "significant_FDR"] == True) & (s < 0)).sum())),
        )
        .reset_index()
        .sort_values(["pairs_sig_fdr", "pairs_total"], ascending=[False, False])
    )
    out_group = Path(f"analysis_imc/ERpos_nhood_enrichment_r{radius}_group_pair_summary.csv")
    group_summary.to_csv(out_group, index=False)

    stromal_to_epithelial = df[(df["origin_group"] == "stromal") & (df["neighbor_group"] == "epithelial")]
    epithelial_to_stromal = df[(df["origin_group"] == "epithelial") & (df["neighbor_group"] == "stromal")]
    fibro_like_to_epithelial = df[df["celltype"].map(lambda x: ("fibroblast" in str(x).lower()) or ("pericyte" in str(x).lower())) & (df["neighbor_group"] == "epithelial")]
    epithelial_to_fibro_like = df[(df["origin_group"] == "epithelial") & df["neighbor"].map(lambda x: ("fibroblast" in str(x).lower()) or ("pericyte" in str(x).lower()))]
    immune_to_epithelial = df[(df["origin_group"] == "immune") & (df["neighbor_group"] == "epithelial")]
    epithelial_to_immune = df[(df["origin_group"] == "epithelial") & (df["neighbor_group"] == "immune")]
    immune_to_immune = df[(df["origin_group"] == "immune") & (df["neighbor_group"] == "immune")]

    hyp_summary = pd.DataFrame(
        [
            summarize_subset(stromal_to_epithelial, "stromal->epithelial (CAF barrier / stromal adjacency)"),
            summarize_subset(epithelial_to_stromal, "epithelial->stromal (CAF barrier / stromal adjacency)"),
            summarize_subset(fibro_like_to_epithelial, "fibro_like->epithelial (fibro/pericyte adjacency)"),
            summarize_subset(epithelial_to_fibro_like, "epithelial->fibro_like (fibro/pericyte adjacency)"),
            summarize_subset(immune_to_epithelial, "immune->epithelial (immune exclusion / infiltration)"),
            summarize_subset(epithelial_to_immune, "epithelial->immune (immune exclusion / infiltration)"),
            summarize_subset(immune_to_immune, "immune->immune (immune aggregates / TLS-like)"),
        ]
    )
    out_hyp = Path(f"analysis_imc/ERpos_nhood_enrichment_r{radius}_curated_hypotheses_summary.csv")
    hyp_summary.to_csv(out_hyp, index=False)

    top_sig = (
        df[df["significant_FDR"] == True]
        .assign(abs_z=lambda d: d["z"].abs())
        .sort_values("abs_z", ascending=False)
        .head(50)
        .drop(columns=["abs_z"])
    )
    out_top = Path(f"analysis_imc/ERpos_nhood_enrichment_top50_significant_r{radius}.csv")
    top_sig.to_csv(out_top, index=False)

    # ER+ breast cancer curated pairs (teaching / interpretation-friendly)
    # We focus on hypotheses that tend to be meaningful in ER+ contexts and are visible at contact scales.
    MIN_ROI_SUPPORT = 25
    sig = df[(df["significant_FDR"] == True) & (df["n_rois_support"] >= MIN_ROI_SUPPORT)].copy()

    def is_fibro_like(x: str) -> bool:
        s = str(x).lower()
        return ("fibroblast" in s) or ("pericyte" in s)

    hypotheses = {
        "Tumor compartmentalization (epithelial↔epithelial)": sig[(sig["origin_group"] == "epithelial") & (sig["neighbor_group"] == "epithelial")],
        "Stromal barrier vs epithelium (stromal↔epithelial)": sig[
            ((sig["origin_group"] == "stromal") & (sig["neighbor_group"] == "epithelial"))
            | ((sig["origin_group"] == "epithelial") & (sig["neighbor_group"] == "stromal"))
        ],
        "Fibro/pericyte vs epithelium (fibro_like↔epithelial)": sig[
            (sig["celltype"].map(is_fibro_like) & (sig["neighbor_group"] == "epithelial"))
            | ((sig["origin_group"] == "epithelial") & sig["neighbor"].map(is_fibro_like))
        ],
        "Immune exclusion/infiltration (immune↔epithelial)": sig[
            ((sig["origin_group"] == "immune") & (sig["neighbor_group"] == "epithelial"))
            | ((sig["origin_group"] == "epithelial") & (sig["neighbor_group"] == "immune"))
        ],
        "Immune aggregates (immune↔immune)": sig[(sig["origin_group"] == "immune") & (sig["neighbor_group"] == "immune")],
        "Vascular niche (endothelial↔pericyte/stroma/tumor)": sig[
            (sig["origin_group"] == "endothelial")
            | (sig["neighbor_group"] == "endothelial")
        ],
        "Luminal ER+ tumor interface (Luminal ER+ ↔ stroma/fibro)": sig[
            ((sig["celltype"] == "Luminal ER+ tumor cells") & (sig["neighbor_group"].isin(["stromal", "endothelial"])))
            | ((sig["neighbor"] == "Luminal ER+ tumor cells") & (sig["origin_group"].isin(["stromal", "endothelial"])))
        ],
        "Luminal ER+ tumor immune proximity (Luminal ER+ ↔ immune)": sig[
            ((sig["celltype"] == "Luminal ER+ tumor cells") & (sig["neighbor_group"] == "immune"))
            | ((sig["neighbor"] == "Luminal ER+ tumor cells") & (sig["origin_group"] == "immune"))
        ],
    }

    rows = []
    for hyp, sub in hypotheses.items():
        if len(sub) == 0:
            continue
        sub_pos = sub[sub["z"] > 0].sort_values("z", ascending=False).head(10)
        sub_neg = sub[sub["z"] < 0].sort_values("z", ascending=True).head(10)
        for direction, df_top in [("enriched", sub_pos), ("depleted", sub_neg)]:
            for _, r in df_top.iterrows():
                rows.append(
                    {
                        "radius_um": radius,
                        "hypothesis": hyp,
                        "direction": direction,
                        "celltype": r["celltype"],
                        "neighbor": r["neighbor"],
                        "z": float(r["z"]),
                        "FDR": float(r["FDR"]),
                        "n_rois_support": float(r["n_rois_support"]),
                        "pair": r.get("pair", f"{r['celltype']} → {r['neighbor']}"),
                    }
                )

    out_curated_pairs = Path(f"analysis_imc/ERpos_nhood_enrichment_r{radius}_erplus_curated_pairs.csv")
    pd.DataFrame(rows).to_csv(out_curated_pairs, index=False)

    print(f"r={radius}: wrote {out_annot}, {out_group}, {out_hyp}, {out_top}")
    print(f"r={radius}: wrote {out_curated_pairs} (MIN_ROI_SUPPORT={MIN_ROI_SUPPORT})")


def main() -> None:
    for radius in (25, 40):
        process_radius(radius)

    # Cross-radius stability: pairs that are FDR-significant in both radii and have consistent sign.
    a25 = Path("analysis_imc/ERpos_nhood_enrichment_r25_significance_annotated.csv")
    a40 = Path("analysis_imc/ERpos_nhood_enrichment_r40_significance_annotated.csv")
    if a25.exists() and a40.exists():
        df25 = pd.read_csv(a25)
        df40 = pd.read_csv(a40)
        keep25 = df25[df25["significant_FDR"] == True][["celltype", "neighbor", "z", "FDR", "n_rois_support", "origin_group", "neighbor_group"]].copy()
        keep40 = df40[df40["significant_FDR"] == True][["celltype", "neighbor", "z", "FDR", "n_rois_support", "origin_group", "neighbor_group"]].copy()
        merged = keep25.merge(keep40, on=["celltype", "neighbor", "origin_group", "neighbor_group"], suffixes=("_r25", "_r40"))
        merged["sign_consistent"] = np.sign(merged["z_r25"]) == np.sign(merged["z_r40"])
        merged = merged.sort_values(["sign_consistent", "z_r25"], ascending=[False, False])
        out_stable = Path("analysis_imc/ERpos_nhood_enrichment_stable_pairs_r25_r40.csv")
        merged.to_csv(out_stable, index=False)
        print(f"Stable cross-radius pairs: wrote {out_stable} (rows={len(merged)})")


if __name__ == "__main__":
    main()
