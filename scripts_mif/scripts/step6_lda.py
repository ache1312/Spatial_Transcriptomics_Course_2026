#!/usr/bin/env python
"""
Step 6: Spatial LDA (Latent Dirichlet Allocation) for Neighborhood Discovery

ISOLATED PIPELINE - Run from pipeline_leiden_validated/ directory:
  cd /path/to/pipeline_leiden_validated
  conda run -p /path/to/squidpy-env python scripts/step6_lda.py

Purpose:
  Discover recurrent spatial neighborhoods (topics) using LDA on neighbor counts.
  Each "document" is a cell, "words" are neighboring cell types.

Input:
  - data/OUTPUT_imc_ERpos_with_graphs.h5ad (with spatial graphs from Step 3)

Outputs (CSVs):
  - analysis_imc/ERpos_neighbor_counts_r25_coarse_percell.csv (neighbor counts per cell)
  - analysis_imc/ERpos_lda_topics_by_patient.csv (4 topics × patients)
  - analysis_imc/ERpos_lda_topic_components.csv (topic definitions)

Outputs (Figures):
  - figures/step6_lda/01_topic_composition_heatmap.png
  - figures/step6_lda/02_topic_weights_by_patient.png
  - figures/step6_lda/03_topic_mean_barplot.png
  - figures/step6_lda/04_topic_distribution_violin.png
  - figures/step6_lda/05_lda_k_perplexity.png
  - figures/step6_lda/06_lda_k_stability.png

Parameters:
  - N_COMPONENTS = 6 (chosen from K diagnostics; paper uses 8)
  - MAX_CELLS = 5000 (subsample for speed)
  - RADIUS = 25μm (neighborhood radius)
  - K_RANGE = 2..10 (diagnostic sweep for K selection)

Notes:
  - LDA treats each cell as a "document" with neighboring cell types as "words"
  - Topics represent recurrent spatial patterns (e.g., "immune-rich", "tumor nests")
  - Each patient gets a distribution over the 4 topics
  - NOW USES GRAPHS FROM STEP 3 (no external CSV dependency)
"""
import os
import sys
import time
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
from sklearn.decomposition import LatentDirichletAllocation
from scipy.sparse import csr_matrix

# ============================================================================
# PARAMETERS
# ============================================================================
N_COMPONENTS = 6  # Number of LDA topics (data-driven; paper uses 8)
MAX_CELLS = 5000  # Subsample cells for speed
RADIUS = 25  # Neighborhood radius in μm
RANDOM_STATE = 42
K_RANGE = list(range(2, 11))  # K diagnostics range
STABILITY_SEEDS = [42, 52, 62]

# Input/Output paths
INPUT_H5AD = "data/OUTPUT_imc_ERpos_with_graphs.h5ad"
OUTPUT_DIR = "analysis_imc"
FIGURE_DIR = "figures/step6_lda"

# Ensure output directories exist
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def compute_neighbor_counts_from_graph(adata, conn_key, celltype_col, celltypes):
    """
    Compute neighbor counts per cell from spatial connectivity graph.

    Parameters:
    -----------
    adata : AnnData
        Annotated data with spatial connectivity in obsp
    conn_key : str
        Key in adata.obsp with connectivity matrix
    celltype_col : str
        Column in adata.obs with cell type labels
    celltypes : list
        List of cell types to count

    Returns:
    --------
    pd.DataFrame with columns: cell_id, celltype1, celltype2, ...
    """
    print(f"   Computing neighbor counts from graph: {conn_key}")

    if conn_key not in adata.obsp:
        raise ValueError(f"Connectivity matrix '{conn_key}' not found in adata.obsp")

    conn_matrix = adata.obsp[conn_key]

    # Initialize counts dataframe
    counts_dict = {'cell_id': adata.obs.index.values}

    # Count neighbors by cell type
    for celltype in celltypes:
        # Get mask of cells with this type
        type_mask = (adata.obs[celltype_col] == celltype).values

        # Count neighbors: sum connections to cells of this type
        # conn_matrix[i, j] = 1 if cell i is neighbor of cell j
        counts = np.array(conn_matrix[:, type_mask].sum(axis=1)).flatten()

        counts_dict[celltype] = counts

    counts_df = pd.DataFrame(counts_dict)

    print(f"   ✓ Neighbor counts computed: {counts_df.shape}")
    print(f"     Cell types: {celltypes}")
    print(f"     Average neighbors per cell: {counts_df[celltypes].sum(axis=1).mean():.1f}")

    return counts_df


def _normalize_rows(mat):
    row_sums = mat.sum(axis=1, keepdims=True) + 1e-12
    return mat / row_sums


def mean_max_cosine_similarity(components_a, components_b):
    """
    Compute mean max cosine similarity between topic components.
    Greedy matching is used for a fast stability proxy.
    """
    a = _normalize_rows(components_a)
    b = _normalize_rows(components_b)
    denom = (np.linalg.norm(a, axis=1, keepdims=True) + 1e-12) * (np.linalg.norm(b, axis=1, keepdims=True).T + 1e-12)
    sim = (a @ b.T) / denom

    used = set()
    scores = []
    for i in range(sim.shape[0]):
        best_j = None
        best_val = -1.0
        for j in range(sim.shape[1]):
            if j in used:
                continue
            if sim[i, j] > best_val:
                best_val = sim[i, j]
                best_j = j
        if best_j is not None:
            used.add(best_j)
            scores.append(best_val)

    return float(np.mean(scores)) if scores else 0.0


def run_k_diagnostics(X, output_dir):
    """
    Evaluate K choice using perplexity and topic stability across seeds.
    """
    print("\n4b. LDA K diagnostics (perplexity + stability)...")

    rng = np.random.RandomState(RANDOM_STATE)
    idx = rng.permutation(X.shape[0])
    split = int(0.8 * X.shape[0])
    X_train = X[idx[:split]]
    X_test = X[idx[split:]]

    perplexities = []
    stabilities = []

    for k in K_RANGE:
        lda = LatentDirichletAllocation(
            n_components=k,
            random_state=RANDOM_STATE,
            learning_method="online",
            max_iter=30,
            batch_size=1024,
            verbose=0
        )
        lda.fit(X_train)
        perplexities.append(lda.perplexity(X_test))

        # Stability across seeds
        comps = []
        for seed in STABILITY_SEEDS:
            lda_seed = LatentDirichletAllocation(
                n_components=k,
                random_state=seed,
                learning_method="online",
                max_iter=30,
                batch_size=1024,
                verbose=0
            )
            lda_seed.fit(X_train)
            comps.append(lda_seed.components_)

        pair_sims = []
        for i in range(len(comps)):
            for j in range(i + 1, len(comps)):
                pair_sims.append(mean_max_cosine_similarity(comps[i], comps[j]))
        stabilities.append(float(np.mean(pair_sims)) if pair_sims else 0.0)

        print(f"   K={k}: perplexity={perplexities[-1]:.2f}, stability={stabilities[-1]:.3f}")

    # Plot perplexity vs K
    fig, ax = plt.subplots(figsize=(7, 4), dpi=300)
    ax.plot(K_RANGE, perplexities, marker='o', linewidth=2, color='slateblue')
    ax.set_xlabel("K (topics)")
    ax.set_ylabel("Perplexity (lower is better)")
    ax.set_title("LDA K Selection - Perplexity")
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "05_lda_k_perplexity.png"), dpi=300, bbox_inches="tight")
    plt.close()

    # Plot stability vs K
    fig, ax = plt.subplots(figsize=(7, 4), dpi=300)
    ax.plot(K_RANGE, stabilities, marker='o', linewidth=2, color='teal')
    ax.set_xlabel("K (topics)")
    ax.set_ylabel("Topic stability (mean cosine similarity)")
    ax.set_title("LDA K Selection - Stability")
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "06_lda_k_stability.png"), dpi=300, bbox_inches="tight")
    plt.close()


# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    print("=" * 80)
    print("STEP 6: Spatial LDA - Discovering Neighborhood Topics")
    print("=" * 80)
    print(f"\nParameters:")
    print(f"  N_COMPONENTS (topics): {N_COMPONENTS}")
    print(f"  MAX_CELLS (subsample): {MAX_CELLS}")
    print(f"  RADIUS: {RADIUS}μm")
    print(f"  RANDOM_STATE: {RANDOM_STATE}\n")
    print(f"  K_RANGE (diagnostics): {K_RANGE}")
    print(f"  STABILITY_SEEDS: {STABILITY_SEEDS}\n")

    t0 = time.time()

    # ========================================================================
    # 1. Load h5ad with spatial graphs from Step 3
    # ========================================================================
    print("1. Loading h5ad with spatial graphs...")
    adata = sc.read_h5ad(INPUT_H5AD)

    print(f"   Total cells: {adata.n_obs:,}")
    print(f"   Unique patients: {adata.obs['Patient'].nunique()}")
    print(f"   Available graphs: {[k for k in adata.obsp.keys() if 'spatial' in k]}")

    # Filter to ER+ IMC only
    adata = adata[(adata.obs['Platform'] == 'IMC') & (adata.obs['subtype'] == 'ER+')].copy()
    print(f"   Filtered to ER+ IMC: {adata.n_obs:,} cells")

    # ========================================================================
    # 2. Compute neighbor counts from spatial graph
    # ========================================================================
    print(f"\n2. Computing neighbor counts from graph (r={RADIUS}μm)...")

    # Use COARSE cell types (5 categories)
    celltype_col = 'leidencelltype5'
    celltypes = ['endothelial', 'epithelial', 'fibroblast', 'immune', 'stromal']

    conn_key = f'spatial_connectivities_r{RADIUS}'

    counts_df = compute_neighbor_counts_from_graph(
        adata,
        conn_key=conn_key,
        celltype_col=celltype_col,
        celltypes=celltypes
    )

    # Save neighbor counts per cell (intermediate output)
    percell_csv = f"{OUTPUT_DIR}/ERpos_neighbor_counts_r{RADIUS}_coarse_percell.csv"
    counts_df.to_csv(percell_csv, index=False)
    print(f"   ✓ Saved neighbor counts: {percell_csv}")

    # ========================================================================
    # 3. Merge with patient metadata
    # ========================================================================
    print("\n3. Merging with patient metadata...")

    # Add patient info to counts
    patient_info = adata.obs[['Patient']].reset_index()
    patient_info.columns = ['cell_id', 'Patient']

    merged = counts_df.merge(patient_info, on='cell_id', how='inner')
    print(f"   Merged cells: {len(merged):,}")

    # Subsample if needed
    if len(merged) > MAX_CELLS:
        print(f"   Subsampling to {MAX_CELLS:,} cells (random_state={RANDOM_STATE})...")
        merged = merged.sample(n=MAX_CELLS, random_state=RANDOM_STATE)

    # ========================================================================
    # 4. Run LDA
    # ========================================================================
    print(f"\n4. Running LDA (n_components={N_COMPONENTS})...")
    X = merged[celltypes].values

    lda = LatentDirichletAllocation(
        n_components=N_COMPONENTS,
        random_state=RANDOM_STATE,
        learning_method="online",
        max_iter=30,
        batch_size=1024,
        verbose=0
    )

    doc_topics = lda.fit_transform(X)
    print(f"   LDA fit complete. Shape: {doc_topics.shape}")

    run_k_diagnostics(X, FIGURE_DIR)

    # ========================================================================
    # 5. Save topic assignments per cell
    # ========================================================================
    print("\n5. Aggregating topics by patient...")
    topic_cols = [f"topic_{i}" for i in range(N_COMPONENTS)]

    topic_df = pd.DataFrame(doc_topics, columns=topic_cols)
    topic_df['Patient'] = merged['Patient'].values

    # Aggregate by patient (mean topic weights)
    topic_by_patient = topic_df.groupby('Patient')[topic_cols].mean()

    # Save
    out_path = f"{OUTPUT_DIR}/ERpos_lda_topics_by_patient.csv"
    topic_by_patient.to_csv(out_path)
    print(f"   ✓ Saved: {out_path}")
    print(f"     Shape: {topic_by_patient.shape} (patients × topics)")

    # ========================================================================
    # 6. Save topic-celltype components
    # ========================================================================
    print("\n6. Extracting topic-celltype compositions...")
    components_df = pd.DataFrame(
        lda.components_,
        columns=celltypes,
        index=topic_cols
    )

    out_path = f"{OUTPUT_DIR}/ERpos_lda_topic_components.csv"
    components_df.to_csv(out_path)
    print(f"   ✓ Saved: {out_path}")
    print(f"     Shape: {components_df.shape} (topics × celltypes)")

    # Print dominant cell type per topic
    print("\n   Topic Interpretations (dominant cell type):")
    for topic in topic_cols:
        dominant = components_df.loc[topic].idxmax()
        weight = components_df.loc[topic].max()
        print(f"     {topic}: {dominant}-rich (weight={weight:.1f})")

    # ========================================================================
    # 7. Generate Visualizations
    # ========================================================================
    print("\n7. Generating visualizations...")

    # ------------------------------------------------------------------------
    # Plot 1: Topic Composition Heatmap (topics × cell types)
    # ------------------------------------------------------------------------
    print("   Generating Plot 1: Topic composition heatmap...")
    plt.figure(figsize=(8, 6), dpi=200)

    # Normalize by row for better visualization
    comp_norm = components_df.div(components_df.sum(axis=1), axis=0)

    sns.heatmap(
        comp_norm,
        cmap='YlOrRd',
        annot=True,
        fmt='.2f',
        cbar_kws={'label': 'Proportion', 'shrink': 0.5},
        linewidths=0.5
    )
    plt.title(f'LDA Topic Composition (n={len(merged):,} cells, k={N_COMPONENTS})')
    plt.xlabel('Cell Type')
    plt.ylabel('Topic')
    plt.tight_layout()
    plt.savefig(f"{FIGURE_DIR}/01_topic_composition_heatmap.png", bbox_inches='tight')
    plt.close()
    print(f"     ✓ Saved: {FIGURE_DIR}/01_topic_composition_heatmap.png")

    # ------------------------------------------------------------------------
    # Plot 2: Topic Weights by Patient (heatmap)
    # ------------------------------------------------------------------------
    print("   Generating Plot 2: Topic weights by patient...")
    plt.figure(figsize=(10, 12), dpi=200)

    sns.heatmap(
        topic_by_patient.T,  # Transpose: topics as rows, patients as cols
        cmap='viridis',
        cbar_kws={'label': 'Topic Weight', 'shrink': 0.5},
        yticklabels=topic_cols,
        xticklabels=False  # Too many patients
    )
    plt.title(f'Topic Distribution Across Patients (n={topic_by_patient.shape[0]} patients)')
    plt.xlabel(f'Patients (n={topic_by_patient.shape[0]})')
    plt.ylabel('Topic')
    plt.tight_layout()
    plt.savefig(f"{FIGURE_DIR}/02_topic_weights_by_patient.png", bbox_inches='tight')
    plt.close()
    print(f"     ✓ Saved: {FIGURE_DIR}/02_topic_weights_by_patient.png")

    # ------------------------------------------------------------------------
    # Plot 3: Mean Topic Weights (barplot)
    # ------------------------------------------------------------------------
    print("   Generating Plot 3: Mean topic weights...")
    plt.figure(figsize=(7, 4), dpi=200)

    topic_means = topic_by_patient.mean(axis=0).sort_values(ascending=False)

    sns.barplot(
        x=topic_means.index,
        y=topic_means.values,
        palette='coolwarm'
    )
    plt.title(f'Average Topic Prevalence (n={topic_by_patient.shape[0]} patients)')
    plt.xlabel('Topic')
    plt.ylabel('Mean Weight')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{FIGURE_DIR}/03_topic_mean_barplot.png", bbox_inches='tight')
    plt.close()
    print(f"     ✓ Saved: {FIGURE_DIR}/03_topic_mean_barplot.png")

    # ------------------------------------------------------------------------
    # Plot 4: Topic Distribution (violin plot)
    # ------------------------------------------------------------------------
    print("   Generating Plot 4: Topic distribution violin plot...")
    plt.figure(figsize=(10, 5), dpi=200)

    # Melt for violin plot
    topic_melted = topic_by_patient.reset_index().melt(
        id_vars='Patient',
        var_name='Topic',
        value_name='Weight'
    )

    sns.violinplot(
        data=topic_melted,
        x='Topic',
        y='Weight',
        palette='Set2',
        inner='box'
    )
    plt.title(f'Topic Weight Distribution Across Patients (n={topic_by_patient.shape[0]})')
    plt.xlabel('Topic')
    plt.ylabel('Weight')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{FIGURE_DIR}/04_topic_distribution_violin.png", bbox_inches='tight')
    plt.close()
    print(f"     ✓ Saved: {FIGURE_DIR}/04_topic_distribution_violin.png")

    # ========================================================================
    # SUMMARY
    # ========================================================================
    elapsed = time.time() - t0
    print("\n" + "=" * 80)
    print("STEP 6 COMPLETE")
    print("=" * 80)
    print(f"Execution time: {elapsed:.1f}s")
    print(f"\nOutputs:")
    print(f"  CSVs:")
    print(f"    - ERpos_neighbor_counts_r{RADIUS}_coarse_percell.csv (neighbor counts per cell)")
    print(f"    - ERpos_lda_topics_by_patient.csv ({topic_by_patient.shape[0]} patients × {N_COMPONENTS} topics)")
    print(f"    - ERpos_lda_topic_components.csv ({N_COMPONENTS} topics × {len(celltypes)} celltypes)")
    print(f"  Figures:")
    print(f"    - 01_topic_composition_heatmap.png")
    print(f"    - 02_topic_weights_by_patient.png")
    print(f"    - 03_topic_mean_barplot.png")
    print(f"    - 04_topic_distribution_violin.png")
    print(f"    - 05_lda_k_perplexity.png")
    print(f"    - 06_lda_k_stability.png")
    print("\nNext Step: Step 7 (Dispersion Metrics)")
    print("=" * 80)


if __name__ == "__main__":
    main()
