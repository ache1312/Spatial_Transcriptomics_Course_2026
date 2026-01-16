#!/usr/bin/env python3
"""
STEP 0: Dataset Description & Structure Overview

Purpose:
    Educational script for workshop students to understand the IMC dataset structure.

What it does:
    1. Loads the pre-processed h5ad file with embeddings and clinical data
    2. Prints detailed description of dataset components
    3. Generates markdown documentation (STEP0_DATASET_GUIDE.md)

Input:
    - data/imc_ERpos_with_embeddings_and_full_names.h5ad

Output:
    - Terminal: Detailed dataset description
    - docs/STEP0_DATASET_GUIDE.md: Complete documentation in English

"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime


def print_section(title, char="="):
    """Print a formatted section header"""
    print(f"\n{char * 80}")
    print(f"{title}")
    print(f"{char * 80}")


def describe_basic_structure(adata):
    """Describe basic dataset dimensions"""
    print_section("1. BASIC DATASET STRUCTURE", "=")

    print(f"\n📊 Dataset Dimensions:")
    print(f"   • Total cells: {adata.n_obs:,}")
    print(f"   • Total markers: {adata.n_vars}")
    print(f"   • Patients: {adata.obs['Patient'].nunique()}")
    print(f"   • Memory size: ~{adata.X.nbytes / 1024 / 1024:.1f} MB")

    print(f"\n🔬 Data Matrix (X):")
    print(f"   • Shape: {adata.X.shape}")
    print(f"   • Type: {type(adata.X).__name__}")
    print(f"   • Dtype: {adata.X.dtype}")

    # Get expression matrix for stats
    X_dense = adata.X if not hasattr(adata.X, 'toarray') else adata.X.toarray()
    non_zero_pct = 100 * np.count_nonzero(X_dense) / X_dense.size

    print(f"   • Non-zero values: {non_zero_pct:.1f}%")
    print(f"   • Value range: [{np.nanmin(X_dense):.2f}, {np.nanmax(X_dense):.2f}]")


def describe_markers(adata):
    """Describe the protein markers"""
    print_section("2. PROTEIN MARKERS (var)", "=")

    print(f"\n🧬 Available Markers ({adata.n_vars} total):")

    # Group markers by biological function
    marker_groups = {
        'Epithelial': ['panCK', 'Ecad', 'CK7', 'CK8', 'CK19', 'CK14', 'CK5', 'ER', 'PgR', 'HER2'],
        'Immune': ['CD45', 'CD3', 'CD20', 'CD68'],
        'Stromal': ['Vim', 'SMA', 'FN'],
        'Endothelial': ['CD31'],
        'Proliferation': ['Ki67'],
        'Functional': ['EGFR', 'CD44', 'CK7', 'CK5'],
        'Nuclear': ['DAPI_X', 'DAPI_Y'],
    }

    found_markers = {}
    other_markers = []

    for group, markers in marker_groups.items():
        found = [m for m in markers if m in adata.var_names]
        if found:
            found_markers[group] = found

    # Find markers not in predefined groups
    all_grouped = set([m for markers in marker_groups.values() for m in markers])
    other_markers = [m for m in adata.var_names if m not in all_grouped]

    # Print grouped markers
    for group, markers in found_markers.items():
        print(f"\n   {group} markers ({len(markers)}):")
        print(f"      {', '.join(markers)}")

    if other_markers:
        print(f"\n   Other markers ({len(other_markers)}):")
        # Print in rows of 8
        for i in range(0, len(other_markers), 8):
            print(f"      {', '.join(other_markers[i:i+8])}")


def describe_cell_metadata(adata):
    """Describe cell-level metadata (obs)"""
    print_section("3. CELL METADATA (obs)", "=")

    print(f"\n📋 Metadata Columns ({len(adata.obs.columns)} total):")

    # Categorize columns
    categories = {
        'Identifiers': ['Patient', 'slide_scene', 'Platform'],
        'Cell Types': ['leiden', 'celltype', 'leidencelltype5', 'celltype_full'],
        'Spatial': ['X_centroid', 'Y_centroid', 'area'],
        'Clinical': ['subtype', 'Survival', 'Survival_time', 'Recurrence', 'Recurrence_time'],
        'Cohort': ['Cohort', 'cohort']
    }

    for category, cols in categories.items():
        found_cols = [c for c in cols if c in adata.obs.columns]
        if found_cols:
            print(f"\n   {category}:")
            for col in found_cols:
                dtype = adata.obs[col].dtype
                n_unique = adata.obs[col].nunique() if dtype == 'object' or dtype.name == 'category' else 'N/A'

                if n_unique != 'N/A' and n_unique < 20:
                    print(f"      • {col:25} ({dtype}, {n_unique} unique values)")
                elif n_unique != 'N/A':
                    print(f"      • {col:25} ({dtype}, {n_unique} unique values)")
                else:
                    range_info = f"[{adata.obs[col].min():.1f}, {adata.obs[col].max():.1f}]"
                    print(f"      • {col:25} ({dtype}, range: {range_info})")

    # Print other columns
    all_categorized = set([c for cols in categories.values() for c in cols])
    other_cols = [c for c in adata.obs.columns if c not in all_categorized]

    if other_cols:
        print(f"\n   Other columns ({len(other_cols)}):")
        for col in other_cols[:10]:  # Show first 10
            print(f"      • {col}")
        if len(other_cols) > 10:
            print(f"      ... and {len(other_cols) - 10} more")


def describe_embeddings(adata):
    """Describe dimensionality reduction embeddings"""
    print_section("4. DIMENSIONALITY REDUCTION (obsm)", "=")

    if 'X_pca' in adata.obsm:
        print(f"\n🔵 PCA (Principal Component Analysis):")
        print(f"   • Embedding shape: {adata.obsm['X_pca'].shape}")
        print(f"   • Components: {adata.obsm['X_pca'].shape[1]}")

        if 'pca' in adata.uns:
            variance_ratio = adata.uns['pca']['variance_ratio']
            cumsum = np.cumsum(variance_ratio)

            # Find components for different variance thresholds
            pc_90 = np.argmax(cumsum >= 0.90) + 1
            pc_95 = np.argmax(cumsum >= 0.95) + 1
            pc_99 = np.argmax(cumsum >= 0.99) + 1

            print(f"   • Variance explained:")
            print(f"      - First 5 PCs: {cumsum[4] * 100:.1f}%")
            print(f"      - First 10 PCs: {cumsum[9] * 100:.1f}%")
            print(f"      - First 20 PCs: {cumsum[19] * 100:.1f}%")
            print(f"   • Components needed for:")
            print(f"      - 90% variance: {pc_90} PCs")
            print(f"      - 95% variance: {pc_95} PCs")
            print(f"      - 99% variance: {pc_99} PCs")
    else:
        print(f"\n⚠️  PCA not found in obsm")

    if 'X_umap' in adata.obsm:
        print(f"\n🔴 UMAP (Uniform Manifold Approximation):")
        print(f"   • Embedding shape: {adata.obsm['X_umap'].shape}")
        print(f"   • Dimensions: {adata.obsm['X_umap'].shape[1]}D")
        print(f"   • Range X: [{adata.obsm['X_umap'][:, 0].min():.2f}, {adata.obsm['X_umap'][:, 0].max():.2f}]")
        print(f"   • Range Y: [{adata.obsm['X_umap'][:, 1].min():.2f}, {adata.obsm['X_umap'][:, 1].max():.2f}]")
    else:
        print(f"\n⚠️  UMAP not found in obsm")


def describe_cell_types(adata):
    """Describe cell type annotations"""
    print_section("5. CELL TYPE ANNOTATIONS", "=")

    if 'celltype_full' in adata.obs.columns:
        print(f"\n🏷️  Full Cell Type Names (celltype_full):")
        celltype_counts = adata.obs['celltype_full'].value_counts()

        for celltype, count in celltype_counts.items():
            pct = 100 * count / len(adata)
            print(f"   • {celltype:40} {count:>8,} cells ({pct:>5.2f}%)")

    if 'leidencelltype5' in adata.obs.columns:
        print(f"\n🔢 Leiden Clusters (leidencelltype5):")
        leiden_counts = adata.obs['leidencelltype5'].value_counts()
        print(f"   • Total clusters: {leiden_counts.shape[0]}")
        print(f"   • Largest cluster: {leiden_counts.iloc[0]:,} cells ({100 * leiden_counts.iloc[0] / len(adata):.1f}%)")
        print(f"   • Smallest cluster: {leiden_counts.iloc[-1]:,} cells ({100 * leiden_counts.iloc[-1] / len(adata):.1f}%)")


def describe_cohorts(adata):
    """Describe cohort distribution"""
    print_section("6. COHORT DISTRIBUTION", "=")

    if 'Cohort' in adata.obs.columns or 'cohort' in adata.obs.columns:
        cohort_col = 'Cohort' if 'Cohort' in adata.obs.columns else 'cohort'

        print(f"\n🏥 Cohort Breakdown:")
        cohort_summary = adata.obs.groupby(cohort_col).agg({
            'Patient': 'nunique'
        })
        cohort_summary['n_cells'] = adata.obs.groupby(cohort_col).size()

        for cohort_name, row in cohort_summary.iterrows():
            pct = 100 * row['n_cells'] / len(adata)
            print(f"   • {cohort_name:15} {row['n_cells']:>8,} cells ({pct:>5.1f}%) from {row['Patient']:>3} patients")


def describe_clinical_data(adata):
    """Describe clinical annotations"""
    print_section("7. CLINICAL DATA", "=")

    survival_cols = ['Survival', 'Survival_time', 'Recurrence', 'Recurrence_time']

    print(f"\n🏥 Survival Data Availability:")
    for col in survival_cols:
        if col in adata.obs.columns:
            n_patients_with_data = adata.obs.groupby('Patient')[col].first().notna().sum()
            total_patients = adata.obs['Patient'].nunique()

            if 'time' in col.lower():
                min_time = adata.obs[col].min()
                max_time = adata.obs[col].max()
                median_time = adata.obs[col].median()
                print(f"   • {col:25} Available for {n_patients_with_data}/{total_patients} patients")
                print(f"     {'':27} Range: {min_time:.1f} - {max_time:.1f} months (median: {median_time:.1f})")
            else:
                n_events = adata.obs.groupby('Patient')[col].first().sum()
                print(f"   • {col:25} Available for {n_patients_with_data}/{total_patients} patients ({n_events:.0f} events)")

    if 'subtype' in adata.obs.columns:
        print(f"\n💊 Cancer Subtype:")
        subtypes = adata.obs['subtype'].unique()
        print(f"   • Subtypes in dataset: {', '.join([str(s) for s in subtypes if pd.notna(s)])}")


def describe_unstructured(adata):
    """Describe unstructured annotations"""
    print_section("8. UNSTRUCTURED METADATA (uns)", "=")

    print(f"\n📦 Available Metadata:")
    for key in adata.uns.keys():
        value_type = type(adata.uns[key]).__name__

        if isinstance(adata.uns[key], dict):
            n_items = len(adata.uns[key])
            print(f"   • {key:30} dict with {n_items} items")
        elif isinstance(adata.uns[key], (list, np.ndarray)):
            length = len(adata.uns[key])
            print(f"   • {key:30} {value_type} with {length} elements")
        else:
            print(f"   • {key:30} {value_type}")


def generate_markdown_guide(adata, output_path):
    """Generate comprehensive markdown documentation"""

    md_content = f"""# STEP 0: Dataset Structure & Overview

**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
**Dataset**: ER+ Breast Cancer IMC Data (Eng et al. 2019)
**File**: `data/imc_ERpos_with_embeddings_and_full_names.h5ad`

---

## Overview

This dataset contains **{adata.n_obs:,} cells** from **{adata.obs['Patient'].nunique()} patients** with ER+ breast cancer, profiled using **Imaging Mass Cytometry (IMC)** technology. The data has been pre-processed to include:

- ✅ Quality-filtered single cells
- ✅ Protein marker intensities ({adata.n_vars} markers)
- ✅ Cell type annotations (Leiden clustering)
- ✅ Dimensionality reduction (PCA, UMAP)
- ✅ Clinical survival data
- ✅ Cohort labels (Basel vs Zürich)

---

## 1. Basic Structure

### AnnData Object Components

```python
import scanpy as sc
adata = sc.read_h5ad('data/imc_ERpos_with_embeddings_and_full_names.h5ad')
```

| Component | Description | Shape/Count |
|-----------|-------------|-------------|
| **X** | Expression matrix (marker intensities) | {adata.X.shape[0]:,} cells × {adata.X.shape[1]} markers |
| **obs** | Cell-level metadata | {len(adata.obs.columns)} columns |
| **var** | Marker-level metadata | {adata.n_vars} markers |
| **obsm** | Embeddings (PCA, UMAP) | {len(adata.obsm.keys())} embedding(s) |
| **uns** | Unstructured metadata | {len(adata.uns.keys())} item(s) |

---

## 2. Protein Markers (var)

**Total markers**: {adata.n_vars}

"""

    # Add marker groups
    marker_groups = {
        'Epithelial': ['panCK', 'Ecad', 'CK7', 'CK8', 'CK19', 'CK14', 'CK5', 'ER', 'PgR', 'HER2'],
        'Immune': ['CD45', 'CD3', 'CD20', 'CD68'],
        'Stromal': ['Vim', 'SMA', 'FN'],
        'Endothelial': ['CD31'],
        'Proliferation': ['Ki67'],
        'Functional': ['EGFR', 'CD44'],
        'Nuclear': ['DAPI_X', 'DAPI_Y'],
    }

    md_content += "\n📌 Full marker list and target descriptions: `docs/MARKERS_PRESENT.md`\n"
    md_content += "\n### Complete marker list\n"
    md_content += f"```\n{', '.join(map(str, adata.var_names))}\n```\n"

    for group, markers in marker_groups.items():
        found = [m for m in markers if m in adata.var_names]
        if found:
            md_content += f"\n### {group} Markers\n"
            md_content += f"```\n{', '.join(found)}\n```\n"

    # Cell types
    if 'celltype_full' in adata.obs.columns:
        md_content += "\n---\n\n## 3. Cell Type Annotations\n\n"
        md_content += "### Full Cell Type Names\n\n"
        md_content += "| Cell Type | Count | Percentage |\n"
        md_content += "|-----------|------:|-----------:|\n"

        celltype_counts = adata.obs['celltype_full'].value_counts()
        for celltype, count in celltype_counts.items():
            pct = 100 * count / len(adata)
            md_content += f"| {celltype} | {count:,} | {pct:.2f}% |\n"

    # Embeddings
    md_content += "\n---\n\n## 4. Dimensionality Reduction\n\n"

    if 'X_pca' in adata.obsm:
        md_content += "### PCA (Principal Component Analysis)\n\n"
        md_content += f"- **Components**: {adata.obsm['X_pca'].shape[1]}\n"

        if 'pca' in adata.uns:
            variance_ratio = adata.uns['pca']['variance_ratio']
            cumsum = np.cumsum(variance_ratio)
            pc_95 = np.argmax(cumsum >= 0.95) + 1

            md_content += f"- **Variance explained by first 20 PCs**: {cumsum[19] * 100:.1f}%\n"
            md_content += f"- **PCs needed for 95% variance**: {pc_95}\n\n"

    if 'X_umap' in adata.obsm:
        md_content += "### UMAP (Uniform Manifold Approximation)\n\n"
        md_content += f"- **Dimensions**: {adata.obsm['X_umap'].shape[1]}D\n"
        md_content += f"- **Used for**: Visualization and exploratory analysis\n\n"

    # Cohorts
    cohort_col = None
    if 'Cohort' in adata.obs.columns:
        cohort_col = 'Cohort'
    elif 'cohort' in adata.obs.columns:
        cohort_col = 'cohort'

    if cohort_col is not None:
        md_content += "---\n\n## 5. Cohort Distribution\n\n"

        cohort_summary = adata.obs.groupby(cohort_col).agg({
            'Patient': 'nunique'
        })
        cohort_summary['n_cells'] = adata.obs.groupby(cohort_col).size()

        md_content += "| Cohort | Patients | Cells | Percentage |\n"
        md_content += "|--------|----------|------:|-----------:|\n"

        for cohort_name, row in cohort_summary.iterrows():
            pct = 100 * row['n_cells'] / len(adata)
            md_content += f"| {cohort_name} | {row['Patient']} | {row['n_cells']:,} | {pct:.1f}% |\n"

    # Clinical data
    md_content += "\n---\n\n## 6. Clinical Data\n\n"
    md_content += "### Survival Endpoints\n\n"

    survival_cols = ['Survival', 'Survival_time', 'Recurrence', 'Recurrence_time']
    available_survival = [col for col in survival_cols if col in adata.obs.columns]

    if available_survival:
        md_content += "| Endpoint | Availability |\n"
        md_content += "|----------|-------------|\n"

        for col in available_survival:
            n_patients = adata.obs.groupby('Patient')[col].first().notna().sum()
            total_patients = adata.obs['Patient'].nunique()
            md_content += f"| {col} | {n_patients}/{total_patients} patients |\n"

    # Usage examples
    md_content += """
---

## 7. Quick Start Examples

### Load the dataset
```python
import scanpy as sc

# Load h5ad file
adata = sc.read_h5ad('data/imc_ERpos_with_embeddings_and_full_names.h5ad')

# Basic info
print(f"Cells: {adata.n_obs:,}")
print(f"Markers: {adata.n_vars}")
print(f"Patients: {adata.obs['Patient'].nunique()}")
```

### Access different components
```python
# Expression matrix
X = adata.X  # Marker intensities

# Cell metadata
cell_types = adata.obs['celltype_full']
patients = adata.obs['Patient']
survival = adata.obs['Survival_time']

# Embeddings
pca_coords = adata.obsm['X_pca']
umap_coords = adata.obsm['X_umap']

# Marker names
markers = adata.var_names
```

### Filter data
```python
# Get only immune cells
immune_cells = adata[adata.obs['celltype'] == 'immune']

# Get specific patient
patient_58 = adata[adata.obs['Patient'] == '58']

# Get Basel cohort
basel = adata[adata.obs['Cohort'] == 'Basel']
```

---

## 8. Next Steps in Pipeline

After understanding the dataset structure (Step 0), the pipeline continues:

1. **Step 1**: IMC Quality Control Visualizations
   - Visual inspection of 3 workshop patients
   - 6 biological views per patient

2. **Step 2**: Embeddings & QC Analysis
   - PCA variance analysis
   - UMAP visualization
   - Quality control plots

3. **Step 3+**: Spatial Biology Analyses
   - Neighbor graphs
   - Cell-cell interactions
   - Survival analysis
   - And more...

---

## References

- **Original Paper**: Eng et al. (2019) Nature Medicine
- **Technology**: Imaging Mass Cytometry (IMC)
- **Dataset**: ER+ Breast Cancer TMAs
- **Cohorts**: Basel University Hospital & University Hospital Zürich

"""

    # Write to file
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(md_content)

    print(f"\n✅ Markdown guide generated: {output_path}")


def main():
    """Main function to describe the dataset"""

    print_section("STEP 0: DATASET DESCRIPTION & STRUCTURE", "=")
    print("\n📚 Educational overview for workshop students")
    print("   Purpose: Understand what's in the h5ad file before analysis")

    # Load dataset
    print("\n🔄 Loading dataset...")
    adata_path = Path("data/imc_ERpos_with_embeddings_and_full_names.h5ad")

    if not adata_path.exists():
        print(f"\n❌ ERROR: Dataset not found at {adata_path}")
        print("   Please ensure the h5ad file exists before running Step 0")
        return 1

    adata = sc.read_h5ad(adata_path)
    print(f"✅ Loaded: {adata.n_obs:,} cells × {adata.n_vars} markers")

    # Describe all components
    describe_basic_structure(adata)
    describe_markers(adata)
    describe_cell_metadata(adata)
    describe_embeddings(adata)
    describe_cell_types(adata)
    describe_cohorts(adata)
    describe_clinical_data(adata)
    describe_unstructured(adata)

    # Generate markdown guide
    print_section("9. GENERATING DOCUMENTATION", "=")
    md_output = Path("docs/STEP0_DATASET_GUIDE.md")
    generate_markdown_guide(adata, md_output)

    # Summary
    print_section("STEP 0 COMPLETE", "=")
    print("\n✅ Dataset description complete!")
    print(f"\n📄 Documentation:")
    print(f"   • Terminal output above (detailed description)")
    print(f"   • {md_output} (markdown guide for students)")

    print(f"\n📊 Quick Summary:")
    print(f"   • Cells: {adata.n_obs:,}")
    print(f"   • Markers: {adata.n_vars}")
    print(f"   • Patients: {adata.obs['Patient'].nunique()}")
    print(f"   • Cell types: {adata.obs['celltype_full'].nunique()}")
    print(f"   • Has PCA: {'✓' if 'X_pca' in adata.obsm else '✗'}")
    print(f"   • Has UMAP: {'✓' if 'X_umap' in adata.obsm else '✗'}")
    print(f"   • Has survival data: {'✓' if 'Survival' in adata.obs.columns else '✗'}")

    print("\n🎯 Next step: Step 1 (IMC visualizations)")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    exit(main())
