#!/usr/bin/env python3
"""
Napari Visualization of REAL IMC Data with Segmentation
========================================================

Visualizes actual IMC OME-TIFF images with segmentation masks and cell type metadata.

Data Structure:
- IMC Images: OME-TIFF multicanal (23 channels, ~700×700 px)
- Segmentation Masks: Single-channel TIFF (uint8, cell IDs)
- Cell Type Metadata: AnnData obs (leiden, celltype_detailed)

Features:
1. Load real IMC 23-channel images
2. Display segmentation masks colored by cell type
3. Interactive channel toggling
4. Cell type legend and statistics
5. Export annotated snapshots

Paper: Eng et al. 2025, JCI Insight
Data: Jackson et al., Nature 2020 (IMC ER+ breast cancer)

"""

import numpy as np
import pandas as pd
import anndata as ad
import tifffile
from pathlib import Path
from typing import List, Tuple, Optional, Dict
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb, ListedColormap, LinearSegmentedColormap
import warnings
warnings.filterwarnings('ignore')

# Napari is optional - only needed for interactive mode
try:
    import napari
    NAPARI_AVAILABLE = True
except ImportError:
    NAPARI_AVAILABLE = False

# ============================================================================
# CONFIGURATION
# ============================================================================

# Get base directory (pipeline_leiden_validated/)
BASE_DIR = Path(__file__).parent.parent.resolve()

PATHS = {
    # Use unified h5ad with all metadata
    'adata': BASE_DIR / 'data/imc_ERpos_with_embeddings_and_full_names.h5ad',
    'imc_images': BASE_DIR / 'data/imc_images',
    'cell_masks': BASE_DIR / 'data/segmentation_masks',
    'output': BASE_DIR / 'figures/step1_imc_visualizations'
}

# IMC Panel (81 channels from Jackson et al. Nature 2020)
IMC_MARKER_NAMES = [
    'CD20', 'CD31', 'CD3', 'CD44', 'CD45', 'CD68', 'CAIX', 'CK19', 'CK5', 'CK7',
    'CK8', 'DAPI1', 'EGFR', 'ER', 'Ecad', 'EpCAM', 'Erk12', 'FN', 'GATA3', 'pHH3',
    'CK14', 'Ki67', 'PgR', 'S6', 'SMA', 'Slug', 'Sox9', 'Twist', 'Vim', 'bCatenin',
    'cMyc', 'HER2', 'CC3cPARP', 'mTOR', 'p53', 'panCK', 'Beta-catenin-n', 'CD11b', 'CD11c', 'CD138',
    'CD16', 'CD209', 'CD45RO', 'CD4', 'CD56', 'CD63', 'CD8', 'FoxP3', 'H3K27me3', 'H3K9ac',
    'HLA-Class-1', 'HLA-DR', 'IDO', 'CK17', 'CK6', 'Lag3', 'MPO', 'PD-L1', 'PD1', 'pS6',
    'BMP2', 'ColI', 'ColIV', 'CoxIV', 'Glut1', 'PDGFRa', 'PDPN', 'pAKT', 'AR', 'GRNZB',
    'H3K4', 'HIF1a', 'LamAC', 'LamB1', 'LamB2', 'PCNA', 'R1c2', 'cPARP', 'gH2AX', 'pERK',
    'pRB'
]

# Marker display names (technical → publication-ready)
MARKER_DISPLAY_NAMES = {
    'DAPI1': 'DAPI',  # Remove numbering for cleaner visualization
}

# Cell type colors (from paper)
CELLTYPE_COLORS = {
    'Luminal ER+ t.': '#FF6B6B',
    'Prolif. t.': '#E74C3C',
    'HER2+ ER+ t.': '#C0392B',
    'HER2+ t.': '#8E44AD',
    'Basal t.': '#E67E22',
    'Myoepithelial': '#F39C12',
    'CD3 T cell': '#3498DB',
    'CD20 B cell': '#2ECC71',
    'Macrophage': '#1ABC9C',
    'Quies. str.': '#95A5A6',
    'FN+ FB': '#E8DAEF',
    'Vim+ FB': '#D7BDE2',
    'Pericyte/SMA+ FB': '#BB8FCE',
    'endothelial': '#16A085',
}

# ============================================================================
# BIOLOGICAL COMPOSITE VIEWS - WORKSHOP PRESETS
# ============================================================================

BIOLOGICAL_VIEWS = {
    'tumor_architecture': {
        'name': 'Tumor Architecture',
        'markers': ['panCK', 'ER', 'Ki67', 'HER2'],
        'colors': ['cyan', 'red', 'green', 'magenta'],
        'indices': [35, 13, 21, 31],  # Will be verified at runtime
        'description': 'Tumor compartment identity and subtypes',
        'connects_to': 'Celltype assignment (step 1), Leiden clustering',
        'look_for': [
            'ER+ tumor cells (red + cyan overlap)',
            'Proliferating regions (green = Ki67+)',
            'HER2+ subpopulations (magenta)',
            'Tumor spatial organization and nests'
        ]
    },
    'immune_landscape': {
        'name': 'Immune Landscape',
        'markers': ['CD45', 'CD3', 'CD68'],  # CD20 removed (noisy)
        'colors': ['blue', 'yellow', 'red'],
        'indices': [4, 2, 5],
        'description': 'Complete immune compartment composition',
        'connects_to': 'Immune infiltration, mixing score (step 10), cell-cell interactions (step 6)',
        'look_for': [
            'T cells (yellow = CD3+)',
            'B cells (green = CD20+)',
            'Macrophages (red = CD68+)',
            'Immune excluded vs infiltrated regions'
        ]
    },
    'stromal_organization': {
        'name': 'Stromal Organization',
        'markers': ['SMA', 'FN', 'Vim', 'panCK'],
        'colors': ['red', 'green', 'blue', 'magenta'],
        'indices': [24, 17, 28, 35],
        'description': 'CAF subtypes and tumor-stroma interface',
        'connects_to': 'Dispersion metrics (step 7), spatial barriers, CAF analysis',
        'look_for': [
            'Activated CAFs (SMA+ = red)',
            'Fibronectin-rich stroma (green)',
            'Mesenchymal cells (Vim+ = blue)',
            'Tumor-stroma boundaries (magenta vs stromal colors)'
        ]
    },
    'invasive_edge': {
        'name': 'Invasive Edge',
        'markers': ['panCK', 'CD44', 'Vim', 'EGFR'],
        'colors': ['cyan', 'red', 'green', 'yellow'],
        'indices': [35, 3, 28, 12],
        'description': 'Tumor plasticity and invasive phenotypes',
        'connects_to': 'EMT programs, survival analysis (step 8), invasive behavior',
        'look_for': [
            'Tumor cells (cyan = panCK+)',
            'Stem-like/invasive tumor (red = CD44+)',
            'EMT-like phenotype (cyan + green = panCK+ Vim+)',
            'EGFR expression (yellow) on invasive edges'
        ]
    },
    'vascular_proximity': {
        'name': 'Vascular Proximity',
        'markers': ['CD31', 'panCK', 'CD44', 'DAPI1'],
        'colors': ['red', 'green', 'yellow', 'blue'],
        'indices': [1, 35, 3, 11],
        'description': 'Tumor-vascular interface and prognostic signals',
        'connects_to': 'Survival analysis (step 8), vascular proximity metrics, prognostic biomarkers',
        'look_for': [
            'Blood vessels (red = CD31+)',
            'Tumor cells (green = panCK+)',
            'Invasive tumor near vessels (yellow + green near red)',
            'Tumor-vascular proximity (prognostic)'
        ]
    },
    'proliferation_receptors': {
        'name': 'Proliferation & Receptors',
        'markers': ['Ki67', 'ER', 'PgR', 'HER2'],
        'colors': ['green', 'red', 'cyan', 'magenta'],
        'indices': [21, 13, 22, 31],
        'description': 'Clinical subtypes and proliferative index',
        'connects_to': 'Clinical stratification, survival analysis (step 8), receptor status',
        'look_for': [
            'Proliferating cells (green = Ki67+)',
            'ER+ cells (red)',
            'PgR+ cells (cyan)',
            'HER2+ cells (magenta)',
            'Hormone receptor positive regions (red + cyan)'
        ]
    }
}

# ============================================================================
# COLORMAP FUNCTIONS
# ============================================================================

def create_black_to_color_cmap(color_name: str):
    """
    Create custom colormap from black to specified color

    Parameters:
    -----------
    color_name : Color name ('blue', 'red', 'green', 'yellow', 'cyan', 'magenta')

    Returns:
    --------
    cmap : LinearSegmentedColormap from black (0,0,0) to color
    """
    color_map = {
        'blue': (0, 0, 1),
        'red': (1, 0, 0),
        'green': (0, 1, 0),
        'yellow': (1, 1, 0),
        'cyan': (0, 1, 1),
        'magenta': (1, 0, 1)
    }

    if color_name.lower() not in color_map:
        # Default to grayscale
        end_color = (1, 1, 1)
    else:
        end_color = color_map[color_name.lower()]

    # Create colormap from black (0,0,0) to end_color
    cmap = LinearSegmentedColormap.from_list(
        f'black_to_{color_name}',
        [(0, 0, 0), end_color]
    )

    return cmap


def get_display_name(marker_name: str) -> str:
    """
    Get publication-ready display name for a marker

    Parameters:
    -----------
    marker_name : Technical marker name from dataset

    Returns:
    --------
    display_name : Clean name for visualization (e.g., "DAPI1" → "DAPI")
    """
    return MARKER_DISPLAY_NAMES.get(marker_name, marker_name)


# ============================================================================
# DATA LOADING FUNCTIONS
# ============================================================================

def load_imc_anndata(path: Path) -> ad.AnnData:
    """Load AnnData with cell metadata"""
    print(f"\n{'='*70}")
    print("📂 Loading AnnData (Cell Metadata)")
    print(f"{'='*70}")

    adata = ad.read_h5ad(path)
    print(f"  ✓ Cells: {adata.n_obs:,}")
    print(f"  ✓ Markers: {adata.n_vars}")
    print(f"  ✓ Patients: {adata.obs['Patient'].nunique()}")

    # Check for cell type annotations
    celltype_col = 'celltype_detailed' if 'celltype_detailed' in adata.obs else 'leidencelltype5'
    print(f"  ✓ Cell type column: {celltype_col}")

    unique_types = adata.obs[celltype_col].nunique()
    print(f"  ✓ Unique cell types: {unique_types}")

    return adata


def find_matching_files(
    roi_id: str,
    imc_dir: Path,
    mask_dir: Path
) -> Tuple[Optional[Path], Optional[Path]]:
    """
    Find matching IMC image and mask files for a given ROI

    ROI ID format: "X5Y2_124" from slide_scene like "2_20170905_44_24_X5Y2_124_a0"
    """
    print(f"\n🔍 Searching for ROI: {roi_id}")

    # Search for image
    imc_files = list(imc_dir.glob(f"*{roi_id}*.tiff")) + list(imc_dir.glob(f"*{roi_id}*.tif"))
    mask_files = list(mask_dir.glob(f"*{roi_id}*.tiff")) + list(mask_dir.glob(f"*{roi_id}*.tif"))

    if not imc_files:
        print(f"  ⚠️  No IMC image found for {roi_id}")
        return None, None

    if not mask_files:
        print(f"  ⚠️  No mask found for {roi_id}")
        return imc_files[0], None

    print(f"  ✓ IMC image: {imc_files[0].name}")
    print(f"  ✓ Mask: {mask_files[0].name}")

    return imc_files[0], mask_files[0]


def load_imc_image(path: Path) -> np.ndarray:
    """
    Load IMC OME-TIFF multichannel image

    Returns: (C, H, W) array with C channels
    """
    print(f"\n📷 Loading IMC image: {path.name}")

    img = tifffile.imread(path)

    print(f"  ✓ Shape: {img.shape} (Channels, Height, Width)")
    print(f"  ✓ Dtype: {img.dtype}")
    print(f"  ✓ Channels: {img.shape[0]}")
    print(f"  ✓ Range: [{img.min():.2f}, {img.max():.2f}]")

    return img


def load_segmentation_mask(path: Path) -> np.ndarray:
    """
    Load segmentation mask (cell IDs)

    Returns: (H, W) array with cell IDs
    """
    print(f"\n🎭 Loading segmentation mask: {path.name}")

    mask = tifffile.imread(path)

    n_cells = len(np.unique(mask)) - 1  # -1 for background (0)

    print(f"  ✓ Shape: {mask.shape}")
    print(f"  ✓ Dtype: {mask.dtype}")
    print(f"  ✓ Cells: {n_cells}")
    print(f"  ✓ Cell ID range: {mask.min()}-{mask.max()}")

    return mask


def create_cell_contours(mask: np.ndarray) -> np.ndarray:
    """
    Create cell boundary/contour overlay from segmentation mask

    Parameters:
    -----------
    mask : np.ndarray
        Segmentation mask (H, W) with cell IDs

    Returns:
    --------
    contours : (H, W) boolean array with cell boundaries
    """
    from scipy.ndimage import binary_erosion

    # Create binary mask (cells vs background)
    binary_mask = mask > 0

    # Find boundaries using erosion
    eroded = binary_erosion(binary_mask)
    boundaries = binary_mask & ~eroded

    return boundaries


def create_celltype_colored_mask(
    mask: np.ndarray,
    adata: ad.AnnData,
    roi_id: str,
    slide_scene: str = None,
    celltype_col: str = 'celltype_detailed',
    contours_only: bool = False
) -> Tuple[np.ndarray, Dict[int, str], Dict[str, np.ndarray]]:
    """
    Create visualization of segmented cells

    NEW: Now defaults to showing only cell contours (outlines) instead of
    filled cell types, which shows ALL cells in the mask regardless of QC filtering.

    Parameters:
    -----------
    mask : np.ndarray
        Segmentation mask (H, W) with cell IDs
    adata : ad.AnnData
        AnnData with cell metadata (not used if contours_only=True)
    roi_id : str
        ROI identifier (e.g., "X3Y1")
    slide_scene : str, optional
        Full slide_scene value from adata (e.g., "257_X3Y1-28")
    contours_only : bool
        If True, show only cell boundaries (default). If False, use cell type colors.

    Returns:
    --------
    colored_mask : (H, W, 3) RGB array
    cell_id_to_type : dict mapping cell ID → cell type (empty if contours_only)
    celltype_colors_used : dict mapping cell type → RGB color (empty if contours_only)
    """
    if contours_only:
        # Simple contour-based visualization (no cell types needed)
        print(f"\n🎨 Creating cell contour overlay (showing ALL {len(np.unique(mask))-1} cells)")

        # Get cell boundaries
        contours = create_cell_contours(mask)

        # Create colored contour overlay (green contours on transparent background)
        h, w = mask.shape
        colored_mask = np.zeros((h, w, 3), dtype=np.float32)
        colored_mask[contours] = [0, 1, 0]  # Bright green contours

        print(f"  ✓ Created contours for all cells")
        print(f"  ✓ Contour pixels: {np.sum(contours):,}")

        return colored_mask, {}, {}

    # Original cell type coloring logic (if contours_only=False)
    print(f"\n🎨 Creating cell type colored mask")

    # Filter adata to this ROI
    if slide_scene:
        roi_cells = adata.obs[adata.obs['slide_scene'] == slide_scene]
        print(f"  🔍 Filtering by slide_scene: {slide_scene}")
    else:
        roi_cells = adata.obs[adata.obs['slide_scene'].str.contains(roi_id, na=False)]
        print(f"  🔍 Filtering by ROI pattern: {roi_id}")

    if len(roi_cells) == 0:
        print(f"  ⚠️  No cells found in adata for ROI {roi_id}")
        # Return grayscale mask
        colored = np.stack([mask.astype(float)/mask.max()]*3, axis=-1)
        return colored, {}, {}

    print(f"  ✓ Found {len(roi_cells)} cells in adata for this ROI")

    # Get cell types
    if celltype_col not in roi_cells.columns:
        celltype_col = 'leidencelltype5'

    celltypes = roi_cells[celltype_col].values

    # Create mapping: cell_id → cell_type
    # Use spatial coordinates to match cells between mask and adata
    unique_mask_ids = np.unique(mask)[1:]  # Exclude 0 (background)
    total_cells_in_mask = len(unique_mask_ids)

    # Map cells using spatial coordinates if available
    cell_id_to_type = {}
    if 'DAPI_X' in roi_cells.columns and 'DAPI_Y' in roi_cells.columns:
        # For each cell in adata, find its position and estimate cell ID
        # Since cell IDs are sequential, we can approximate based on position
        roi_cells_sorted = roi_cells.sort_index()

        # Simple approach: assume cells in adata are a QC-filtered subset of mask cells
        # Map the first len(celltypes) cells by assuming sequential processing
        for i, celltype in enumerate(celltypes):
            if i < total_cells_in_mask:
                cell_id = unique_mask_ids[i]
                cell_id_to_type[cell_id] = celltype
    else:
        # Fallback: sequential mapping
        for i, cell_id in enumerate(unique_mask_ids):
            if i < len(celltypes):
                cell_id_to_type[cell_id] = celltypes[i]

    # Create colored mask
    h, w = mask.shape
    colored_mask = np.zeros((h, w, 3), dtype=np.float32)

    # Color each cell randomly 
    print(f"  🎨 Assigning random colors to all {len(unique_mask_ids)} cells...")

    # Use a colormap for diverse, distinct colors
    from matplotlib import cm
    import matplotlib.pyplot as plt

    # Generate distinct colors using a categorical colormap
    # Use tab20 + tab20b + tab20c for up to 60 distinct colors
    colormap = plt.cm.get_cmap('tab20')
    colormap_b = plt.cm.get_cmap('tab20b')
    colormap_c = plt.cm.get_cmap('tab20c')

    for i, cell_id in enumerate(unique_mask_ids):
        # Cycle through colormaps for distinct colors
        if i < 20:
            color_rgb = colormap(i % 20)[:3]
        elif i < 40:
            color_rgb = colormap_b((i - 20) % 20)[:3]
        elif i < 60:
            color_rgb = colormap_c((i - 40) % 20)[:3]
        else:
            # For cells beyond 60, generate random bright colors
            np.random.seed(cell_id)  # Consistent color per cell ID
            color_rgb = np.random.rand(3) * 0.7 + 0.3  # Bright colors (0.3-1.0 range)

        # Color all pixels belonging to this cell
        cell_pixels = mask == cell_id
        colored_mask[cell_pixels] = color_rgb

    print(f"  ✓ Colored {len(unique_mask_ids)} cells with random distinct colors")

    return colored_mask, {}, {}  # Return empty dicts for compatibility


# ============================================================================
# NAPARI VISUALIZATION
# ============================================================================

def visualize_imc_with_segmentation(
    imc_img: np.ndarray,
    mask: np.ndarray,
    adata: ad.AnnData,
    roi_id: str,
    slide_scene: str = None,
    channels_to_show: List[int] = [11, 13, 35, 2],  # DAPI1, ER, panCK, CD3
    channel_colors: List[str] = ['blue', 'red', 'green', 'yellow'],  # Custom colors
    title: str = "IMC with Segmentation"
):
    """
    Visualize IMC image with segmentation masks in Napari

    Parameters:
    -----------
    imc_img : (C, H, W) multichannel IMC image
    mask : (H, W) segmentation mask with cell IDs
    adata : AnnData with cell metadata
    roi_id : ROI identifier (e.g., "X3Y1")
    slide_scene : Full slide_scene value from adata (e.g., "257_X3Y1-28")
    channels_to_show : List of channel indices to display [DAPI1, ER, panCK, CD3]
    channel_colors : List of colormap names for each channel
    """
    print(f"\n{'='*70}")
    print(f"🔬 Opening Napari: {title}")
    print(f"{'='*70}")

    # Create Napari viewer
    viewer = napari.Viewer(title=title)

    # Add IMC channels with custom colors
    channel_names = IMC_MARKER_NAMES[:imc_img.shape[0]]

    print("\n📺 Adding IMC channels with custom colors...")
    for idx, (ch_idx, colormap) in enumerate(zip(channels_to_show, channel_colors)):
        if ch_idx >= imc_img.shape[0]:
            continue

        channel = imc_img[ch_idx]
        marker_name = channel_names[ch_idx] if ch_idx < len(channel_names) else f"Channel_{ch_idx}"
        display_name = get_display_name(marker_name)

        # Normalize channel for display
        channel_norm = channel.astype(float)
        p2, p98 = np.percentile(channel_norm, [2, 98])
        if p98 > p2:
            channel_norm = np.clip((channel_norm - p2) / (p98 - p2), 0, 1)

        viewer.add_image(
            channel_norm,
            name=display_name,
            colormap=colormap,
            blending='additive',
            contrast_limits=[0, 1],
            visible=True  # Show all channels
        )
        print(f"  ✓ {display_name} (channel {ch_idx}) → {colormap.upper()}")

    # Add segmentation mask (uncolored)
    print("\n🎭 Adding segmentation mask...")
    viewer.add_labels(
        mask,
        name='Cell Segmentation',
        opacity=0.3,
        visible=True
    )

    # Add cell type colored mask
    print("\n🎨 Adding cell type overlay...")
    colored_mask, cell_map, colors_used = create_celltype_colored_mask(
        mask, adata, roi_id, slide_scene=slide_scene
    )

    viewer.add_image(
        colored_mask,
        name='Segmented Cells',
        opacity=0.6,
        visible=True,
        rgb=True
    )

    # Print legend
    print("\n📊 Cell Type Legend:")
    for celltype, color in sorted(colors_used.items()):
        print(f"  ● {celltype}: RGB{tuple(color)}")

    print(f"\n{'='*70}")
    print("✨ Napari viewer opened!")
    print("='*70}")
    print("\n💡 Controls:")
    print("  - Toggle layers with eye icon (👁️)")
    print("  - Adjust opacity with slider")
    print("  - Zoom: mouse wheel")
    print("  - Pan: click + drag")
    print("  - Layer controls: left panel")
    print("\n  Close window to continue...")
    print(f"{'='*70}\n")

    napari.run()


def create_rgb_merge(
    imc_img: np.ndarray,
    channels_to_merge: List[int],
    channel_colors: List[str]
) -> np.ndarray:
    """
    Create RGB composite merge of multiple channels with additive blending

    Parameters:
    -----------
    imc_img : (C, H, W) multichannel image
    channels_to_merge : List of channel indices [DAPI1, ER, panCK, CD3]
    channel_colors : List of colors ['blue', 'red', 'green', 'yellow']

    Returns:
    --------
    rgb_merge : (H, W, 3) RGB composite image
    """
    h, w = imc_img.shape[1], imc_img.shape[2]
    rgb_merge = np.zeros((h, w, 3), dtype=np.float32)

    # Color mapping (RGB channels)
    color_map = {
        'blue': [0, 0, 1],
        'red': [1, 0, 0],
        'green': [0, 1, 0],
        'yellow': [1, 1, 0],  # R + G
        'cyan': [0, 1, 1],    # G + B
        'magenta': [1, 0, 1]  # R + B
    }

    for ch_idx, color in zip(channels_to_merge, channel_colors):
        if ch_idx >= imc_img.shape[0]:
            continue

        # Get channel and normalize
        channel = imc_img[ch_idx]
        p2, p98 = np.percentile(channel, [2, 98])
        if p98 > p2:
            channel_norm = np.clip((channel - p2) / (p98 - p2), 0, 1)
        else:
            channel_norm = channel / (channel.max() + 1e-6)

        # Add to RGB channels with additive blending
        rgb_weights = color_map.get(color, [1, 1, 1])
        for i in range(3):
            rgb_merge[:, :, i] += channel_norm * rgb_weights[i]

    # Clip to [0, 1] range
    rgb_merge = np.clip(rgb_merge, 0, 1)

    return rgb_merge


def export_annotated_snapshot(
    imc_img: np.ndarray,
    mask: np.ndarray,
    adata: ad.AnnData,
    roi_id: str,
    slide_scene: str,
    channels_to_show: List[int],
    output_path: Path,
    channel_colors: List[str] = ['Blues', 'Reds', 'Greens', 'YlOrBr']
):
    """
    Export annotated snapshot (non-interactive)

    Parameters:
    -----------
    channel_colors : List of matplotlib colormap names (e.g., 'Blues', 'Reds', 'Greens')
    """
    print(f"\n📸 Exporting snapshot: {output_path.name}")

    # Create figure (channels + cell types overlay)
    n_channels = len(channels_to_show)
    fig, axes = plt.subplots(1, n_channels + 1, figsize=(4*(n_channels+1), 4))

    # Plot channels with specific colormaps
    for i, ch_idx in enumerate(channels_to_show):
        if ch_idx >= imc_img.shape[0]:
            continue

        channel = imc_img[ch_idx]

        # Normalize using percentiles (better contrast)
        p2, p98 = np.percentile(channel, [2, 98])
        if p98 > p2:
            channel_norm = np.clip((channel - p2) / (p98 - p2), 0, 1)
        else:
            channel_norm = channel / (channel.max() + 1e-6)

        # Use custom colormap (black to color)
        color_name = channel_colors[i] if i < len(channel_colors) else 'gray'
        cmap = create_black_to_color_cmap(color_name)

        axes[i].imshow(channel_norm, cmap=cmap, interpolation='nearest')
        marker_name = IMC_MARKER_NAMES[ch_idx] if ch_idx < len(IMC_MARKER_NAMES) else f"Ch{ch_idx}"
        display_name = get_display_name(marker_name)
        axes[i].set_title(display_name, fontsize=12, fontweight='bold')
        axes[i].axis('off')

    # Plot cell type overlay (segmentation already shown here by color)
    colored, _, _ = create_celltype_colored_mask(mask, adata, roi_id, slide_scene=slide_scene)
    axes[n_channels].imshow(colored)
    axes[n_channels].set_title('Segmented Cells', fontsize=12, fontweight='bold')
    axes[n_channels].axis('off')

    # Set black background and white text (fluorescence microscopy style)
    fig.patch.set_facecolor('black')
    plt.suptitle(f'IMC ROI: {roi_id} (Patient {slide_scene.split("_")[0]})',
                fontsize=16, fontweight='bold', color='white')

    # Set title color for each subplot
    for ax in axes:
        ax.title.set_color('white')

    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='black')
    plt.close()

    print(f"  ✓ Saved: {output_path}")


def find_interesting_region(
    imc_img: np.ndarray,
    mask: np.ndarray,
    zoom_size: int = 250
) -> Tuple[int, int, int, int]:
    """
    Find an interesting region with high cell density for zoom

    Parameters:
    -----------
    imc_img : (C, H, W) multichannel image
    mask : (H, W) segmentation mask
    zoom_size : Size of zoom region (width and height)

    Returns:
    --------
    (x, y, width, height) : Coordinates of interesting region
    """
    h, w = mask.shape

    # Create density map by counting cells in sliding windows
    step = 50
    max_density = 0
    best_x, best_y = 0, 0

    for y in range(0, h - zoom_size, step):
        for x in range(0, w - zoom_size, step):
            # Count unique cells in this region
            region = mask[y:y+zoom_size, x:x+zoom_size]
            n_cells = len(np.unique(region)) - 1  # -1 for background

            if n_cells > max_density:
                max_density = n_cells
                best_x, best_y = x, y

    return best_x, best_y, zoom_size, zoom_size


def export_zoom_figures(
    imc_img: np.ndarray,
    mask: np.ndarray,
    adata: ad.AnnData,
    roi_id: str,
    slide_scene: str,
    patient_id: str,
    channels_to_show: List[int],
    channel_colors: List[str],
    output_dir: Path,
    zoom_coords: Tuple[int, int, int, int] = None
):
    """
    Export zoomed versions of annotated and merge figures

    Parameters:
    -----------
    zoom_coords : (x, y, width, height) or None for auto-detect
    """
    print(f"\n{'='*70}")
    print("🔍 Generating Zoom Views")
    print("="*70)

    # Find or use provided zoom region
    if zoom_coords is None:
        print("\n  Auto-detecting interesting region...")
        x, y, w, h = find_interesting_region(imc_img, mask, zoom_size=250)
        print(f"  ✓ Found region: x={x}, y={y}, size={w}×{h}")
    else:
        x, y, w, h = zoom_coords
        print(f"\n  Using specified region: x={x}, y={y}, size={w}×{h}")

    # Extract zoom regions
    imc_zoom = imc_img[:, y:y+h, x:x+w]
    mask_zoom = mask[y:y+h, x:x+w]

    # Count cells in zoom region
    n_cells_zoom = len(np.unique(mask_zoom)) - 1
    print(f"  ✓ Zoom region contains {n_cells_zoom} cells")

    # === ZOOM 1: Annotated panels ===
    print(f"\n📸 Creating zoom of annotated panels...")

    n_channels = len(channels_to_show)
    fig, axes = plt.subplots(1, n_channels + 1, figsize=(4*(n_channels+1), 4))

    # Plot zoomed channels
    for i, ch_idx in enumerate(channels_to_show):
        if ch_idx >= imc_zoom.shape[0]:
            continue

        channel = imc_zoom[ch_idx]
        p2, p98 = np.percentile(channel, [2, 98])
        if p98 > p2:
            channel_norm = np.clip((channel - p2) / (p98 - p2), 0, 1)
        else:
            channel_norm = channel / (channel.max() + 1e-6)

        # Use custom colormap (black to color)
        color_name = channel_colors[i] if i < len(channel_colors) else 'gray'
        cmap = create_black_to_color_cmap(color_name)
        axes[i].imshow(channel_norm, cmap=cmap, interpolation='nearest')
        marker_name = IMC_MARKER_NAMES[ch_idx] if ch_idx < len(IMC_MARKER_NAMES) else f"Ch{ch_idx}"
        display_name = get_display_name(marker_name)
        axes[i].set_title(display_name, fontsize=12, fontweight='bold')
        axes[i].axis('off')

    # Plot zoomed cell types
    colored_zoom, _, _ = create_celltype_colored_mask(mask_zoom, adata, roi_id, slide_scene=slide_scene)
    axes[n_channels].imshow(colored_zoom)
    axes[n_channels].set_title('Segmented Cells', fontsize=12, fontweight='bold')
    axes[n_channels].axis('off')

    # Set black background and white text
    fig.patch.set_facecolor('black')
    plt.suptitle(f'IMC ROI: {roi_id} | Patient: {patient_id} | ZOOM ({n_cells_zoom} cells)',
                fontsize=16, fontweight='bold', color='white')

    # Set title color for each subplot
    for ax in axes:
        ax.title.set_color('white')

    plt.tight_layout()

    zoom_annotated_path = output_dir / f"roi_{roi_id}_patient_{patient_id}_annotated_ZOOM.png"
    plt.savefig(zoom_annotated_path, dpi=300, bbox_inches='tight', facecolor='black')
    plt.close()
    print(f"  ✓ Saved: {zoom_annotated_path.name}")

    # === ZOOM 2: RGB Merge + Segmented Cells ===
    print(f"\n🎨 Creating zoom of RGB merge with cell types...")

    rgb_merge_zoom = create_rgb_merge(imc_zoom, channels_to_show, channel_colors)

    # Create figure with 2 subplots: RGB merge + Cell types
    fig, axes = plt.subplots(1, 2, figsize=(20, 10))

    # Panel 1: RGB Merge
    axes[0].imshow(rgb_merge_zoom)
    axes[0].axis('off')

    channel_names_list = [get_display_name(IMC_MARKER_NAMES[i]) if i < len(IMC_MARKER_NAMES) else f"Ch{i}"
                          for i in channels_to_show]
    legend_text = " + ".join([f"{name} ({color})"
                             for name, color in zip(channel_names_list, channel_colors)])
    axes[0].set_title(f'RGB Merge\n{legend_text}', fontsize=12, pad=10, color='white')

    # Panel 2: Segmented Cells
    axes[1].imshow(colored_zoom)
    axes[1].axis('off')
    axes[1].set_title('Segmented Cells', fontsize=12, pad=10, color='white')

    # Set black background and white text
    fig.patch.set_facecolor('black')
    plt.suptitle(f'IMC Composite | ZOOM\nROI: {roi_id} | Patient: {patient_id} | {n_cells_zoom} cells',
                fontsize=16, fontweight='bold', color='white')
    plt.tight_layout()

    zoom_merge_path = output_dir / f"roi_{roi_id}_patient_{patient_id}_merge_ZOOM.png"
    plt.savefig(zoom_merge_path, dpi=300, bbox_inches='tight', facecolor='black')
    plt.close()
    print(f"  ✓ Saved: {zoom_merge_path.name}")

    print(f"\n{'='*70}")
    print(f"✅ Zoom views generated")
    print(f"{'='*70}")
    print(f"  Region: x={x}, y={y}, size={w}×{h}")
    print(f"  Cells in zoom: {n_cells_zoom}")


def export_merge_figure(
    imc_img: np.ndarray,
    channels_to_merge: List[int],
    channel_colors: List[str],
    output_path: Path,
    roi_id: str,
    patient_id: str
):
    """
    Export RGB composite merge figure

    Parameters:
    -----------
    imc_img : (C, H, W) multichannel image
    channels_to_merge : List of channel indices
    channel_colors : List of color names
    output_path : Path to save figure
    roi_id : ROI identifier
    patient_id : Patient ID
    """
    print(f"\n🎨 Creating RGB merge composite: {output_path.name}")

    # Create RGB merge
    rgb_merge = create_rgb_merge(imc_img, channels_to_merge, channel_colors)

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    ax.imshow(rgb_merge)
    ax.axis('off')

    # Create channel legend
    channel_names_list = [get_display_name(IMC_MARKER_NAMES[i]) if i < len(IMC_MARKER_NAMES) else f"Ch{i}"
                          for i in channels_to_merge]
    legend_text = " + ".join([f"{name} ({color})"
                             for name, color in zip(channel_names_list, channel_colors)])

    # Set black background and white text
    fig.patch.set_facecolor('black')
    plt.suptitle(f'IMC RGB Composite\nROI: {roi_id} | Patient: {patient_id}',
                fontsize=16, fontweight='bold', color='white')
    ax.set_title(legend_text, fontsize=12, pad=10, color='white')

    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='black')
    plt.close()

    print(f"  ✓ Saved RGB merge: {output_path}")
    print(f"  ✓ Channels merged: {legend_text}")


def export_biological_view(
    imc_img: np.ndarray,
    mask: np.ndarray,
    adata: ad.AnnData,
    roi_id: str,
    slide_scene: str,
    patient_id: str,
    view_key: str,
    output_dir: Path,
    export_zoom: bool = True
):
    """
    Export a predefined biological composite view

    Parameters:
    -----------
    view_key : str
        Key from BIOLOGICAL_VIEWS dict (e.g., 'tumor_architecture')
    export_zoom : bool
        Whether to also export zoomed version
    """
    if view_key not in BIOLOGICAL_VIEWS:
        print(f"⚠️  Unknown view: {view_key}")
        return

    view = BIOLOGICAL_VIEWS[view_key]

    print(f"\n{'='*70}")
    print(f"📊 Exporting View: {view['name']}")
    print(f"{'='*70}")
    print(f"  Markers: {', '.join(view['markers'])}")
    print(f"  Colors: {', '.join(view['colors'])}")
    print(f"  Connects to: {view['connects_to']}")

    channels = view['indices']
    colors = view['colors']

    # === FULL VIEW: Annotated panels ===
    print(f"\n📸 Creating annotated panels...")

    n_channels = len(channels)
    fig, axes = plt.subplots(1, n_channels + 1, figsize=(4*(n_channels+1), 4))

    # Plot channels
    for i, ch_idx in enumerate(channels):
        if ch_idx >= imc_img.shape[0]:
            continue

        channel = imc_img[ch_idx]
        p2, p98 = np.percentile(channel, [2, 98])
        if p98 > p2:
            channel_norm = np.clip((channel - p2) / (p98 - p2), 0, 1)
        else:
            channel_norm = channel / (channel.max() + 1e-6)

        # Use custom colormap
        color_name = colors[i]
        cmap = create_black_to_color_cmap(color_name)
        axes[i].imshow(channel_norm, cmap=cmap, interpolation='nearest')

        marker_name = view['markers'][i]
        axes[i].set_title(f"{marker_name}\n({color_name})",
                         fontsize=12, fontweight='bold', color='white')
        axes[i].axis('off')

    # Plot cell types
    colored, _, _ = create_celltype_colored_mask(mask, adata, roi_id, slide_scene=slide_scene)
    axes[n_channels].imshow(colored)
    axes[n_channels].set_title('Segmented Cells', fontsize=12, fontweight='bold', color='white')
    axes[n_channels].axis('off')

    # Set black background
    fig.patch.set_facecolor('black')
    plt.suptitle(f'{view["name"]}\nROI: {roi_id} | Patient: {patient_id}',
                fontsize=16, fontweight='bold', color='white')

    plt.tight_layout()

    output_path = output_dir / f"roi_{roi_id}_{view_key}_annotated.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='black')
    plt.close()
    print(f"  ✓ Saved: {output_path.name}")

    # === FULL VIEW: RGB Merge ===
    print(f"\n🎨 Creating RGB merge...")

    rgb_merge = create_rgb_merge(imc_img, channels, colors)

    fig, axes = plt.subplots(1, 2, figsize=(20, 10))

    # Panel 1: RGB Merge
    axes[0].imshow(rgb_merge)
    axes[0].axis('off')
    legend_text = " + ".join([f"{m} ({c})" for m, c in zip(view['markers'], colors)])
    axes[0].set_title(f'Composite\n{legend_text}', fontsize=12, pad=10, color='white')

    # Panel 2: Segmented Cells
    axes[1].imshow(colored)
    axes[1].axis('off')
    axes[1].set_title('Segmented Cells', fontsize=12, pad=10, color='white')

    fig.patch.set_facecolor('black')
    plt.suptitle(f'{view["name"]}\nROI: {roi_id} | Patient: {patient_id}',
                fontsize=16, fontweight='bold', color='white')
    plt.tight_layout()

    merge_path = output_dir / f"roi_{roi_id}_{view_key}_merge.png"
    plt.savefig(merge_path, dpi=300, bbox_inches='tight', facecolor='black')
    plt.close()
    print(f"  ✓ Saved: {merge_path.name}")

    # === ZOOM VERSION ===
    if export_zoom:
        print(f"\n🔍 Creating zoom view...")

        # Find interesting region
        x, y, w, h = find_interesting_region(imc_img, mask, zoom_size=250)
        imc_zoom = imc_img[:, y:y+h, x:x+w]
        mask_zoom = mask[y:y+h, x:x+w]
        n_cells_zoom = len(np.unique(mask_zoom)) - 1

        # Zoom annotated
        fig, axes = plt.subplots(1, n_channels + 1, figsize=(4*(n_channels+1), 4))

        for i, ch_idx in enumerate(channels):
            if ch_idx >= imc_zoom.shape[0]:
                continue

            channel = imc_zoom[ch_idx]
            p2, p98 = np.percentile(channel, [2, 98])
            if p98 > p2:
                channel_norm = np.clip((channel - p2) / (p98 - p2), 0, 1)
            else:
                channel_norm = channel / (channel.max() + 1e-6)

            color_name = colors[i]
            cmap = create_black_to_color_cmap(color_name)
            axes[i].imshow(channel_norm, cmap=cmap, interpolation='nearest')
            axes[i].set_title(f"{view['markers'][i]}\n({color_name})",
                             fontsize=12, fontweight='bold', color='white')
            axes[i].axis('off')

        colored_zoom, _, _ = create_celltype_colored_mask(mask_zoom, adata, roi_id, slide_scene=slide_scene)
        axes[n_channels].imshow(colored_zoom)
        axes[n_channels].set_title('Segmented Cells', fontsize=12, fontweight='bold', color='white')
        axes[n_channels].axis('off')

        fig.patch.set_facecolor('black')
        plt.suptitle(f'{view["name"]} | ZOOM\nROI: {roi_id} | Patient: {patient_id} | {n_cells_zoom} cells',
                    fontsize=16, fontweight='bold', color='white')
        plt.tight_layout()

        zoom_path = output_dir / f"roi_{roi_id}_{view_key}_annotated_ZOOM.png"
        plt.savefig(zoom_path, dpi=300, bbox_inches='tight', facecolor='black')
        plt.close()
        print(f"  ✓ Saved: {zoom_path.name}")

        # Zoom merge
        rgb_merge_zoom = create_rgb_merge(imc_zoom, channels, colors)

        fig, axes = plt.subplots(1, 2, figsize=(20, 10))
        axes[0].imshow(rgb_merge_zoom)
        axes[0].axis('off')
        axes[0].set_title(f'Composite\n{legend_text}', fontsize=12, pad=10, color='white')

        axes[1].imshow(colored_zoom)
        axes[1].axis('off')
        axes[1].set_title('Segmented Cells', fontsize=12, pad=10, color='white')

        fig.patch.set_facecolor('black')
        plt.suptitle(f'{view["name"]} | ZOOM\nROI: {roi_id} | Patient: {patient_id} | {n_cells_zoom} cells',
                    fontsize=16, fontweight='bold', color='white')
        plt.tight_layout()

        zoom_merge_path = output_dir / f"roi_{roi_id}_{view_key}_merge_ZOOM.png"
        plt.savefig(zoom_merge_path, dpi=300, bbox_inches='tight', facecolor='black')
        plt.close()
        print(f"  ✓ Saved: {zoom_merge_path.name}")


def batch_export_all_views(
    imc_img: np.ndarray,
    mask: np.ndarray,
    adata: ad.AnnData,
    roi_id: str,
    slide_scene: str,
    patient_id: str,
    output_dir: Path,
    export_zoom: bool = True
):
    """
    Export all 6 biological composite views for a single ROI

    This creates a complete visualization set for workshop demonstrations,
    showing all major biological axes relevant to downstream analysis.
    """
    print(f"\n{'='*70}")
    print(f"🎯 BATCH EXPORT: All Biological Views")
    print(f"{'='*70}")
    print(f"  ROI: {roi_id}")
    print(f"  Patient: {patient_id}")
    print(f"  Output: {output_dir}/")
    print(f"  Views: {len(BIOLOGICAL_VIEWS)}")

    for i, (view_key, view_info) in enumerate(BIOLOGICAL_VIEWS.items(), 1):
        print(f"\n[{i}/{len(BIOLOGICAL_VIEWS)}] Processing: {view_info['name']}")

        export_biological_view(
            imc_img, mask, adata, roi_id, slide_scene, patient_id,
            view_key, output_dir, export_zoom=export_zoom
        )

    print(f"\n{'='*70}")
    print(f"✅ BATCH EXPORT COMPLETE")
    print(f"{'='*70}")
    print(f"  Total views exported: {len(BIOLOGICAL_VIEWS)}")
    print(f"  Files per view: 4 (annotated, merge, 2 zooms)" if export_zoom else "  Files per view: 2")
    print(f"  Output directory: {output_dir}/")


# ============================================================================
# MAIN
# ============================================================================

def main():
    """
    Main visualization workflow for 3 workshop patients

    Generates 6 biological views for patients 58, 108, 144
    Used as STEP 1 in run_pipeline_complete_v2.sh
    """

    print("""
╔════════════════════════════════════════════════════════════════════╗
║                                                                    ║
║  STEP 1: IMC Visualizations (QC Check)                            ║
║  Workshop Patients: 58, 108, 144                                  ║
║                                                                    ║
║  Generates 6 biological views per patient:                        ║
║  • Tumor Architecture  • Immune Landscape  • Stromal Organization ║
║  • Invasive Edge       • Vascular Proximity • Proliferation       ║
║                                                                    ║
╚════════════════════════════════════════════════════════════════════╝
    """)

    # Workshop patients (3 selected for demonstration)
    WORKSHOP_PATIENTS = [
        {'id': '58', 'roi': 'X15Y1', 'slide': 'SP41', 'scene': '58_X15Y1-60', 'cells': 4347, 'coverage': 70.9},
        {'id': '108', 'roi': 'X13Y8', 'slide': 'SP43', 'scene': '108_X13Y8-183', 'cells': 6320, 'coverage': 95.6},
        {'id': '144', 'roi': 'X15Y1', 'slide': 'SP43', 'scene': '144_X15Y1-195', 'cells': 6267, 'coverage': 69.3},
    ]

    # Load AnnData once (shared across all patients)
    adata_path = PATHS['adata']  
    print("\nLoading AnnData...")
    adata = load_imc_anndata(adata_path)
    print(f"  ✓ Loaded {len(adata.obs):,} cells")


    # Process each patient
    for i, patient_info in enumerate(WORKSHOP_PATIENTS, 1):
        patient_id = patient_info['id']
        roi_id = patient_info['roi']
        slide = patient_info['slide']
        scene = patient_info['scene']
        cells = patient_info['cells']
        coverage = patient_info['coverage']

        print(f"\n{'='*70}")
        print(f"[{i}/{len(WORKSHOP_PATIENTS)}] Patient {patient_id} ({slide}, {roi_id})")
        print(f"{'='*70}")
        print(f"  Coverage: {coverage}%")
        print(f"  Expected cells: {cells:,}")

        # Find image and mask files
        print(f"\n🔍 Searching for {slide} Patient {patient_id} files...")
        imc_files = list(PATHS['imc_images'].glob(f"*{slide}*{patient_id}*{roi_id}*.tiff"))
        mask_files = list(PATHS['cell_masks'].glob(f"*{slide}*{patient_id}*{roi_id}*_maks.tiff"))

        if not imc_files or not mask_files:
            print(f"  ⚠️  Files not found!")
            continue

        imc_path = imc_files[0]
        mask_path = mask_files[0]

        print(f"  ✓ Found image: {imc_path.name}")
        print(f"  ✓ Found mask: {mask_path.name}")

        # Load image and mask
        imc_img = load_imc_image(imc_path)
        mask = load_segmentation_mask(mask_path)

        # Verify dimensions
        if imc_img.shape[1:] != mask.shape:
            print(f"\n⚠️  Dimension mismatch!")
            print(f"  Image: {imc_img.shape[1:]} (H, W)")
            print(f"  Mask: {mask.shape} (H, W)")
            continue

        print(f"  ✓ Dimensions match: {mask.shape}")

        # Create patient-specific output directory
        output_dir = PATHS['output'] / f'patient_{patient_id}'
        output_dir.mkdir(parents=True, exist_ok=True)

        # Export all biological views
        print(f"\n📊 Generating 6 biological views...")
        batch_export_all_views(
            imc_img, mask, adata, roi_id, scene, patient_id,
            output_dir, export_zoom=True
        )

        print(f"\n  ✅ Patient {patient_id} complete!")
        print(f"  📁 Output: {output_dir.relative_to(BASE_DIR)}/")

    print(f"\n\n{'='*70}")
    print(f"✅ STEP 1 COMPLETE: ALL VISUALIZATIONS GENERATED")
    print(f"{'='*70}")
    print(f"\n📊 Summary:")
    print(f"  Total patients: {len(WORKSHOP_PATIENTS)}")
    print(f"  Views per patient: 6")
    print(f"  Files per patient: 24 (6 views × 4 files)")
    print(f"  Total files: {len(WORKSHOP_PATIENTS) * 24} PNG images")
    print(f"\n📁 Output directory: {PATHS['output'].relative_to(BASE_DIR)}/")
    print(f"\n💡 Next pipeline steps:")
    print(f"  → STEP 2: Spatial neighbor graphs (erpos_step2_neighbors.py)")
    print(f"  → STEP 3: Targeted neighborhood composition")
    print(f"  → ... (see run_pipeline_complete_v2.sh for full pipeline)")


if __name__ == "__main__":
    main()
