import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Optional

# --- Publication-Quality Typography & Styling ---
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica', 'Liberation Sans'],
    'font.size': 9,                  # Standard paper font size (8-10 pt)
    'axes.titlesize': 10,
    'axes.titleweight': 'bold',
    'axes.labelsize': 9,
    'figure.titlesize': 11,
    'pdf.fonttype': 42,               # Embeds fonts as TrueType in PDF/EPS for publishers
    'ps.fonttype': 42
})

def normalize_volume(data: np.ndarray, method: str = 'minmax', ignore_zeros: bool = True) -> np.ndarray:
    """Normalizes NIfTI volume data safely."""
    data = data.astype(np.float32)
    mask = data > 0 if ignore_zeros else np.ones_like(data, dtype=bool)
    
    if not np.any(mask):
        return data

    if method == 'minmax':
        min_val, max_val = np.min(data[mask]), np.max(data[mask])
        if max_val - min_val > 0:
            normalized = (data - min_val) / (max_val - min_val)
            normalized[~mask] = 0
            return np.clip(normalized, 0, 1)
        return data
    elif method == 'zscore':
        mean_val, std_val = np.mean(data[mask]), np.std(data[mask])
        if std_val > 0:
            normalized = (data - mean_val) / std_val
            normalized[~mask] = 0
            return normalized
        return data
    else:
        raise ValueError("Unsupported method. Use 'minmax' or 'zscore'.")


def crop_to_match(data1: np.ndarray, data2: np.ndarray, center: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """Crops two arrays to their shared minimum dimensions."""
    target_shape = [min(s1, s2) for s1, s2 in zip(data1.shape, data2.shape)]

    def get_slices(arr_shape):
        slices = []
        for current_dim, target_dim in zip(arr_shape, target_shape):
            start = (current_dim - target_dim) // 2 if center else 0
            slices.append(slice(start, start + target_dim))
        return tuple(slices)

    return data1[get_slices(data1.shape)], data2[get_slices(data2.shape)]


def plot_multislice_nifti_comparison(
    file1_path: str,
    file2_path: str,
    slice_indices: List[int],
    method1_name: str = "Recon Method A",
    method2_name: str = "Recon Method B",
    axis: int = 2,
    flip_image2: bool = False,
    flip_axis: int = 0,
    normalize: bool = True,
    norm_method: str = 'minmax',
    output_basename: str = "fig_recon_slice_comparison"
):
    """
    Plots a 3-row comparison matrix from two NIfTI reconstruction volumes using 3 different slice indices.
    
    :param file1_path: Path to first NIfTI file (Recon 1)
    :param file2_path: Path to second NIfTI file (Recon 2)
    :param slice_indices: List of exactly 3 integer slice indices (e.g., [70, 90, 110])
    :param method1_name: Name of Method 1 for column header
    :param method2_name: Name of Method 2 for column header
    :param axis: Slice plane (0 = Sagittal, 1 = Coronal, 2 = Axial)
    """
    if len(slice_indices) != 3:
        raise ValueError("Please provide exactly 3 slice indices in 'slice_indices'.")

    # 1. Load NIfTI volumes
    img1_obj = nib.load(file1_path)
    img2_obj = nib.load(file2_path)
    
    data1 = img1_obj.get_fdata()
    data2 = img2_obj.get_fdata()

    data1 = data1[170:232, 120:220, :]
    data2 = data2[:, 121:221, :]

    # 2. Fix mirrored volume if needed
    if flip_image2:
        data2 = np.flip(data2, axis=flip_axis)

    # 3. Apply normalization across 3D volumes
    if normalize:
        data1 = normalize_volume(data1, method=norm_method)
        data2 = normalize_volume(data2, method=norm_method)

    # 4. Crop matching dimensions if shape mismatches exist
    if data1.shape != data2.shape:
        print(f"Cropping shapes from {data1.shape} and {data2.shape} to match.")
        data1, data2 = crop_to_match(data1, data2, center=True)

    # 5. Extract 2D slices & compute signed differences
    processed_slices1 = []
    processed_slices2 = []
    processed_diffs = []
    global_max_diff = 0.0

    for s_idx in slice_indices:
        s1 = np.take(data1, s_idx, axis=axis)
        s2 = np.take(data2, s_idx, axis=axis)
        diff = s1 - s2

        max_diff = np.max(np.abs(diff))
        if max_diff > global_max_diff:
            global_max_diff = max_diff

        processed_slices1.append(s1)
        processed_slices2.append(s2)
        processed_diffs.append(diff)

    if global_max_diff == 0:
        global_max_diff = 1.0

    # 6. Plot Matrix Setup (Two-column paper width ~ 7.2 inches)
    fig, axes = plt.subplots(3, 3, figsize=(7.2, 7.5), dpi=300)
    panel_tags = ["(a)", "(b)", "(c)"]


    for i in range(3):
        s1 = processed_slices1[2-i]
        s2 = processed_slices2[2-i]
        diff = processed_diffs[2-i]
        s_idx = slice_indices[2-i]

        # Recon Method 1
        axes[i, 0].imshow(s1.T, cmap="gray", origin="lower")
        axes[i, 0].axis("off")

        # Recon Method 2
        axes[i, 1].imshow(s2.T, cmap="gray", origin="lower")
        axes[i, 1].axis("off")

        # Signed Difference
        diff_img = axes[i, 2].imshow(
            diff.T, 
            cmap="coolwarm", 
            origin="lower", 
            vmin=-global_max_diff, 
            vmax=global_max_diff
        )
        axes[i, 2].axis("off")

        # Subpanel Label (a), (b), (c) + Slice Number
        axes[i, 0].text(
            -0.12, 0.5, f"{panel_tags[i]}  Slice {s_idx}", 
            transform=axes[i, 0].transAxes, 
            rotation=90, 
            va='center', 
            ha='right', 
            fontsize=9.5, 
            fontweight='bold'
        )

    # Column Headers
    axes[0, 0].set_title(method1_name, pad=8)
    axes[0, 1].set_title(method2_name, pad=8)
    axes[0, 2].set_title(f"Difference\n({method1_name} − {method2_name})", pad=8, fontsize=9)

    # Adjust spacing for publication layout
    plt.subplots_adjust(wspace=0.01, hspace=0.04, right=0.86, left=0.12, top=0.93, bottom=0.02)

    # Shared Colorbar
    cbar_ax = fig.add_axes([0.88, 0.08, 0.022, 0.82])
    cbar = fig.colorbar(diff_img, cax=cbar_ax)
    
    label_unit = "Rescaled Intensity (0..1)" if norm_method == 'minmax' else "Z-score Difference"
    cbar.set_label(f"Difference [{label_unit}]", rotation=270, labelpad=14, fontsize=9, fontweight='bold')
    cbar.ax.tick_params(labelsize=8)

    # Export vector (PDF) and raster (PNG)
    plt.savefig(f"{output_basename}.pdf", format="pdf", bbox_inches='tight', pad_inches=0.02)
    plt.savefig(f"{output_basename}.png", format="png", dpi=300, bbox_inches='tight', pad_inches=0.02)

    plt.show()

# ==========================================
# Example Usage:
# ==========================================
plot_multislice_nifti_comparison(
     file1_path="/home/rhiannon.b/Downloads/ReconData/sub9027/thoracic/niftis/tmp_dcm2bids/sub-9027/normorraw/rep_averaged/rms_4echoes_3D_fixed.nii.gz",
     file2_path="/home/rhiannon.b/Downloads/ReconData/SLR_recon_results/sub9027/thoracic/MID00072_thoracic/combined/echo_rms_1to4.nii.gz",
     slice_indices=[2, 6,10],              # <-- Your 3 slice choices here
     method1_name="Standard",
     method2_name="SLR",
     normalize=True,
     norm_method='minmax',
     output_basename="fig_recon_comparison"
 )