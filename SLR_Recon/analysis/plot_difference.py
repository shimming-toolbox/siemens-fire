from operator import lt

import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt

def crop_to_match(data1: np.ndarray, data2: np.ndarray, center: bool = True):
    """
    Crops two arrays to their shared minimum shape along each axis.
    """
    # 1. Determine minimum dimensions across all axes
    target_shape = [min(s1, s2) for s1, s2 in zip(data1.shape, data2.shape)]

    def get_slices(arr_shape):
        slices = []
        for current_dim, target_dim in zip(arr_shape, target_shape):
            if center:
                # Calculate start index to keep the center aligned
                start = (current_dim - target_dim) // 2
            else:
                # Crop from index 0
                start = 0
            slices.append(slice(start, start + target_dim))
        return tuple(slices)

    cropped_data1 = data1[get_slices(data1.shape)]
    cropped_data2 = data2[get_slices(data2.shape)]

    return cropped_data1, cropped_data2

def normalize_volume(data: np.ndarray, method: str = 'minmax', ignore_zeros: bool = True) -> np.ndarray:
    """
    Normalizes a NIfTI volume.
    
    :param data: Input image array
    :param method: 'minmax' (scales to 0..1) or 'zscore' (mean=0, std=1)
    :param ignore_zeros: Excludes background zero-pixels when computing stats
    """
    data = data.astype(np.float32)
    
    # Create mask for foreground non-zero voxels if requested
    mask = data > 0 if ignore_zeros else np.ones_like(data, dtype=bool)
    
    if not np.any(mask):
        return data  # Handle empty image edge case

    if method == 'minmax':
        min_val = np.min(data[mask])
        max_val = np.max(data[mask])
        if max_val - min_val > 0:
            normalized = (data - min_val) / (max_val - min_val)
            normalized[~mask] = 0  # Preserve background zero
            return np.clip(normalized, 0, 1)
        return data

    elif method == 'zscore':
        mean_val = np.mean(data[mask])
        std_val = np.std(data[mask])
        if std_val > 0:
            normalized = (data - mean_val) / std_val
            normalized[~mask] = 0  # Preserve background zero
            return normalized
        return data

    else:
        raise ValueError("Unsupported method. Use 'minmax' or 'zscore'.")


def plot_nifti_difference(file1_path: str, file2_path: str, slice_idx: int = None, axis: int = 2):
    # 1. Load NIfTI files
    img1_obj = nib.load(file1_path)
    img2_obj = nib.load(file2_path)
    
    data1 = img1_obj.get_fdata()
    data2 = img2_obj.get_fdata()
    data2 - np.fliplr(data2)  # Flip data2 left-right for better visual comparison
    data1 = data1[:,60:200,:]
    data2 = data2[157:223,60:200,:]

    data1 = normalize_volume(data1, method='minmax')
    data2 = normalize_volume(data2, method='minmax')

    # 2. Crop arrays to matching dimensions if shapes differ
    if data1.shape != data2.shape:
        print(f"Shape mismatch detected: {data1.shape} vs {data2.shape}")
        data1, data2 = crop_to_match(data1, data2, center=True)
        print(f"Cropped both images to common shape: {data1.shape}")

    # 3. Select default slice if not provided (middle slice of the cropped volume)
    if slice_idx is None:
        slice_idx = data1.shape[axis] // 2

    # 4. Extract 2D slice along chosen axis
    slice1 = np.take(data1, 1, axis=axis)
    slice2 = np.take(data2, 1, axis=axis)

    # 5. Compute SIGNED difference (no np.abs)
    diff_slice = slice1 - slice2

    # 6. Calculate symmetric range around 0
    max_val = np.max(np.abs(diff_slice))
    if max_val == 0:
        max_val = 1.0

    # 7. Plot side-by-side
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    axes[0].imshow(slice1.T, cmap="gray", origin="lower")
    axes[0].set_title(f" SLR (Thoracic)")
    axes[0].axis("off")

    axes[1].imshow(slice2.T, cmap="gray", origin="lower")
    axes[1].set_title(f" Standard (Thoracic)")
    axes[1].axis("off")


    # Coolwarm colormap centered around 0 = White
    im = axes[2].imshow(
        diff_slice.T, 
        cmap="coolwarm", 
        origin="lower", 
        vmin=-max_val, 
        vmax=max_val
    )
    axes[2].set_title(" SLR - Siemens")
    axes[2].axis("off")

    cbar = fig.colorbar(im, ax=axes[2], fraction=0.046, pad=0.04)
    cbar.set_label("Difference Magnitude (0 = White)", rotation=270, labelpad=15)

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']

    # 2. Set the global font size
    plt.rcParams['font.size'] = 14

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.00)  # Smaller number = less horizontal space (0 = touching)
    plt.show()

# Example usage:/home/rhiannon.b/Downloads/ReconData/RandomBinExperiment/seperatebins_random
plot_nifti_difference("/home/rhiannon.b/Downloads/ReconData/SLR_recon_results/MarkResults_DefaultParams/meas_MID00151_FID35300_gre_spine/combined/echo_1.nii.gz", "/home/rhiannon.b/Downloads/ReconData/Dataset3Niftis/tmp_dcm2bids/sub-datasetthree/1603803612_Dataset3Dicoms_rms_gre_spine_000000_e1.nii.gz", slice_idx=1, axis=2)


    
