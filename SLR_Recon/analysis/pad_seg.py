import nibabel as nib
import numpy as np


def pad_and_orient_mask(
    siemens_path,
    cropped_seg_path,
    out_path,
    crop_indices,
    reverse_z=True,
    flip_y=True,
):
    """Reverses slice order, applies any in-plane flips, and embeds the cropped

    segmentation into the full Siemens volume matrix using manual indices.
    """
    # 1. Load files
    siemens_nii = nib.load(siemens_path)
    seg_nii = nib.load(cropped_seg_path)

    siemens_data = siemens_nii.get_fdata()
    seg_data = seg_nii.get_fdata()

    x_min, x_max = crop_indices

    # 2. Reverse slice order along Z-axis (axis=2)
    if reverse_z:
        seg_data = np.flip(seg_data, axis=2)
        print("  -> Reversed slice order along Z-axis")

    if flip_y:
        seg_data = np.flip(seg_data, axis=0)
        print("  -> Flipped along Y-axis")

    # 4. Create empty full-sized volume
    padded_seg = np.zeros(siemens_data.shape[:3], dtype=seg_data.dtype)

    # 5. Insert segmentation at the designated coordinates
    padded_seg[x_min:x_max+1, :, :] = seg_data

    # 6. Save using full Siemens header and affine
    out_nii = nib.Nifti1Image(
        padded_seg.astype(np.uint8), siemens_nii.affine, siemens_nii.header
    )
    nib.save(out_nii, out_path)
    print(f"--> Saved padded and aligned mask to: {out_path}\n")


# Crop bounding box: (x_min, x_max, y_min, y_max, z_min, z_max)
CROP_INDICES = (161, 221)

pad_and_orient_mask(
    siemens_path="/home/rhiannon.b/Downloads/ReconData/sub9021/lumbar_siemens/niftis/raw/rep_averaged/rms_4echoes_3D.nii.gz",
    cropped_seg_path="/home/rhiannon.b/Downloads/ReconData/SLR_recon_results/sub9021/lumbar/sub9021_lumbar_seg.nii.gz",
    out_path="/home/rhiannon.b/Downloads/ReconData/sub9021/lumbar_siemens/niftis/raw/rep_averaged/sub9021_lumbar_gmseg.nii.gz",
    crop_indices=CROP_INDICES,
    reverse_z=True,  # Reverses slice order (0 -> N becomes N -> 0)
    flip_y=False,  # Flips vertical in-plane axis
)


