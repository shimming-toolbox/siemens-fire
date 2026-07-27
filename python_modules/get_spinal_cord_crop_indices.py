import os
import shutil
import subprocess
import nibabel as nib
import numpy as np


def find_sct_path():
    # 1. Environment variable override
    if os.getenv("SCT_PATH"):
        return os.environ["SCT_PATH"]

    # 2. Check system PATH
    for cmd in ["sct_deepseg", "sct_deepseg_sc"]:
        executable = shutil.which(cmd)
        if executable:
            return executable

    # 3. Check default ~/sct/bin/ path
    for binary in ["sct_deepseg", "sct_deepseg_sc"]:
        home_sct = os.path.expanduser(f"~/sct/bin/{binary}")
        if os.path.exists(home_sct):
            return home_sct

    raise FileNotFoundError(
        "Could not find Spinal Cord Toolbox!\n"
        "Please ensure SCT is installed and added to your PATH, or set "
        "export SCT_PATH='/path/to/sct_deepseg' in your terminal."
    )


def get_indices(
    img_3d,
    affine,
    contrast="t2",
    margin_x=20,
    margin_y=5,
    x_axis=0,
    y_axis=1
):
    sct_cmd = find_sct_path()

    temp_img_path = "temp_sc_recon.nii.gz"
    temp_seg_path = "temp_sc_seg.nii.gz"

    # Save temporary NIfTI
    nifti_img = nib.Nifti1Image(img_3d.astype(np.float32), affine)
    nib.save(nifti_img, temp_img_path)

    # Build CLI argument list based on command detected
    if "sct_deepseg_sc" in sct_cmd:
        # Legacy syntax
        cmd_args = [sct_cmd, "-i", temp_img_path, "-c", contrast, "-o", temp_seg_path]
    else:
        # Modern sct_deepseg syntax (Task as positional argument)
        cmd_args = [sct_cmd, "spinalcord", "-i", temp_img_path, "-o", temp_seg_path]

    try:
        subprocess.run(cmd_args, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        # Clean up files before raising error
        for f in [temp_img_path, temp_seg_path]:
            if os.path.exists(f): os.remove(f)
        print("SCT STDERR:\n", e.stderr)
        raise e

    # Load segmentation mask
    seg_mask = nib.load(temp_seg_path).get_fdata()

    # Clean up temporary NIfTI files
    for f in [temp_img_path, temp_seg_path]:
        if os.path.exists(f): os.remove(f)

    # Calculate coordinate profiles
    x_profile = np.sum(seg_mask, axis=tuple(i for i in range(seg_mask.ndim) if i != x_axis))
    x_nonzero = np.where(x_profile > 0)[0]

    y_profile = np.sum(seg_mask, axis=tuple(i for i in range(seg_mask.ndim) if i != y_axis))
    y_nonzero = np.where(y_profile > 0)[0]

    if len(x_nonzero) == 0 or len(y_nonzero) == 0:
        raise ValueError("SCT failed to detect the spinal cord.")

    min_x = max(0, np.min(x_nonzero) - margin_x)
    max_x = min(img_3d.shape[x_axis], np.max(x_nonzero) + margin_x + 1)

    min_y = max(0, np.min(y_nonzero) - margin_y)
    max_y = min(img_3d.shape[y_axis], np.max(y_nonzero) + margin_y + 1)

    x_idx = np.arange(min_x, max_x)
    y_idx = np.arange(min_y, max_y)

    return x_idx, y_idx