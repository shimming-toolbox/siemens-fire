import os
import time
import tempfile
import subprocess
import nibabel as nib
import numpy as np


def get_indices(
    img_3d,
    affine,
    id_num="unknown",
    r=150,
    margin_x=10,
    margin_y=5,
    x_axis=0,
    y_axis=1,
    task="spinalcord",
):
    t0 = time.time()
    sct_binary = "/home/ge.polymtl.ca/p124290/spinalcordtoolbox_latest/bin/sct_deepseg"

    with tempfile.TemporaryDirectory(dir="/tmp") as tmpdir:
        temp_img_path = os.path.join(tmpdir, f"{id_num}_{r}_input.nii.gz")
        temp_seg_path = os.path.join(tmpdir, f"{id_num}_{r}_seg.nii.gz")

        nifti_img = nib.Nifti1Image(img_3d.astype(np.float32), affine)
        nib.save(nifti_img, temp_img_path)

        run_env = os.environ.copy()
        run_env.pop("PYTHONPATH", None)
        run_env.pop("LD_LIBRARY_PATH", None)
        run_env["CUDA_VISIBLE_DEVICES"] = "0"

        cmd_args = [sct_binary, task, "-i", temp_img_path, "-o", temp_seg_path]

        t_infer_start = time.time()
        try:
            subprocess.run(
                cmd_args,
                check=True,
                capture_output=True,
                text=True,
                env=run_env
            )
        except subprocess.CalledProcessError as e:
            print("SCT STDERR:\n", e.stderr)
            raise e
        t_infer = time.time() - t_infer_start

        seg_mask = nib.load(temp_seg_path).get_fdata() > 0

    x_profile = np.sum(seg_mask, axis=tuple(i for i in range(seg_mask.ndim) if i != x_axis))
    x_nonzero = np.where(x_profile > 0)[0]

    y_profile = np.sum(seg_mask, axis=tuple(i for i in range(seg_mask.ndim) if i != y_axis))
    y_nonzero = np.where(y_profile > 0)[0]

    if len(x_nonzero) == 0 or len(y_nonzero) == 0:
        raise ValueError(f"sct_deepseg failed to detect the {task}.")

    min_x = max(0, np.min(x_nonzero) - margin_x)
    max_x = min(img_3d.shape[x_axis], np.max(x_nonzero) + margin_x + 1)

    min_y = max(0, np.min(y_nonzero) - margin_y)
    max_y = min(img_3d.shape[y_axis], np.max(y_nonzero) + margin_y + 1)

    print(f"  [SCT sc-crop Breakdown] Inference: {t_infer:.2f}s | Total Function: {time.time()-t0:.2f}s")

    return np.arange(min_x, max_x), np.arange(min_y, max_y)
