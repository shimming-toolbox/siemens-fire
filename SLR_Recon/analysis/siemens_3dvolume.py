import os
import re
import glob
import numpy as np
import nibabel as nib


base_dir = "/home/rhiannon.b/Downloads/ReconData/sub9025/siemensniftis/tmp_dcm2bids/sub-9025/raw"  
rep_folders = ["rep1", "rep2", "rep3"]
output_dir = os.path.join(base_dir, "rep_averaged")
os.makedirs(output_dir, exist_ok=True)

def get_echo_and_instance(filename):
    match = re.search(r'_e(\d+)_i(\d+)\.nii\.gz$', filename)
    if match:
        return int(match.group(1)), int(match.group(2))
    return None

file_groups = {}

for rep_idx, folder in enumerate(rep_folders):
    folder_path = os.path.join(base_dir, folder)
    for fpath in glob.glob(os.path.join(folder_path, "*.nii.gz")):
        key = get_echo_and_instance(os.path.basename(fpath))
        if key:
            if key not in file_groups:
                file_groups[key] = [None] * len(rep_folders)
            file_groups[key][rep_idx] = fpath

print("--- Step 1: Repetition Averaging ---")
averaged_files = {}

for (echo, instance), paths in sorted(file_groups.items()):
    if any(p is None for p in paths):
        continue

    imgs = [nib.load(p) for p in paths]
    mean_data = np.mean([img.get_fdata() for img in imgs], axis=0)
    
    ref_img = imgs[0]
    out_path = os.path.join(output_dir, f"avg_e{echo}_i{instance:05d}.nii.gz")
    nib.save(nib.Nifti1Image(mean_data, ref_img.affine, ref_img.header), out_path)

    averaged_files[(echo, instance)] = {
        'data': mean_data,  # Expected shape: (191, 191, 2) or (2, 191, 191)
        'affine': ref_img.affine,
        'header': ref_img.header
    }

print(f"Averaged {len(averaged_files)} repetition files.")

print("\n--- Step 2: Assembling (191, 191, 14) Volumes & RMS ---")

echoes = sorted(list(set(e for e, i in averaged_files.keys())))

if len(echoes) >= 4:
    echo_3d_vols = []

    for e in range(1, 5):
        instances = sorted([inst for (echo_num, inst) in averaged_files.keys() if echo_num == e])
        
        echo_slices_14 = []
        for inst in instances:
            file_data = averaged_files[(e, inst)]['data']
            
            # Robust extraction of sub-slices to ensure exact (191, 191) 2D shape
            if file_data.ndim == 3 and file_data.shape[2] == 2:
                s0 = file_data[:, :, 0]
                s1 = file_data[:, :, 1]
            elif file_data.ndim == 3 and file_data.shape[0] == 2:
                s0 = file_data[0, :, :]
                s1 = file_data[1, :, :]
            else:
                raise ValueError(f"Unexpected array shape for Echo {e}, Instance {inst}: {file_data.shape}")

            echo_slices_14.append(s0)
            echo_slices_14.append(s1)
        vol_3d = np.stack(echo_slices_14, axis=-1)
        echo_3d_vols.append(vol_3d)
        print(f"Echo {e} assembled 3D Volume Shape: {vol_3d.shape}  (X, Y, Z)")

    multi_echo_4d = np.stack(echo_3d_vols, axis=-1)
    print(f"\nFull 4D Multi-Echo Array Shape: {multi_echo_4d.shape}  (X, Y, Z, Echo)")

    rms_4echo = np.sqrt(np.mean(np.square(multi_echo_4d), axis=-1))
    rms_3echo = np.sqrt(np.mean(np.square(multi_echo_4d[..., :3]), axis=-1))

    print(f"\nFinal RMS 4-Echo Volume Shape: {rms_4echo.shape}  (X, Y, Z)")
    print(f"Final RMS 3-Echo Volume Shape: {rms_3echo.shape}  (X, Y, Z)")

    first_key = (1, min(inst for e, inst in averaged_files.keys() if e == 1))
    ref_affine = averaged_files[first_key]['affine'].copy()


    path_rms4 = os.path.join(output_dir, "rms_4echoes_3D.nii.gz")
    path_rms3 = os.path.join(output_dir, "rms_3echoes_3D.nii.gz")

    nib.save(nib.Nifti1Image(rms_4echo, ref_affine), path_rms4)
    nib.save(nib.Nifti1Image(rms_3echo, ref_affine), path_rms3)

    print(f"\nSaved 3D RMS NIfTI files:")
    print(f" -> {path_rms4}")
    print(f" -> {path_rms3}")
else:
    print(f"Error: Found {len(echoes)} echoes, need 4.")