import glob
import os
import re
import nibabel as nib
import numpy as np


def natural_sort_key(s):
    return [
        int(text) if text.isdigit() else text.lower()
        for text in re.split(r'(\d+)', s)
    ]


def compute_rms(data, num_echoes=None):
    if num_echoes is not None:
        selected_data = data[..., :num_echoes]
    else:
        selected_data = data

    return np.sqrt(np.mean(np.square(selected_data), axis=-1))

input_dir = '/home/rhiannon.b/Downloads/ReconData/init_zeros_sub9021_lumbar/SOS_bin'
''
output_dir = os.path.join(input_dir, 'combined')
os.makedirs(output_dir, exist_ok=True)

all_files = glob.glob(os.path.join(input_dir, '**', '*.nii*'), recursive=True)
all_files = sorted(all_files, key=natural_sort_key)

print(f'Total slices found: {len(all_files)}')  # 

if len(all_files) == 0:
    raise FileNotFoundError(
        f'No NIfTI files found in {input_dir}. Check path and extensions.'
    )


ref_nii = nib.load(all_files[0])
ref_affine = ref_nii.affine
ref_header = ref_nii.header

first_data = ref_nii.get_fdata()
spatial_shape = first_data.shape[:2] 
num_echoes = first_data.shape[-1]  
total_slices = len(all_files) 


volume_data = np.zeros(
    (*spatial_shape, total_slices, num_echoes), dtype=first_data.dtype
)

for slice_idx, filepath in enumerate(all_files):
    nii = nib.load(filepath)
    volume_data[:, :, slice_idx, :] = nii.get_fdata()

for echo_idx in range(num_echoes):
    echo_volume = volume_data[:, :, :, echo_idx]
    echo_nii = nib.Nifti1Image(echo_volume, ref_affine, ref_header)

    out_path = os.path.join(output_dir, f'echo_{echo_idx + 1}.nii.gz')
    nib.save(echo_nii, out_path)
    print(f'Saved: {out_path}')

rms_1to3_volume = compute_rms(volume_data, num_echoes=3)
rms_1to3_nii = nib.Nifti1Image(rms_1to3_volume, ref_affine, ref_header)
out_path_1to3 = os.path.join(output_dir, 'echo_rms_1to3.nii.gz')
nib.save(rms_1to3_nii, out_path_1to3)
print(f'Saved: {out_path_1to3}')

rms_all_volume = compute_rms(volume_data)
rms_all_nii = nib.Nifti1Image(rms_all_volume, ref_affine, ref_header)
out_path_all = os.path.join(output_dir, 'echo_rms_1to4.nii.gz')
nib.save(rms_all_nii, out_path_all)
print(f'Saved: {out_path_all}')