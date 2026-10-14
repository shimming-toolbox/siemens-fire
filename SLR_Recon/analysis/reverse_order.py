import nibabel as nib
import numpy as np

def reverse_slice_order(input_nifti_path: str, output_nifti_path: str, axis: int = 2):
    """
    Reverses the slice order along the specified 3rd dimension (axis 2 by default).
    """
    # 1. Load the NIfTI file
    nifti_img = nib.load(input_nifti_path)
    data = nifti_img.get_fdata()

    # 2. Reverse the array along the slice axis (axis=2)
    reversed_data = np.flip(data, axis=axis)

    # 3. Save as a new NIfTI file with original affine and header
    reversed_img = nib.Nifti1Image(reversed_data, nifti_img.affine, nifti_img.header)
    nib.save(reversed_img, output_nifti_path)
    print(f"Reversed NIfTI saved successfully to: {output_nifti_path}")

# Example Usage:
reverse_slice_order("/home/rhiannon.b/Downloads/ReconData/sub9027/thoracic/niftis/tmp_dcm2bids/sub-9027/normorraw/rep_averaged/rms_4echoes_3D.nii.gz", 
                    "/home/rhiannon.b/Downloads/ReconData/sub9027/thoracic/niftis/tmp_dcm2bids/sub-9027/normorraw/rep_averaged/rms_4echoes_3D_fixed.nii.gz", axis=2)