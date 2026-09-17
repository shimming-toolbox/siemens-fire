from t2star_diagnostics import load_echo_niftis, load_cord_mask_nifti, run_t2star_diagnostics
from t2star_multiple import plot_multi_echo_3d_comparison
import nibabel as nib


# these are true for the Zurich dataset of 10 subjects (4 with only thoracic images)
TEs = [6.9, 10.92, 14.94, 18.96]  # echo times in seconds
'''
combined_mag, TEs = load_echo_niftis(
    ["/home/rhiannon.b/Downloads/rank1000_iters1000_dataset3/combined/echo_1.nii.gz", 
     "/home/rhiannon.b/Downloads/rank1000_iters1000_dataset3/combined/echo_2.nii.gz", 
     "/home/rhiannon.b/Downloads/rank1000_iters1000_dataset3/combined/echo_3.nii.gz", 
     "/home/rhiannon.b/Downloads/rank1000_iters1000_dataset3/combined/echo_4.nii.gz"], TEs
)
cord_mask = load_cord_mask_nifti("/home/rhiannon.b/Downloads/ReconData/SLR_recon_results/MarkResults_DefaultParams/meas_MID00151_FID35300_gre_spine/combined/dataset3_thoracic_seg.nii.gz")

results_r150 = run_t2star_diagnostics(combined_mag, cord_mask, TEs, rank=150, out_dir="t2star_dataset3_regular_rank")
'''


reconstruction_files = {
    'SLR joint-echoes': [
        '/home/rhiannon.b/Downloads/ReconData/SLR_recon_results/sub9021/thoracic/thoracic_4bins/combined/echo_1.nii.gz',
        '/home/rhiannon.b/Downloads/ReconData/SLR_recon_results/sub9021/thoracic/thoracic_4bins/combined/echo_2.nii.gz',
        '/home/rhiannon.b/Downloads/ReconData/SLR_recon_results/sub9021/thoracic/thoracic_4bins/combined/echo_3.nii.gz',
        '/home/rhiannon.b/Downloads/ReconData/SLR_recon_results/sub9021/thoracic/thoracic_4bins/combined/echo_4.nii.gz',
    ],
    'SLR separate-echoes': [
        '/home/rhiannon.b/Downloads/seperate_echoes_meas00061/SOS_bin/combined/echo_1.nii.gz',
        '/home/rhiannon.b/Downloads/seperate_echoes_meas00061/SOS_bin/combined/echo_2.nii.gz',
        '/home/rhiannon.b/Downloads/seperate_echoes_meas00061/SOS_bin/combined/echo_3.nii.gz',
        '/home/rhiannon.b/Downloads/seperate_echoes_meas00061/SOS_bin/combined/echo_4.nii.gz',
    ],
}

mask_nifti = nib.load('/home/rhiannon.b/Downloads/ReconData/SLR_recon_results/sub9021/thoracic/sub9021_seg.nii.gz')

# Extract data as a boolean matrix (True inside the cord, False outside)
cord_mask = mask_nifti.get_fdata() > 0

# 3. Call the plotter
plot_multi_echo_3d_comparison(
    recon_dict=reconstruction_files,
    TEs=TEs,
    mask_3d =cord_mask,
    save_path="multi_rank_comparison.png"
)