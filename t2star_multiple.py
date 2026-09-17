import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

# -------------------------------------------------------------------------
# 1. Helper fitting function
# -------------------------------------------------------------------------
def mono_exp_with_noise(TE, S0, T2star, C):
    return S0 * np.exp(-TE / T2star) + C

def fit_decay(TEs, signal):
    """Fits 3-parameter exponential decay to 3D ROI-averaged signal."""
    p0 = [np.max(signal), 25.0, np.min(signal)]
    bounds = ([0, 1.0, 0], [np.inf, 200.0, np.max(signal)])
    try:
        popt, _ = opt.curve_fit(mono_exp_with_noise, TEs, signal, p0=p0, bounds=bounds)
        return popt
    except Exception as e:
        print(f"Fit failed: {e}")
        return [np.nan, np.nan, np.nan]

# -------------------------------------------------------------------------
# 2. Multi-File 3D Comparison Plotter
# -------------------------------------------------------------------------
def plot_multi_echo_3d_comparison(recon_dict, TEs, mask_3d=None, save_path="multi_rank_3d_comparison.png"):
    """
    recon_dict: dict of { 'Label': ['eco_0.nii.gz', 'eco_1.nii.gz', 'eco_2.nii.gz', 'eco_3.nii.gz'] }
    TEs: list or array of echo times in ms
    mask_3d: 3D boolean array of shape (nx, ny, nslices)
    """
    fig, ax = plt.subplots(figsize=(9, 6), dpi=300)
    te_fit = np.linspace(min(TEs), max(TEs), 200)
    colors = plt.cm.Set1(np.linspace(0, 1, len(recon_dict)))

    print("--- 3D Volume Fitting Results ---")

    for (label, echo_files), color in zip(recon_dict.items(), colors):
        echo_volumes = []
        
        # 1. Load each echo file (15 slices)
        for file_path in echo_files:
            img_data = nib.load(file_path).get_fdata()
            
            # If image has a 4th bin dimension: (nx, ny, nbins, nslices) -> average over bins
            if img_data.ndim == 4:
                img_data = np.nanmean(img_data, axis=2)
                
            echo_volumes.append(img_data)  # Shape per echo: (nx, ny, nslices)

        # Stack into a 4D array: (nx, ny, nslices, neco)
        data_4d = np.stack(echo_volumes, axis=-1)

        # 2. Extract ROI mean signal across all 15 slices
        if mask_3d is not None:
            # Ensure mask is boolean and broadcast across echoes
            mask_bool = mask_3d > 0
            
            # Extract all cord voxels across all slices: resulting shape (N_cord_voxels, neco)
            cord_voxels = data_4d[mask_bool, :]
            mean_signal = np.nanmean(cord_voxels, axis=0)  # Shape: (neco,)
        else:
            # Global mean across all 3 spatial dimensions if no mask provided
            mean_signal = np.nanmean(data_4d, axis=(0, 1, 2))

        # 3. Fit exponential model to the high-SNR 3D-averaged signal
        S0, T2star, C = fit_decay(TEs, mean_signal)
        print(f"[{label}] T2* = {T2star:.2f} ms | C = {C:.3f} | S0 = {S0:.1f}")

        # 4. Plot raw ROI data points
        ax.plot(TEs, mean_signal, 'o', color=color, markersize=7, label=f"{label} Data")

        # 5. Plot fitted curve
        if not np.isnan(T2star):
            curve = mono_exp_with_noise(te_fit, S0, T2star, C)
            ax.plot(
                te_fit, 
                curve, 
                '--', 
                color=color, 
                linewidth=2,
                label=f"{label} Fit ($T_2^*={T2star:.1f}$ ms, $C={C:.2f}$)"
            )

    # Formatting
    ax.set_xlabel("Echo Time $TE$ (ms)", fontsize=12, fontweight='bold')
    ax.set_ylabel("Mean ROI Signal Intensity", fontsize=12, fontweight='bold')
    ax.set_title("Spinal Cord $T_2^*$ Decay Across 3D Reconstructions (15 Slices)", fontsize=13, fontweight='bold')
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True, fontsize=10)
    
    plt.tight_layout()
    plt.savefig(save_path, bbox_inches='tight')
    plt.show()
    print(f"Plot saved to: {save_path}")