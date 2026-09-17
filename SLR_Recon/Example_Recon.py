#!/usr/bin/env python
# coding: utf-8

# Joint multi-echo navigator/respiratory resolved structured low-rank (SLR) reconstruction
#
# Mark Chiew (mark.chiew@utoronto.ca)



# Imports
import numpy as np

# Try importing CuPy and checking for an active CUDA GPU
HAS_GPU = False

try:
    import cupy as cp
    # Verify that an active NVIDIA device is accessible
    _ = cp.cuda.Device(0).compute_capability
    HAS_GPU = True
    import gpuSLR
    print(f"NVIDIA GPU detected ({cp.cuda.runtime.getDeviceProperties(0)['name'].decode()}). Running on GPU.")
except ImportError:
    print("CuPy is not installed. Falling back to CPU.")
    import SLR
except Exception as e:
    print(f"GPU initialization failed ({e}). Falling back to CPU.")
    import SLR

import argparse
import numpy as np
from matplotlib import pyplot as plt
import nibabel as nib
import twixtools # for reading/loading raw twix data
import sklearn   # used only for k-means 
import grappa    # grappa for computing initialization
import get_spinal_cord_crop_indices
import time
from scipy.linalg import hankel

print("running new recon")

# Define some helper functions

# handles centric k-space with shifting, along arbitrary dimensions
def fftdim(x, dims=None):
    return np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(x), axes=dims, norm="ortho"))

def ifftdim(x, dims=None):
    return np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(x), axes=dims, norm="ortho"))

# root-sum-of-squares combination
def sos(x, axis=-1):
    return np.sqrt(np.sum(np.abs(x)**2, axis=axis))

def plot_navigator_clusters(navigator_data, cluster_indices, slice_idx, nbins):
    flat_com = navigator_data.reshape(-1, navigator_data.shape[-1])
    flat_labels = cluster_indices.flatten()

    phases = np.angle(flat_com)
    unit_vectors = np.exp(1j * phases)
    mean_unit_vectors = np.mean(unit_vectors, axis=-1)
    circ_mean_phases = np.angle(mean_unit_vectors)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw = {'height_ratios': [2, 1]})
    cmap = plt.get_cmap('tab10', nbins)
    line_indices = np.arange(len(flat_labels))

    scatter = ax1.scatter(line_indices, circ_mean_phases, c=flat_labels, cmap=cmap, s=15, alpha=0.8, edgecolors='none')
    ax1.set_title(f"Navigator Phase Clustering (Slice {slice_idx})")
    ax1.set_xlabel("Line Index")
    ax1.set_ylabel("Circular Mean Phase (radians)")
    ax1.set_yticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
    ax1.set_yticklabels([r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$'])
    ax1.grid(True, alpha=0.3)
    cbar = plt.colorbar(scatter, ax=ax1, ticks=np.arange(nbins))
    cbar.set_label("Cluster Index", rotation=270, labelpad=15)

    for b in range(nbins):
        mask = flat_labels == b
        if np.any(mask):
            ax2.hist(
                circ_mean_phases[mask],
                bins=25,
                range=(-np.pi, np.pi),
                alpha=0.5,
                label=f'Bin {b}',
                color=cmap(b / max(1, nbins - 1)),
            )

    ax2.set_xlabel('Circular Mean Phase (rad)')
    ax2.set_ylabel('Line Count')
    ax2.set_xticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
    ax2.set_xticklabels([r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$'])
    ax2.grid(True, linestyle='--', alpha=0.5)
    ax2.legend(loc='upper right', fontsize=8)

    plt.tight_layout()

    # Save image and close figure to free memory
    filename = f'slice_{slice_idx}_kmeans_navigator_phase_angle.png'
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved navigator plot to {filename}')



# define arguments 
parser = argparse.ArgumentParser(description='Reconstruct a selected slice from the GRE data.')
slice_group = parser.add_mutually_exclusive_group(required=True)
slice_group.add_argument(
    "--slc", 
    type=int, 
    help="index of the single slice to reconstruct"
)
slice_group.add_argument(
    "--all-slices", 
    action="store_true", 
    help="reconstruct all slices in the volume"
)
parser.add_argument('--r', type=int, default=150, help='rank parameter for the SLR reconstruction')
parser.add_argument('--nbins', type=int, default=4, help='number of navigator bins for the reconstruction')
parser.add_argument('--manual_crop', action="store_false", help='use sct segmentation cropping instead of hardcoded indices')
parser.add_argument('--i', type=str, required=True, help='.dat file containing the raw kspace data to be reconstructed')
parser.add_argument('--niters', type=int, default=100, help='number of iterations for the reconstruction')
parser.add_argument('--mode', type=str, default=None, help='reconstruction mode (hard or soft). default is hard.')
parser.add_argument('--lam', type=float, default=5.7, help='lambda parameter for soft thresholding. default is 5.7.')
parser.add_argument('--plot_nav', type=bool, default=False, help='reconstruction mode (hard or soft). default is hard.')
args = parser.parse_args()



# Load data

# Map twix file
print("loading in and whitening data")
map = twixtools.map_twix(args.i)

# Image data
map[-1]['image'].flags['zf_missing_lines'] = True
map[-1]['image'].flags['remove_os'] = True
img = map[-1]['image'][:].squeeze()

# Get image dimensions
nrep, neco, nslc, ny, nc, nx = img.shape

hdr = map[-1]['hdr']

slice0 = hdr['MeasYaps']['sSliceArray']['asSlice'][0]
fov_read = float(slice0.get('dReadoutFOV', 192.0))
fov_phase = float(slice0.get('dPhaseFOV', 192.0))
thickness = float(slice0.get('dThickness', 3.0))
gap = float(slice0.get('dGap', 0.0))

dx = fov_read / nx
dy = fov_phase / ny
dz = thickness + gap

affine = np.array([
    [dx, 0, 0, -(nx * dx) / 2],
    [0, dy, 0, -(ny * dy) / 2],
    [0, 0, dz, -(nslc * dz) / 2],
    [0, 0, 0, 1.0]
], dtype=np.float32)

# Get acceleration factor (assume integer)
R = int(map[-1]['hdr']['Config']['AFLin'])

# Navigator
map[-1]['phasestab'].flags['zf_missing_lines'] = True
map[-1]['phasestab'].flags['remove_os'] = True
img_nav = map[-1]['phasestab'][:].squeeze()

# Reference data
map[-1]['refscan'].flags['skip_empty_lead'] = True
map[-1]['refscan'].flags['remove_os'] = True
ref = map[-1]['refscan'][:].squeeze()
nref = ref.shape[-3]

# Reference navigators
map[-1]['ref_ps'].flags['skip_empty_lead'] = True
map[-1]['ref_ps'].flags['remove_os'] = True
ref_nav = map[-1]['ref_ps'][:].squeeze()

# place slices in anatomical order instead of order in which they were acquired
hdr_slice_series = map[-1]['hdr']['MeasYaps']['sSliceArray']
slc_arr = map[-1]['hdr']['MeasYaps']['sSliceArray']['asSlice']

if 'ucMode' in hdr_slice_series or 'aInSlice' in hdr_slice_series:
    pass

meas = map[-1]

slice_positions = []
for i in range(nslc):
    tra = slc_arr[i].get('sPosition', {}).get('dTra', 0.0)
    slice_positions.append(tra)

# Look for Siemens Interleaved Mode flag (ucMode: 0x1 = Ascending, 0x2 = Descending, 0x4 = Interleaved)
mode = hdr_slice_series.get('ucMode', None)

if mode == 4 or mode == '0x4' or (nslc > 1 and sort_idx[0] == 0): 
    # Interleaved slice ordering (0, 2, 4... 1, 3, 5...)
    # We construct the chronological-to-anatomical index map
    even_slices = list(range(0, nslc, 2))
    odd_slices = list(range(1, nslc, 2))
    acq_to_anat_map = np.array(even_slices + odd_slices)
    
    sort_idx = np.argsort(acq_to_anat_map)
else:
    sort_idx = np.array(range(nslc))

print("Corrected anatomical sort indices:", sort_idx)
img = img[:, :, sort_idx, ...]
img_nav = img_nav[:, sort_idx, ...]
ref = ref[:, :, sort_idx, ...]
ref_nav = ref_nav

# Compute and apply noise pre-whitening transform
# This is optional - but if you don't apply the tensordot transformation, you need to manually move the channel dimension to the end
W = np.linalg.cholesky(np.linalg.inv(np.cov(np.reshape(map[0]['noise'][:].squeeze().transpose((1,0,2)),(nc,-1)))))

img = np.tensordot(img, W, axes=((-2,),(0)))
print(img.shape)
ref = np.tensordot(ref, W, axes=((-2,),(0)))
img_nav = np.tensordot(img_nav, W, axes=((-2,),(0)))
ref_nav = np.tensordot(ref_nav, W, axes=((-2,),(0)))



# Slice-specific prep 
# The data and binning are done on a slice-specific basis. Needs to be run once per slice

# Define parameters for reconstruction
# echoes to reconstruct
eco = np.arange(neco)

# repetitions to reconstruct
rep = np.arange(nrep)

# number of bins to partition data into, based on navigator k-means clustering
# increasing this number increases the number of resolved dynamic states, but also increases computation time and memory
# on CPU, I would recommend nbins ≤ 8, anything beyond that gets pretty slow
# on GPU, I've tested up to nbins = 16 and it works reasonably fast, probably diminishing returns
nbins = args.nbins

# need grappa initial guess to get indices
# GRAPPA recon for initialization
# we average across repetitions for this
# for an acceleration factor 2, we use a grappa kernel of size (3,2)
# Initialize container for all slices: (nc, nx, ny, nslc, neco)
init_all = np.zeros((nslc, nc, nx, ny, neco), dtype=np.complex64)
for s in range(nslc):
    for i in range(neco):
        init_all[s,:,:,:,i] = grappa.grappa(np.mean(img[:,i,s,:,:,:], axis=0).transpose((2,1,0)), np.mean(ref[:,i,s,:,:,:], axis=0).transpose((2,1,0)), (1,2), (3,2))


img_space = np.fft.fftshift(
    np.fft.ifft2(
        np.fft.ifftshift(init_all, axes=(2, 3)), 
        axes=(2, 3)
    ), 
    axes=(2, 3)
)
img_3d_full = sos(img_space, axis=1) # sos combine the channel dimension
img_3d_full = np.abs(img_3d_full[:,:,:,1]).transpose((1,2,0))  # just use first echo for segmentation
img_3d_full = (img_3d_full / np.percentile(img_3d_full, 99.9)) * 1000.0 # normalize for sct

# 3. Clip negative or extreme out-of-range values
img_3d_full = np.clip(img_3d_full, 0, None)


# select RO indices near the spinal cord, as we only care about that region
sct_crop = args.manual_crop
print("cropping around spinal cord")
if sct_crop:
    x_idx, y_idx = get_spinal_cord_crop_indices.get_indices(
    img_3d_full,
    affine,
    "test",
    args.r,
)

    # Use x_idx for navigator selection
    sc_idx = x_idx
    print("cropping indices:", sc_idx)
# maybe should just get rid of this option?
else:
    sc_idx = (164, 220) # these indices were used for the original Zurich data that we sent mark


# slice to reconstruct
if args.all_slices:
    slices_to_process = range(nslc)
else:
    slices_to_process = [args.slc]

for slc in slices_to_process:
    init = init_all[slc, :, :, :, :]
    print("processing slice:", slc)
    # reference everything relative to the first line, and average the relative signals across coil channels
    # concatenate navigators, and inverse FFT navigators along RO dimension
    # use only sampled lines, to avoid an extra trivial cluster from the empty lines
    nav = ifftdim(np.concatenate((img_nav[:,slc,::R,:,:], ref_nav[:,slc,:,:,:]), axis=1), dims=(-2,))

    if args.plot_nav == True:
        print('plotting navigators...')
        line_idx = 0  # Index of the repetition/line you want to inspect

        # 1. Navigator line already in image domain (from ifftdim along RO dimension)
        # 1. Navigator line (averaged across channels)
        nav_profile = np.angle(
            np.mean(nav[line_idx, line_idx, sc_idx, :], axis=-1)
        )

        # Extract matching image k-space line
        img_line_kspace = img[line_idx, 0, slc, 0, sc_idx, 0]
        img_line_imdom = np.angle(
            np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(img_line_kspace)))
        )

        # 3. Plot both magnitude profiles
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.5))

        ax1.plot(nav_profile, color='crimson')
        ax1.set_title(f'Navigator Line {line_idx} (Image Domain)')
        ax1.set_xlabel('Readout (RO) Index')
        ax1.set_ylabel('Phase')
        ax1.grid(True, alpha=0.3)

        ax2.plot(img_line_imdom, color='navy')
        ax2.set_title(f'Image Line {line_idx} (Image Domain)')
        ax2.set_xlabel('Readout (RO) Index')
        ax2.set_ylabel('Phase')
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(
            f'slice_{slc}_line_{line_idx}_nav_vs_img.png',
            dpi=150,
            bbox_inches='tight',
        )
        plt.close(fig)
        print(
            f'Saved line profile comparison to slice_{slc}_line_{line_idx}_nav_vs_img.png'
        )

    tmp = nav[:, :, sc_idx, :]
    tmp_complex = np.squeeze(np.mean(tmp*np.conj(tmp[:,[0],:,:]), axis=-1))

    # because sk-learn k-means requires real-valued input, I've concatenated the real and imaginary parts of the navigator
    tmp = np.concatenate((np.real(tmp_complex), np.imag(tmp_complex)), axis=-1)

    # alternatively, you could try extracting just the phase of the navigator
    #tmp = np.angle(tmp_complex)
    # get k-means cluster indices, with nbins clusters
    idx = sklearn.cluster.KMeans(n_clusters=nbins, random_state=42).fit(tmp.reshape((-1,tmp.shape[-1]))).labels_.reshape((nrep,-1))
    plot_navigator_clusters(tmp_complex, idx, slc, nbins)
    #Prep binned data and initialization

    # sort data into new bin dimension using k-means indices
    # the data across all repetitions is being used here, as well as the reference data
    # the cnt array just keeps track in case the same line appears in the same bin across repetitions
    # if this happens, we simply average the lines
    dat = np.zeros((nx, ny, nbins, neco, nc), dtype='complex64')

    cnt = np.zeros((nx, ny, nbins, neco, nc))
    ref_offset = ny//2 - nref//2
    for i in range(nrep): 
        for j in range(ny//R):
            dat[:, R*j, idx[i,j], :, :] += img[i, :, slc, R*j, :, :].transpose((1,0,2))
            cnt[:, R*j, idx[i,j], :, :] += 1
        for k in range(nref):
            dat[:, ref_offset+k, idx[i,ny//R+k], :, :] += ref[i, :, slc, k, :, :].transpose((1,0,2))
            cnt[:, ref_offset+k, idx[i,ny//R+k], :, :] += 1    

    dat[dat!=0] = dat[dat!=0]/cnt[dat!=0]

    # for the initalization, we just copy the GRAPPA recon across the bin dimension
    init  = np.tile(init[:,:,:,None,:].transpose((1,2,3,4,0)),(1,1,nbins,1,1))

    # reshape input data and initialization to combine bin, eco and channel dimensions
    dat = dat.reshape((nx, ny, -1))
    init = init.reshape((nx, ny, -1))


    #Crop RO dimension
    # choose RO indices to keep
    # this is not strictly necessary, but I recommend it, particularly if you used the cropped navigator for k-means clustering
    # image quality near the spinal cord will be better, because the SLR reconstruction doesn't have to "fit" the entire FOV all at once
    # it is possible to do this, of course, but it requires a bit more tweaking of hyperparameters (kernel size, rank, etc.)
    # if you do want a full FOV image, I would actually recommend trying to generate it with a series R0 cropped reconstructions, and combining afterwards
    # an example of this is provided in Full_FOV_Recon.ipynb
    xidx = sc_idx
    nx_crop = len(xidx)

    # ifft to x-dimension, crop, the fft back to kx
    # doing this both for the prepared data and the initialization
    dat = fftdim(ifftdim(dat, dims=(0,))[xidx, :, :], dims=(0,))
    init = fftdim(ifftdim(init, dims=(0,))[xidx, :, :], dims=(0,))


    # ## Reconstruction
    # There is a CPU version (SLR) and a GPU version (gpuSLR). Otherwise reconstruction function calls are very similar
    
    # The reconstruction uses an alternating direction method of multipliers (ADMM) optimization to solve the structured low-rank constrained reconstruction.
    
    # The third input parameter is the type of structured low-rank matrix formulation to use. There are several options:\
    # `c_matrix`: the most basic, straightforward phase smoothness and limited image support constraint\
    # `s_matrix`: everything c_matrix does, but additionally exploits some conjugate symmetry properties\
    # `vcc_matrix`: similar to s_matrix, but formulated differently, using the virtual conjugate coil framework

    # SLR reconstruction

    # set rank parameter
    # this is a bit tricky to tune - lower numbers will result in greater regularization
    # set too low, signal loss in the output can results
    # set too high, nothing really happens
    # also, this number interacts with kernel size and type of SLR matrix. A larger kernel may require a larger r value to prevent over-regularization
    r = args.r
    lam = args.lam

    # set number of iterations
    niters = args.niters

    # slr kernel size
    kernel = (5,5)
    # example gpu reconstruction using the c_matrix
    start_time = time.perf_counter()
    if HAS_GPU:
        if args.mode is not None:
            out = gpuSLR.ADMM(dat, gpuSLR.c_matrix, kernel, lam, niters=niters, init=init, mode=args.mode)
        else:   
            print("running SLR recon with zeros as initialization...")
            out = gpuSLR.ADMM(dat,              # input data
                            gpuSLR.c_matrix,    # type of structured low-rank matrix. options are `c_matrix`, `s_matrix` or `vcc_matrix`
                            kernel,             # SLR kernel size
                            r,                  # rank (d
                            niters=niters)      # initialization (defaults to array of zeros)
    else:
        if args.mode is not None:
            out = SLR.ADMM(dat, SLR.c_matrix, kernel, lam, niters=niters, init=init, mode=args.mode)
        else:
            # similar reconstruction using cpu
            out = SLR.ADMM(dat, SLR.c_matrix, kernel, r, niters=niters, init=init)

    elapsed_time = time.perf_counter() - start_time
    print(f'ADMM Converged in {elapsed_time:.2f} seconds ({niters} iterations)')

    # Plot results
    # reshape and ifftdim output 
    # use the reconstructed result, not the initialization, so nbins changes are visible
    mag = ifftdim(out.reshape((nx_crop, ny, nbins, neco, nc)), dims=(0,1))

    for k in range(nbins):
        bin_mag = sos(mag[:, :, k, :, :], axis=-1)
        
        # Save using proper orientation affine
        nib.save(nib.Nifti1Image(bin_mag, affine), f'{slc}_bin_{k}_result_{nbins}_{r}_{niters}.nii.gz')

    # typically for magnitude images, we would sos-combine the bin and channel dimensions
    # this is not necessary, you can keep the bin-dimension uncombined and do something else if you like
    # the bin dimension resolves the different navigator states
    mag = sos(mag.transpose((0,1,3,2,4)).reshape((nx_crop, ny, neco, -1))) 
    y_vis = np.arange(96,224)

    # plot all magnitude of all echoes
    _, ax = plt.subplots(1, neco, figsize=(12,12*(2/neco))) 
    for i in range(neco):
        ax[i].imshow(np.rot90(mag[:,y_vis,i]), vmin=0, vmax=np.max(mag)*.8, cmap='gray')
        ax[i].set_title(f'Recon Echo {i}')

    # save results as nifti
    final = nib.Nifti1Image(mag, np.eye(4))
    nib.save(final, f'{slc}_recon_result_{nbins}_{r}_{niters}.nii.gz')

print("reconstruction finished!")