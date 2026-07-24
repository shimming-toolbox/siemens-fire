#!/usr/bin/env python
# coding: utf-8

# # Joint multi-echo navigator/respiratory resolved structured low-rank (SLR) reconstruction
# 
# Mark Chiew (mark.chiew@utoronto.ca)

# ---
# ## Setup & Load Data
# ---
# This only needs to be done once per twix file

# Imports
import argparse
import numpy as np
from matplotlib import pyplot as plt
import nibabel as nib
import twixtools # for reading/loading raw twix data
import sklearn   # used only for k-means 
#import gpuSLR    # GPU version of SLR
import SLR       # Structured low-rank methods for joint recon
import grappa    # grappa for computing initialization
import subprocess
import os
# Define some helper functions

# handles centric k-space with shifting, along arbitrary dimensions
def fftdim(x, dims=None):
    return np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(x), axes=dims, norm="ortho"))

def ifftdim(x, dims=None):
    return np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(x), axes=dims, norm="ortho"))

# root-sum-of-squares combination
def sos(x, axis=-1):
    return np.sqrt(np.sum(np.abs(x)**2, axis=axis))

# get spinal cord indices along the readout dimension, using SCT deep learning segmentation
def get_spinal_cord_crop_indices(
    img_3d, 
    affine, 
    contrast="t2", 
    margin_x=20,  
    margin_y=5,   
    x_axis=0, 
    y_axis=1
):
    temp_img_path = "temp_sc_recon.nii.gz"
    temp_seg_path = "temp_sc_seg.nii.gz"
    
    nifti_img = nib.Nifti1Image(img_3d.astype(np.float32), affine)
    nib.save(nifti_img, temp_img_path)
    
    subprocess.run([
        "/home/rhiannon.b/sct/bin/sct_deepseg_sc",
        "-i", temp_img_path,
        "-c", contrast,
        "-o", temp_seg_path
    ], check=True, stdout=subprocess.DEVNULL)
    
    seg_mask = nib.load(temp_seg_path).get_fdata()

    for f in [temp_img_path, temp_seg_path]:
        if os.path.exists(f): os.remove(f)

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
# ### Load and whiten data


# Map twix file
map = twixtools.map_twix('meas_MID00151_FID35300_gre_spine.dat')

# Image data
map[-1]['image'].flags['zf_missing_lines'] = True
map[-1]['image'].flags['remove_os'] = True
img = map[-1]['image'][:].squeeze()

# Get image dimensions
nrep, neco, nslc, ny, nc, nx = img.shape

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


# Compute and apply noise pre-whitening transform
# This is optional - but if you don't apply the tensordot transformation, you need to manually move the channel dimension to the end
W = np.linalg.cholesky(np.linalg.inv(np.cov(np.reshape(map[0]['noise'][:].squeeze().transpose((1,0,2)),(nc,-1)))))

img = np.tensordot(img, W, axes=((-2,),(0)))
ref = np.tensordot(ref, W, axes=((-2,),(0)))
img_nav = np.tensordot(img_nav, W, axes=((-2,),(0)))
ref_nav = np.tensordot(ref_nav, W, axes=((-2,),(0)))


# ---
# ## Slice-specific prep 
# ---
# 
# The data and binning are done on a slice-specific basis. Needs to be run once per slice

# ### Define parameters for reconstruction

parser = argparse.ArgumentParser(description='Reconstruct a selected slice from the GRE data.')
parser.add_argument('--slc', type=int, default=2, help='slice index to reconstruct')
parser.add_argument('--r', type=int, default=150, help='rank parameter for the SLR reconstruction')
parser.add_argument('--nbins', type=int, default=4, help='number of navigator bins for the reconstruction')
args = parser.parse_args()

# slice to reconstruct
slc = args.slc
if not 0 <= slc < nslc:
    raise ValueError(f'slc={slc} is out of range for the available {nslc} slices')

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
init_all = np.zeros((nslc, nc, nx, ny, neco), dtype=np.complex64)

for s in range(nslc):
    for i in eco:
        # GRAPPA recon per slice
        init_all[s, :, :, :, i] = grappa.grappa(
            np.mean(img[:, i, s, :, :, :], axis=0).transpose((2, 1, 0)),
            np.mean(ref[:, i, s, :, :, :], axis=0).transpose((2, 1, 0)),
            (1, 2), (3, 2)
        )
img_3d_full = np.mean(np.abs(init_all), axis=(1, 4)).transpose((1, 2, 0)) # shape (nx, ny, nslc)

# ### Navigator binning with k-means


# concatenate navigators, and inverse FFT navigators along RO dimension
# use only sampled lines, to avoid an extra trivial cluster from the empty lines
nav = ifftdim(np.concatenate((img_nav[:,slc,::R,:,:], ref_nav[:,slc,:,:,:]), axis=1), dims=(-2,))

# select RO indices near the spinal cord, as we only care about that region
affine = np.eye(4)
x_idx, y_idx = get_spinal_cord_crop_indices(
    img_3d_full, 
    affine, 
    contrast="t2", 
    margin_x=20,  
    margin_y=5,   
    x_axis=0, 
    y_axis=1
)

# Use x_idx for navigator selection
sc_idx = x_idx
init = init_all[slc]

# reference everything relative to the first line, and average the relative signals across coil channels
tmp = nav[:, :, sc_idx, :]
tmp = np.squeeze(np.mean(tmp*np.conj(tmp[:,[0],:,:]), axis=-1))

# because sk-learn k-means requires real-valued input, I've concatenated the real and imaginary parts of the navigator
tmp = np.concatenate((np.real(tmp), np.imag(tmp)), axis=-1)

# alternatively, you could try extracting just the phase of the navigator
#tmp = np.angle(tmp)

# get k-means cluster indices, with nbins clusters
idx = sklearn.cluster.KMeans(n_clusters=nbins).fit(tmp.reshape((-1,tmp.shape[-1]))).labels_.reshape((nrep,-1))


# ### Prep binned data and initialization



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


# ### Crop RO dimension


# choose RO indices to keep
# this is not strictly necessary, but I recommend it, particularly if you used the cropped navigator for k-means clustering
# image quality near the spinal cord will be better, because the SLR reconstruction doesn't have to "fit" the entire FOV all at once
# it is possible to do this, of course, but it requires a bit more tweaking of hyperparameters (kernel size, rank, etc.)
# if you do want a full FOV image, I would actually recommend trying to generate it with a series R0 cropped reconstructions, and combining afterwards
# an example of this is provided in Full_FOV_Recon.ipynb
xidx = sc_idx
nx = len(xidx)

# ifft to x-dimension, crop, the fft back to kx
# doing this both for the prepared data and the initialization
dat = fftdim(ifftdim(dat, dims=(0,))[xidx, :, :], dims=(0,))
init = fftdim(ifftdim(init, dims=(0,))[xidx, :, :], dims=(0,))


# ---
# ## Reconstruction
# ---
# There is a CPU version (SLR) and a GPU version (gpuSLR). Otherwise reconstruction function calls are very similar
# 
# The reconstruction uses an alternating direction method of multipliers (ADMM) optimization to solve the structured low-rank constrained reconstruction.
# 
# The third input parameter is the type of structured low-rank matrix formulation to use. There are several options:\
# `c_matrix`: the most basic, straightforward phase smoothness and limited image support constraint\
# `s_matrix`: everything c_matrix does, but additionally exploits some conjugate symmetry properties\
# `vcc_matrix`: similar to s_matrix, but formulated differently, using the virtual conjugate coil framework

# ### SLR reconstruction



# set rank parameter
# this is a bit tricky to tune - lower numbers will result in greater regularization
# set too low, signal loss in the output can results
# set too high, nothing really happens
# also, this number interacts with kernel size and type of SLR matrix. A larger kernel may require a larger r value to prevent over-regularization
r = args.r

# set number of iterations
niters = 100

# slr kernel size
kernel = (5,5)

# example gpu reconstruction using the c_matrix
'''
out = gpuSLR.ADMM(dat,                # input data
                 gpuSLR.c_matrix,    # type of structured low-rank matrix. options are `c_matrix`, `s_matrix` or `vcc_matrix`
                  kernel,             # SLR kernel size
                  r,                  # rank (d
                  niters=niters,      # number of iterations (default 100)
                  init=init)          # initialization (defaults to array of zeros)
'''
# similar reconstruction using cpu
out = SLR.ADMM(dat, SLR.c_matrix, kernel, r, niters=niters, init=init)

# example gpu reconstruction with no initialization
# this will work, but requires more iterations to converge
# out = gpuSLR.ADMM(dat, gpuSLR.c_matrix, kernel, r, niters=niters*10)


# ### Plot results

# reshape and ifftdim output 
# use the reconstructed result, not the initialization, so nbins changes are visible
mag = ifftdim(out.reshape((nx, ny, nbins, neco, nc)), dims=(0,1))

# typically for magnitude images, we would sos-combine the bin and channel dimensions
# this is not necessary, you can keep the bin-dimension uncombined and do something else if you like
# the bin dimension resolves the different navigator states
mag = sos(mag.transpose((0,1,3,2,4)).reshape((nx, ny, neco, -1))) #maybe need to change this to neco

y_vis = np.arange(96,224)

# plot all magnitude of all echoes
_, ax = plt.subplots(1, neco, figsize=(12,12*(2/neco))) 
for i in range(neco):
    ax[i].imshow(np.rot90(mag[:,y_vis,i]), vmin=0, vmax=np.max(mag)*.8, cmap='gray')
    ax[i].set_title(f'Recon Echo {i}')


# save results as nifti
img = nib.Nifti1Image(mag, np.eye(4))
nib.save(img, f'{slc}_gpu_recon_result.nii.gz')
