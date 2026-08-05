"""
GRAPPA 

Mark Chiew (mark.chiew@utoronto.ca)
Python port of https://github.com/mchiew/grappa-tools/tree/main/grappa

"""
import numpy as np

def grappa(data, calib, R, kernel):
    """
    inputs: 
        data    -   (nc, nx, ny, nz) complex undersampled k-space data
        calib   -   (nc, cx, cy, cz) complex fully-sampled calibration k-space data
        R       -   [1, Ry] or [1, Ry, Rz] acceleration factors
                    requires that the first spatial dimension (x) be fully sampled
        kernel  -   [kx, ky] or [kx, ky, kz] kernel size 

    output:
        recon   -   (nc, nx, ny, nz) complex reconstructed k-space data


    example usage:
        # for a 2D dataset (nc, nx, ny), Ry=2, kernel 3 points in kx, 2 points in ky
        recon = grappa.grappa(data, calib, (1,2), (3,2))
    """

    if data.ndim == 3:
        data = data[:,:,:,None]

    if calib.ndim == 3:    
        calib = calib[:,:,:,None]

    if len(R) == 2:
        R = R + (1,)

    if len(kernel) == 2:
        kernel = kernel + (1,)

    # initialize output
    recon = data.copy()

    # find sampling mask
    sampling_mask = data != 0

    # generate different kernel geometries
    kernels = build_kernels(R, kernel)

    # loop over possible kernel types
    for ker in kernels:

        # get src, trg points
        src, trg = get_training_data(calib, ker)

        # fit/estimate weights
        ker.weights = trg@np.linalg.pinv(src)

        # apply weights to reconstruct missing data 
        apply_kernel(recon, sampling_mask, ker)

    return recon.squeeze()

class GrappaKernel:
    """
    One kernel type corresponding to a particular missing point
    within the Ry x Rz acceleration cell.
    """

    def __init__(self, source_offsets, target_offset):
        self.source_offsets = np.asarray(source_offsets, dtype=int)
        self.target_offset = np.asarray(target_offset, dtype=int)
        self.weights = None


def build_kernels(R, kernel_size):
    """
    Create all GRAPPA kernel geometries.
    """

    Rx, Ry, Rz = R
    kx, ky, kz = kernel_size

    if Rx != 1:
        raise ValueError("This implementation assumes x is fully sampled")

    cx = kx // 2
    cy = ky // 2
    cz = kz // 2

    source_offsets = []

    for dx in range(-cx, cx + 1):
        for dy in range(-cy, cy + 1):
            for dz in range(-cz, cz + 1):
                source_offsets.append((dx * Rx, dy * Ry, dz * Rz))

    kernels = []

    for ry in range(Ry):
        for rz in range(Rz):

            if ry == 0 and rz == 0:
                continue

            target_offset = np.array([0, ry, rz])
            kernels.append(GrappaKernel(source_offsets, target_offset))

    return kernels


def apply_kernel(data, sampling_mask, kernel):
    """
    Reconstruct one kernel type.
    """

    nc, nx, ny, nz = data.shape

    src_offsets = kernel.source_offsets
    trg_offset = kernel.target_offset

    maxx = max(abs(o[0]) for o in src_offsets)
    maxy = max(abs(o[1]) for o in src_offsets)
    maxz = max(abs(o[2]) for o in src_offsets)

    W = kernel.weights

    for x in range(maxx, nx - maxx):
        for y in range(maxy, ny - maxy):
            for z in range(maxz, nz - maxz):

                tx = x + trg_offset[0]
                ty = y + trg_offset[1]
                tz = z + trg_offset[2]

                if sampling_mask[0, tx, ty, tz]:
                    continue

                src = []

                for dx, dy, dz in src_offsets:
                    src.append(data[:,x + dx, y + dy, z + dz])

                src = np.concatenate(src)

                data[:, tx, ty, tz] = W @ src

def get_training_data(calib, kernel):
    """
    Returns

    src : source points
    trg : target points
    """

    nc, nx, ny, nz = calib.shape

    src_offsets = kernel.source_offsets
    trg_offset = kernel.target_offset

    maxx = max(abs(o[0]) for o in src_offsets)
    maxy = max(abs(o[1]) for o in src_offsets)
    maxz = max(abs(o[2]) for o in src_offsets)

    src = []
    trg = []

    for x in range(maxx, nx - maxx):
        for y in range(maxy, ny - maxy):
            for z in range(maxz, nz - maxz):

                A = []

                for dx, dy, dz in src_offsets:
                    A.append(calib[:, x + dx, y + dy, z + dz,])

                A = np.concatenate(A)

                tx = x + trg_offset[0]
                ty = y + trg_offset[1]
                tz = z + trg_offset[2]

                B = calib[:, tx, ty, tz]

                src.append(A)
                trg.append(B)

    src = np.asarray(src).T
    trg = np.asarray(trg).T

    return src, trg
