import ismrmrd
from ismrmrdtools import coils
import os
import itertools
import logging
import traceback
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
import xml.dom.minidom
import base64
import ctypes
import re
import ismrmrd_server.mrdhelper as mrdhelper
import ismrmrd_server.constants as constants
from time import perf_counter
from tempfile import mkdtemp
from pathlib import Path
import nibabel as nib
import pandas as pd
import gc
from dipy.denoise.localpca import mppca

from common import SiemensRAW, KSPACE_LAYOUT
from grappa import grappa

# Folder for debug output files
debugFolder = "/tmp/share/debug"

def process(connection, config, mrdHeader):
    logging.info("Config: \n%s", config)

    # mrdHeader should be xml formatted MRD header, but may be a string
    # if it failed conversion earlier
    try:
        logging.info("Incoming dataset contains %d encodings", len(mrdHeader.encoding))
        logging.info("First encoding is of type '%s', with a matrix size of (%s x %s x %s) and a field of view of (%s x %s x %s)mm^3", 
            mrdHeader.encoding[0].trajectory, 
            mrdHeader.encoding[0].encodedSpace.matrixSize.x, 
            mrdHeader.encoding[0].encodedSpace.matrixSize.y, 
            mrdHeader.encoding[0].encodedSpace.matrixSize.z, 
            mrdHeader.encoding[0].encodedSpace.fieldOfView_mm.x, 
            mrdHeader.encoding[0].encodedSpace.fieldOfView_mm.y, 
            mrdHeader.encoding[0].encodedSpace.fieldOfView_mm.z)

    except:
        logging.info("Improperly formatted MRD header: \n%s", mrdHeader)

    # Continuously parse incoming data parsed from MRD messages
    try:
        raw = SiemensRAW(mrdHeader)
        for item in connection:
            # ----------------------------------------------------------
            # Raw k-space data messages
            # ----------------------------------------------------------
            if isinstance(item, ismrmrd.Acquisition):
                raw.add_acq(item)
            
                # Process one repetition at the time (for now?)
                if item.is_flag_set(ismrmrd.ACQ_LAST_IN_REPETITION):
                    images = process_raw(raw, mrdHeader)
                    connection.send_image(images)
                    raw.reset_acq()

            elif item is None:
                break

            else:
                logging.error("Unsupported data type %s", type(item).__name__)

    except Exception as e:
        logging.error(traceback.format_exc())
        connection.send_logging(constants.MRD_LOGGING_ERROR, traceback.format_exc())

    finally:
        connection.send_close()

def build_reference_volume(kspace, acs_mask, nKx, nKy, nKx_recon, fov_x, fov_y, output_dir, raw, CROP_SIZE):
    """
    Reconstruct reference volume from echo 0 and save as NIfTI.
    Slices are stored in anatomical order (inf→sup) for SCT centerline detection.

    Parameters
    ----------
    kspace    : (1, 1, 1, 1, nEcho, nSlice, 1, nKy, nKx, nCoils)
    acs_mask  : same shape as kspace
    output_dir : Path — where to save the NIfTI
    raw        : SiemensRAW object — needed to extract physical slice positions

    Returns
    -------
    ref_path : Path — path to saved NIfTI, or None if failed
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    ref_path = output_dir / "ref_echo0.nii.gz"
    ref_path_full = output_dir / "ref_echo0_full.nii.gz"

    if ref_path.exists():
        print(f"Reference volume already present : {ref_path}")
        return ref_path

    print("Building reference volume from echo 0...")

    pixel_size_x = fov_x / nKx_recon
    pixel_size_y = fov_y / nKy

    # --------------------------------------------------
    # GRAPPA RECONSTRUCTION ON ECHO 0
    # --------------------------------------------------
    kspace_ec0   = kspace[:, :, :, :, [0], :, :, :, :, :]

    # squeeze (rep, phase, set, segment, echo, kz) → (nSlice, nKy, nKx, nCoils)
    squeeze_axes = (0, 1, 2, 3, 4, 6)
    kspace_ec0    = np.squeeze(kspace_ec0,   axis=squeeze_axes)

    print(f"  kspace_sq.shape : {kspace_ec0.shape}")

    kspace_grappa = grappa_reconstruction(kspace_ec0, acs_mask)
    print(f"  kspace_grappa.shape : {kspace_grappa.shape}")
    del kspace_ec0

    # Restore to (rep=1, echo=1, slice, y, x, coils) for raw_to_image
    kspace_grappa = kspace_grappa[np.newaxis, np.newaxis, :, :, :, :]
    print(f"  kspace_restored.shape : {kspace_grappa.shape}")

    data = np.flip(kspace_grappa, 4) # x axis inversion
    data = reconstruct_image(data)
    del kspace_grappa

    # Coil combination (RMS, coils_axis=-1)
    data = np.sqrt(np.sum(np.abs(data)**2, axis=-1))

    # Remove readout oversampling by cropping
    images = remove_oversampling(data, nKx, nKx_recon, 4)

    images = mag_images(images)
    print(f"  images.shape : {images.shape}")

    print("nKy=", nKy)
    print("nKx_recon=", nKx_recon)

    # images : (1, 1, nSlice, nKy_recon, nKx_recon)
    images_slices = images[0, 0]   # (nSlice, nKy_recon, nKx_recon)
    print(f"  images_slices.shape : {images_slices.shape}")

    # Crop parameters
    CROP_HALF = CROP_SIZE // 2
    x_lo  = (nKx_recon // 2) - CROP_HALF
    y_lo  = (nKy // 2) - CROP_HALF

    # --------------------------------------------------
    # GET PHYSICAL SLICE POSITIONS → ANATOMICAL ORDER
    # --------------------------------------------------
    slice_z_positions = {}
    for acq in raw.acquisitions:
        sl = acq.idx.slice
        if sl not in slice_z_positions:
            slice_z_positions[sl] = acq.position[2]

    sorted_sl_indices = sorted(slice_z_positions, key=slice_z_positions.get)

    print(f"\n  Anatomical slice order :")
    for anat_idx, sl in enumerate(sorted_sl_indices):
        print(f"    anat[{anat_idx}] = ISMRMRD sl={sl}  "
              f"z={slice_z_positions[sl]:.2f} mm")

    # --------------------------------------------------
    # FILL ref_volume IN ANATOMICAL ORDER
    # --------------------------------------------------
    n_sorted   = len(sorted_sl_indices)
    ref_volume = np.zeros((n_sorted, CROP_SIZE, CROP_SIZE), dtype=np.float32)
    ref_volume_full = np.zeros((n_sorted, nKx_recon, nKy), dtype=np.float32)

    # dz from physical slice positions — more reliable than FOV z / nSlice
    # which gives slice thickness (5mm) not slice spacing
    if len(sorted_sl_indices) > 1:
        z0 = slice_z_positions[sorted_sl_indices[0]]
        z1 = slice_z_positions[sorted_sl_indices[1]]
        dz = abs(z1 - z0)
    else:
        dz = 5.0   # fallback

    for anat_idx, sl in enumerate(sorted_sl_indices):
        if sl >= images_slices.shape[0]:
            print(f"    WARNING : sl={sl} out of bounds")
            continue

        img_sl = images_slices[sl]               # (nKy_recon, nKx_recon)
        img_sl = img_sl / (img_sl.max() + 1e-8)  # normalize to [0,1]
        ref_volume[anat_idx] = img_sl[
            y_lo : y_lo + CROP_SIZE,
            x_lo : x_lo + CROP_SIZE
        ]
        ref_volume_full[anat_idx] = img_sl

    # --------------------------------------------------
    # FLIP y + SIMPLE DIAGONAL AFFINE (MRINavigator.jl convention)
    # y flip : SCT expects y=0 at top (radiological convention)
    # affine : voxel sizes only, no patient orientation
    # --------------------------------------------------
    affine = np.diag([-pixel_size_x, pixel_size_y, dz, 1.0]).astype(np.float32)

    nii = nib.Nifti1Image(ref_volume.T, affine)
    nib.save(nii, str(ref_path))

    nii_full = nib.Nifti1Image(ref_volume_full.T, affine)
    nib.save(nii_full, str(ref_path_full))

    print(f"\nReference volume saved : {ref_path}")
    print(f"  shape={ref_volume.shape}  "
          f"range=[{ref_volume.min():.4f}, {ref_volume.max():.4f}]")

    return ref_path


def run_sct_deepseg(output_dir, nKy, nKx_recon, nSlice, CROP_SIZE):
    """
    Segment the spinal cord with SCT deepseg on the cropped reference volume
    (ref_echo0.nii.gz) and extract the center of mass of the
    segmentation as the spinal cord center.

    Parameters
    ----------
    output_dir : Path
        Directory containing ref_echo0.nii.gz.
    nKy : int
        Size of the full image in the y direction.
    nKx_recon : int
        Size of the full image in the x direction.
    nSlice : int
        Number of slices.
    CROP_SIZE : int
        Size of the cropped reference image.

    Returns
    -------
    csv_path : Path or None
        Path to the centerline CSV, or None if SCT failed.
    """
    import subprocess

    output_dir  = Path(output_dir)
    ref_path    = output_dir / "ref_echo0.nii.gz"
    seg_path = output_dir / "ref_echo0_seg.nii.gz"
    csv_path    = output_dir / "ref_echo0_centerline.csv"
    csv_path_crop = output_dir / "ref_echo0_centerline_crop.csv"

    if csv_path.exists():
        print(f"Centerline CSV already present : {csv_path}")
        return csv_path

    if not ref_path.exists():
        print(f"ERROR : reference volume not found : {ref_path}")
        return None

    print(f"Running sct_deepseg spinalcord ...")
    print(f"  Input  : {ref_path}")
    print(f"  Output : {seg_path}")

    result = subprocess.run(
        [   "sct_deepseg",
            "spinalcord",
            "-i", str(ref_path.resolve()),
            "-o", str(seg_path.resolve()),
        ],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print("SCT deepseg failed -> check SCT installation and PATH")
        return None
    
    # Find SCT output segmentation
    if not seg_path.exists():
        print(f"SCT segmentation not found : {seg_path}")
        return None

    # Load segmentation
    seg_img = nib.load(seg_path)
    seg_mask = seg_img.get_fdata()

    # Check that segmentation and reference have same shape
    if seg_mask.shape != nib.load(ref_path).shape:
        print("ERROR : segmentation and reference volume shapes differ.")
        return None

    # Extract centerline from segmentation
    print("\n=== Extracting spinal cord centerline ===")
    centerline = []

    for z in range(nSlice):

        # Binary spinal cord mask for current slice
        mask = seg_mask[:, :, z] > 0

        if not np.any(mask):
            print(f"  WARNING : no spinal cord detected at slice z={z}")
            continue

        x_coords, y_coords = np.where(mask)

        # Center of mass of the segmented spinal cord
        x_center = np.mean(x_coords)
        y_center = np.mean(y_coords)

        centerline.append((x_center, y_center, z))

        print(f" z={z:2d} → center_crop=({x_center:.1f}, {y_center:.1f})")

    if len(centerline) == 0:
        print("ERROR : no spinal cord detected in segmentation.")
        return None
    
    # Crop coordinates
    CROP_HALF = CROP_SIZE // 2
    x_lo_crop = (nKx_recon // 2) - CROP_HALF
    y_lo_crop = (nKy // 2) - CROP_HALF

    print("\n=== Crop → full image coordinate transformation ===")
    print(f"  CROP_SIZE = {CROP_SIZE}")
    print(f"  x_lo_crop = {x_lo_crop}")
    print(f"  y_lo_crop = {y_lo_crop}")

    # --------------------------------------------------
    # CSV 1 — cropped image space
    # --------------------------------------------------
    with open(csv_path_crop, "w") as f:
        f.write("x_crop,y_crop,z_anat\n")

        for x_crop, y_crop, z in centerline:
            f.write(
                f"{x_crop:.4f},{y_crop:.4f},{z}\n"
            )

    print(f"\nCenterline in cropped image space saved :{csv_path_crop}")

    # --------------------------------------------------
    # CSV 2 — full image space
    # --------------------------------------------------
    with open(csv_path, "w") as f:

        for x_crop, y_crop, z in centerline:

            # No flip is applied here because ref_echo0.nii.gz
            # was not flipped before being given to SCT.
            x_img = x_crop + x_lo_crop
            y_img = y_crop + y_lo_crop

            f.write(
                f"{x_img:.4f},{y_img:.4f},{z}\n"
            )

            print(
                f"  z={z:2d} → "
                f"crop({x_crop:.1f},{y_crop:.1f}) → "
                f"img({x_img:.1f},{y_img:.1f})"
            )

    print(f"\nCenterline in full image space saved : {csv_path}")

    return csv_path

def apply_nav_mask_from_centerline(S, center_x_per_slice, nSlice, nKx, nKx_recon, width=35):
    """
    Apply spatial masks to all navigator lines around the spinal cord centerline.
    Navigators are FFT'd to spatial domain, masked per slice, and returned masked.

    S (rep, samples, ky, sl, coil): complex — navigator kspace lines
    center_x_per_slice {sl: x_img}:  centerline x in image coords (0 to nKx_recon)
    nKx_recon (int): image size along x (after oversampling removal)
    width (int): half-width of mask in navigator spatial domain
    """

    # FFT to spatial domain: one batched FFT over the samples axis
    nav_spatial = np.fft.fftshift(np.fft.fft(S, axis=1), axes=1)
   
    # Build (samples, sl) mask
    mask = np.zeros((nKx, nSlice))
    for sl, x_img in center_x_per_slice.items():
        # Scale x_img (image coords) to navigator spatial domain
        x_center = int(x_img / nKx_recon * nKx)
        # Rectangular mask
        lo = max(0, x_center - width)
        hi = min(nKx, x_center + width)
        mask[lo:hi, sl] = 1

    # Broadcast (samples, sl) → (rep, samples, ky, sl, coil)
    nav_masked = nav_spatial * mask[None, :, None, :, None]

    return nav_masked

def process_raw(raw, mrdHeader, use_memmap=False):

    # Metadata from ISMRMRD header
    acq_metadata=SiemensRAW(mrdHeader)
    nEcho = acq_metadata.n_echo
    nSlice = acq_metadata.n_slice
    nKx = acq_metadata.n_kx
    nKy =  acq_metadata.n_ky
    nKx_recon = acq_metadata.n_kx_recon
    echo_times = np.array(acq_metadata.echo_times, dtype=np.float32) * 1e-3 
    FOV_x = acq_metadata.FOV_x
    FOV_y = acq_metadata.FOV_y
    FOV_z = acq_metadata.FOV_z

    # Metadata from acquisitions
    rep_index = raw.acquisitions[0].idx.repetition
    dt = raw.acquisitions[0].sample_time_us * 1e-6

    # Navigator echo time is not inside the header. It should be under user_int 
    # from acquisition metadata but is not present when using FIRE
    navigator_te = 24e-3   

    # Preprocessing acquisitions
    # First repetition will contain a noise acq. Extract it and keep it for all reps.
    if raw.noise is None:
        raw.extract_noise()
    
    # Remove useless acquisitions at the beginning of each rep
    raw.remove_phase_stabilization_references()

    # Build kspace and navigator data structure
    print("Building kspace array...")
    kspace, navigator, acs_mask = raw.build_kspace(use_memmap=use_memmap)
    kspace = kspace[[rep_index], ...]           # (1, 1, 1, 1, nEcho=4, nSlice=15, 1, nKy=384, nKx=768, nCoils=4)
    navigator = navigator[[rep_index], ...]     # (1, 1, 1, 1, 1,       nSlice=15, 1, nKy=384, nKx=768, nCoils=4)

    # Save it for tests
    #raw.save_kspace("ice_data.npz")
    
    # Load precomputed kspace
    #raw.load_kspace("ice_data.npz")

    # --------------------------------------------------
    # Build reference volume and save it under NifTi
    # TODO: Change this path to something else than workspaces
    CENTERLINE_DIR = Path("/workspaces/siemens-fire/Icesimu_output/sct_centerline")
    CROP_SIZE = 150

    ref_path = build_reference_volume(kspace, acs_mask, nKx, nKy, nKx_recon, FOV_x, FOV_y, CENTERLINE_DIR, raw, CROP_SIZE)
    print(f"Reference volume ready : {ref_path}")
    print("SCT centerline detection will now be run on this volume.")

    # Run SCT centerline detection — results saved to CENTERLINE_DIR for inspection
    csv_path = run_sct_deepseg(CENTERLINE_DIR, nKy, nKx_recon, nSlice, CROP_SIZE)

    if csv_path is not None:
        print(f"Centerline detection successful : {csv_path}")
    else:
        print(f"Centerline detection failed — continuing without masking")

    # --------------------------------------------------
    # LOAD CENTERLINE CSV
    # --------------------------------------------------
    center_x_per_slice = {}

    if csv_path is not None:
        df_cl = pd.read_csv(csv_path, header=None, names=["x", "y", "z"])
        center_x_per_slice = {int(row["z"]): float(row["x"]) for _, row in df_cl.iterrows()}
        print(f"Centerline loaded : {len(center_x_per_slice)} slices")
        use_mask = True
    else:
        print("No centerline — using full navigator line")
        use_mask = False

    # --------------------------------------------------
    # NAVIGATOR PREPARATION
    # S : (rep, samples=nKx, lines=nKy, slices=nSlice, coils)
    # --------------------------------------------------
    # Axis order is derived from KSPACE_LAYOUT (source used
    # by _get_kspace_dims), so a reordering there propagates here automatically.

    axis_names = list(KSPACE_LAYOUT) + ["kx", "coil"]

    # axes kept in S: rep, slice, ky, kx, coil: all others must be singleton
    KEEP = {"repetition", "slice", "kspace_encoding_step_1", "kx", "coil"}

    # drop singleton axes
    squeeze_axes = tuple(i for i, n in enumerate(axis_names) if n not in KEEP)
    S = navigator.squeeze(axis=squeeze_axes)

    # reorder to (rep, kx, ky, slice, coil)
    remaining = [n for n in axis_names if n in KEEP]          # order after squeeze
    target    = ["repetition", "kx", "kspace_encoding_step_1", "slice", "coil"]
    S = np.transpose(S, [remaining.index(n) for n in target])

    # --------------------------------------------------
    # CENTERLINE MASKING on navigator lines
    # S[rep, samples, ky, sl, coil] — samples axis = nav readout
    # apply_nav_mask_from_centerline works on (nKx,) 1D line
    # center_x_per_slice is in image coordinates (0 to nKx_recon)
    # --------------------------------------------------
    if use_mask:
        print("Applying centerline mask to navigator...")
        S = apply_nav_mask_from_centerline(S, center_x_per_slice, nSlice, nKx, nKx_recon, width=35)
        print("Centerline masking applied.")

    # --------------------------------------------------
    # PHASE EXTRACTION — unchanged from original
    # --------------------------------------------------
    print("Computing corrections...")
    phase_extractions = np.stack([phase_extraction(s, raw.noise.data.T) for s in S])
    field_estimates = field_conversion(phase_extractions, navigator_te) # rad/s

    # Apply navigator correction
    print("Applying corrections...")
    print("Start processing slices sequentially to reduce memory usage")
    # Store corrected images for all slices

    if use_memmap:
        filename = os.path.join(mkdtemp(), "corrected.dat")
        corrected = np.memmap(filename, dtype=np.complex64, mode='w+', shape=kspace.shape)

    print("K-Space correction...")
    corrected[:] = kspace_correction(kspace, field_estimates, nKx, echo_times, dt) # (1, 4, 1, 384, 768, 24)

    # Use GRAPPA to fill in missing kspace lines
    print("GRAPPA correction...")
    corrected[:] = grappa_reconstruction(corrected, acs_mask)
    print("Reconstruction...")
    img = raw_to_image(corrected, nKx, nKx_recon)
    img = mag_images(img)

    # Save images for faster testing
    # np.save("images.npy", images)

    print("Denoising...")
    # Reshape for MPPCA : (nEcho, nRep, nSlice, nKy, nKx)
    imgs_for_denoise = corrected[0]   # (4, 15, 384, 384) = (echo, slice, y, x)
    imgs_for_denoise = imgs_for_denoise[:, np.newaxis, :, :, :] # (4, 1, 15, 384, 384) = (echo, rep, slice, y, x)

    imgs_denoised = denoise_mppca(imgs_for_denoise, patch_radius=2)

    # Reshape back
    corrected = imgs_denoised[:, 0, :, :, :][np.newaxis, :, :, :, :] # (1, 4, 15, 384, 384) = (rep, echo, slice, y, x)

    field_of_view = (ctypes.c_float(FOV_x), ctypes.c_float(FOV_y), ctypes.c_float(FOV_z))

    header_map = {}
    for acq in raw.acquisitions:
        key = (acq.idx.contrast, acq.idx.slice)
        if key not in header_map:
            header_map[key] = acq.getHead()
    acq_headers = [header_map[(c, s)] for c in range(nEcho) for s in range(nSlice)]
    for i, h in enumerate(acq_headers):
        print(f"Image {i}: contrast={h.idx.contrast}, slice={h.idx.slice}")
    ismrmrd_images = convert_to_ismrmrd_images(corrected, acq_headers, field_of_view)
    
    print("Sending images...")
    return ismrmrd_images

def mag_images(images):
    return np.abs(images)

def convert_to_int16(data):
    return np.around(data * (2**12 - 1)/data.max()).astype(np.int16)

def convert_to_ismrmrd_images(images, acq_headers, fov):
    images_out = []
    *_, y, x = images.shape
    data = convert_to_int16(images)
    for i, img in enumerate(data.reshape(-1, y, x)):
        ismrmrd_image = ismrmrd.Image.from_array(img, transpose=False)
        ismrmrd_image.setHead(mrdhelper.update_img_header_from_raw(ismrmrd_image.getHead(), acq_headers[i]))
        ismrmrd_image.field_of_view = fov
        ismrmrd_image.image_index = i
        tmp_meta = ismrmrd.Meta()
        tmp_meta['DataRole']               = 'Image'
        tmp_meta['ImageProcessingHistory'] = ['FIRE', 'PYTHON']
        tmp_meta['Keep_image_geometry']    = 1
        ismrmrd_image.attribute_string = tmp_meta.serialize()
        images_out.append(ismrmrd_image)
    return images_out

def raw_to_image(raw, nKx, nKx_recon, noise_data=None):
    # assumed shape : (repetitions, echoes, slices, y, x, coils)
    #                 (0          , 1     , 2,    , 3, 4, 5) 

    data = np.flip(raw, (3, 4)) # inverser les données en x et y for some reason
    data = reconstruct_image(data)
    data *= np.prod(data.shape) # FFT scaling, for consistency with ICE apparently

    # Walsh coil combination
    print("  Using Inati coil combination with prewhitening...")
    data = coil_combination_Inati(data, noise_data=noise_data)

    # Remove readout oversampling by cropping
    data = remove_oversampling(data, nKx, nKx_recon, 4)

    return data

def remove_oversampling(data, nKx, nKx_recon, readout_axis):
    start = (nKx - nKx_recon) // 2
    
    return data.take(np.arange(start, start+nKx_recon), axis=readout_axis)

def coil_combination_Inati(data, noise_data=None, smoothing=5, niter=3):
    """
    Coil combination using Inati iterative method with prewhitening.
    
    Pipeline :
        1. Prewhitening: decorrelates coil noise
        2. Inati CSM estimation: estimates sensitivity maps from image itself
        3. Sensitivity-weighted combination: optimal SNR combination
    
    Parameters
    ----------
    data (rep, echo, slice, y, x, coils): complex coil images
    noise_data (nCoils, nSamples): raw.noise.data
    smoothing (int): smoothing kernel for Walsh CSM (default 5)
    niter (int): Walsh power iterations (default 3)

    Returns
    -------
    combined (rep, echo, slice, y, x): magnitude combined image
    """
    data = data[0, :, 0, ...]
    echo,  y, x, _ = data.shape

    # Prewhitening matrix
    if noise_data is not None:
        dmtx = coils.calculate_prewhitening(noise_data)
        print(f"  Prewhitening matrix computed: shape={np.asarray(dmtx).shape}")
    else:
        dmtx = None

    combined = np.zeros((echo, y, x), dtype=np.float32)

    for e in range(echo):
        # Extract coil images : (y, x, nCoils) → (nCoils, y, x)
        img_coils = np.moveaxis(data[e], -1, 0)

        # Prewhitening
        if dmtx is not None:
            img_coils = coils.apply_prewhitening(img_coils, dmtx)

        # Inati CSM estimation + combination
        _, combined_complex = coils.calculate_csm_inati_iter(
            img_coils,
            smoothing=smoothing,
            niter=niter,
            thresh=1e-3
        )

        combined[e] = np.abs(combined_complex)

    return combined[None, :, None, ...]

def reconstruct_image(kspace, axes=(3, 4)):
    # First ifftshift, because numpy assumes the DC component to be at index 0.
    # Physically, the acquisition has the DC component at its center and the high frequencies at its edges

    # Preallocate an array on disk for our results
    image = np.memmap(os.path.join(mkdtemp(), 'reconstruction.dat'), dtype=np.complex64, mode='w+', shape=kspace.shape)

    image[:] = np.fft.ifftshift(kspace, axes=axes)
    # Inverse FFT to get the image
    image[:] = np.fft.ifft2(image, axes=axes)
    # Pour replacer l'objet au centre de l'image?
    image[:] = np.fft.fftshift(image, axes=axes)

    return image

def phase_extraction(navigator, noise):
    """
    Arguments:
    navigator -- shape : (j, l, p, c) -> (samples, lines, slices, coils)
    noise     -- shape : (j, c)       -> (samples, coils)
    """

    # subtract first navigator phase to remove static phase contributions
    delta_S = (navigator * np.exp(-1j*np.angle(navigator[:, [0], :, :]))).astype(np.complex64)  # (samples=768, lines=384, slices=15, coils=(4,8,...))

    w = np.abs(delta_S) / np.std(noise, axis=0)   # (768, 384, 15, coils=(4,8,...))  (TODO: check if should need to raise to power 2)
    # RuntimeWarning here because of dividing by zero
    w_tilde = w/np.sum(w, axis=(0, 3), keepdims=True) # (768, 384, 15, coils=(4,8,...))
    # Replace resulting NaNs by zero
    w_tilde[np.isnan(w_tilde)] = 0.0

    delta_S = np.sum(w_tilde * delta_S, axis=(0, 3)) # (384, 15)  (lines, slices)

    delta_phi_mean = np.angle(np.mean(delta_S, axis=0)) # (15,)  (slices,)

    delta_S_tilde = delta_S * np.exp(-1j * delta_phi_mean) # (384, 15)  (lines, slices)

    delta_phi = np.angle(delta_S_tilde) # (384, 15)  (lines, slices)

    return delta_phi

def field_conversion(nav_phases, te_nav):
    return nav_phases/te_nav    # rad/s

def kspace_correction(raw_data, field_estimates, n_samples, echo_times, dt):
    t = np.array([(j - n_samples/2)*dt for j in range(n_samples)], dtype=np.float32) # (768,) (samples,)
    t = np.repeat(t[:, np.newaxis], echo_times.shape[0], axis=1) # (768, 4) (samples, echo)
    t += echo_times # will probably crash  # (768, 4) (samples, echo)

    demodulation = np.exp(-1j * np.einsum('rlp,je->replj', field_estimates, t)).astype(np.complex64)
    print("demodulation.shape", demodulation.shape) # (1, 4, 15, 384, 768) (1, echo, slices, lines, samples)

    # TODO: handle all the dimensions correctly instead of squeezing
    corrected = raw_data*demodulation[..., np.newaxis] # (1, 4, 1, 384, 768, coils=(4,8,...)) (1, echo, slices, lines, samples, coils)
    return corrected
 
def grappa_reconstruction(kspace, acs_lines, R=2, kernel_size=(5, 5)):
    
    *leading, y, x, c = kspace.shape

    kspace_r   = kspace.reshape(-1, y, x, c)

    results = np.memmap(
        os.path.join(mkdtemp(), 'grappa.dat'),
        dtype=np.complex64, mode='w+',
        shape=kspace_r.shape
    )

    for i in range(kspace_r.shape[0]):
        # kspace_r[i]   : (nKy, nKx, nCoils) → reorder to (nCoils, nKy, nKx)
        k_coils   = np.moveaxis(kspace_r[i],   -1, 0)   # (nCoils, nKy, nKx)

        # Extract ACS line indices from mask
        # acs_lines[0] is boolean mask of shape (nKy,)
        # ACS lines are those where the entire ky line is non-zero
        acs_idx = np.where(acs_lines)[0]

        if len(acs_idx) == 0:
            print(f"  WARNING : no ACS lines found for batch {i}")
            results[i] = kspace_r[i]
            continue

        # Apply GRAPPA — returns (nCoils, nKy, nKx)
        k_filled = apply_grappa(k_coils, acs_idx, R=R, kernel_size=kernel_size)

        # Reorder back to (nKy, nKx, nCoils)
        results[i] = np.moveaxis(k_filled, 0, -1)
    del kspace_r
    results.flush()
    return results.reshape(*leading, y, x, c)

def apply_grappa(kspace_2d, acs_lines, R, kernel_size=(5, 5)):
    
    if R is None:
        raise ValueError("R cannot be None")

    k = np.swapaxes(kspace_2d, 1, 2).astype(np.complex64)  # (nCoil, nKx, nKy)
    calib = k[:, :, acs_lines].copy()                       # (nCoil, nKx, nACS)

    k_filled = grappa(
        data   = k,
        calib  = calib,
        R      = (1, int(R)),    # int() pour garantir le type
        kernel = kernel_size
    )

    return np.swapaxes(k_filled, 1, 2)   # (nCoil, nKy, nKx)

def denoise_mppca(images, patch_radius=2):
    """
    MPPCA denoising on multi-echo multi-rep magnitude images.
    Exploits redundancy across echoes AND repetitions.

    images (nEcho, nRep, nSlice, nKy, nKx): magnitude
    patch_radius (int): spatial patch size = (2r+1)³

    Returns
    -------
    denoised : same shape as images
    """    

    nEcho, nRep, nSlice, nKy, nKx = images.shape

    # MPPCA works on 4D volume (x, y, z, volumes)
    # Reshape to (nKy, nKx, nSlice, nEcho*nRep)
    imgs_4d = images.reshape(nEcho*nRep, nSlice, nKy, nKx)
    imgs_4d = np.moveaxis(imgs_4d, 0, -1)   # (nSlice, nKy, nKx, nVol)
    imgs_4d = np.moveaxis(imgs_4d, 0, 2)    # (nKy, nKx, nSlice, nVol)
    imgs_4d = imgs_4d.astype(np.float64)

    print(f"  MPPCA input shape : {imgs_4d.shape}")

    denoised_4d, sigma = mppca(imgs_4d, patch_radius=patch_radius,
                                return_sigma=True)

    print(f"  Estimated noise sigma : {sigma.mean():.4f}")

    # Reshape back
    denoised_4d = np.moveaxis(denoised_4d, 2, 0)   # (nSlice, nKy, nKx, nVol)
    denoised_4d = np.moveaxis(denoised_4d, -1, 0)  # (nVol, nSlice, nKy, nKx)
    denoised    = denoised_4d.reshape(nEcho, nRep, nSlice, nKy, nKx)

    return denoised