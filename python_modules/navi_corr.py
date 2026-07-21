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
import mrdhelper
import constants
from time import perf_counter
from tempfile import mkdtemp
from pathlib import Path
import nibabel as nib
import pandas as pd

from common import SiemensRAW
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

def build_reference_volume(kspace, acs_mask, mrdHeader, output_dir, raw):
    """
    Reconstruct reference volume from echo 0 and save as NIfTI.
    Slices are stored in anatomical order (inf→sup) for SCT centerline detection.

    Parameters
    ----------
    kspace    : (1, 1, 1, 1, nEcho, nSlice, 1, nKy, nKx, nCoils)
    acs_mask  : same shape as kspace
    mrdHeader : ISMRMRD header
    output_dir : Path — where to save the NIfTI
    raw        : SiemensRAW object — needed to extract physical slice positions

    Returns
    -------
    ref_path : Path — path to saved NIfTI, or None if failed
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    ref_path = output_dir / "ref_echo0.nii.gz"

    if ref_path.exists():
        print(f"Reference volume already present : {ref_path}")
        return ref_path

    print("Building reference volume from echo 0...")

    nSlice    = kspace.shape[5]
    nKy       = kspace.shape[7]
    enc       = mrdHeader.encoding[0]
    fov_x     = enc.reconSpace.fieldOfView_mm.x
    fov_y     = enc.reconSpace.fieldOfView_mm.y
    nKx_recon = enc.reconSpace.matrixSize.x
    pixel_size_x = fov_x / nKx_recon
    pixel_size_y = fov_y / nKy
    # dz           = enc.reconSpace.fieldOfView_mm.z / nSlice \
    #                if enc.reconSpace.fieldOfView_mm.z > 0 else 5.0

    # --------------------------------------------------
    # GRAPPA RECONSTRUCTION ON ECHO 0
    # --------------------------------------------------
    kspace_ec0   = kspace[:, :, :, :, [0], :, :, :, :, :]
    acs_mask_ec0 = acs_mask[:, :, :, :, [0], :, :, :, :, :]

    # squeeze (rep, phase, set, segment, echo, kz) → (nSlice, nKy, nKx, nCoils)
    squeeze_axes = (0, 1, 2, 3, 4, 6)
    kspace_sq    = np.squeeze(kspace_ec0,   axis=squeeze_axes)
    acs_mask_sq  = np.squeeze(acs_mask_ec0, axis=squeeze_axes)

    print(f"  kspace_sq.shape : {kspace_sq.shape}")

    kspace_grappa = grappa_reconstruction(kspace_sq, kspace_sq * acs_mask_sq)
    print(f"  kspace_grappa.shape : {kspace_grappa.shape}")

    # Restore to (rep=1, echo=1, slice, y, x, coils) for raw_to_image
    kspace_restored = kspace_grappa[np.newaxis, np.newaxis, :, :, :, :]
    print(f"  kspace_restored.shape : {kspace_restored.shape}")

    images = raw_to_image(kspace_restored)
    images = mag_images(images)
    print(f"  images.shape : {images.shape}")

    # images : (1, 1, nSlice, nKy_recon, nKx_recon)
    images_slices = images[0, 0]   # (nSlice, nKy_recon, nKx_recon)
    print(f"  images_slices.shape : {images_slices.shape}")

    nKy_recon  = images_slices.shape[1]
    nKx_recon_ = images_slices.shape[2]

    # Crop parameters
    CROP_SIZE = 150
    CROP_HALF = CROP_SIZE // 2
    cx    = nKx_recon_ // 2
    cy    = nKy_recon  // 2
    x_lo  = cx - CROP_HALF
    y_lo  = cy - CROP_HALF

    # --------------------------------------------------
    # GET PHYSICAL SLICE POSITIONS → ANATOMICAL ORDER
    # --------------------------------------------------
    slice_z_positions = {}
    for acq in raw.acquisitions:
        sl = acq.idx.slice
        if sl not in slice_z_positions:
            slice_z_positions[sl] = acq.position[2]

    sorted_sl_indices = sorted(
        slice_z_positions.keys(),
        key=lambda sl: slice_z_positions[sl]
    )

    print(f"\n  Anatomical slice order :")
    for anat_idx, sl in enumerate(sorted_sl_indices):
        print(f"    anat[{anat_idx}] = ISMRMRD sl={sl}  "
              f"z={slice_z_positions[sl]:.2f} mm")

    # --------------------------------------------------
    # FILL ref_volume IN ANATOMICAL ORDER
    # --------------------------------------------------
    n_sorted   = len(sorted_sl_indices)
    ref_volume = np.zeros((n_sorted, CROP_SIZE, CROP_SIZE), dtype=np.float32)

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
        
        # dz from physical slice positions — more reliable than FOV z / nSlice
        # which gives slice thickness (5mm) not slice spacing
        if len(sorted_sl_indices) > 1:
            z0 = slice_z_positions[sorted_sl_indices[0]]
            z1 = slice_z_positions[sorted_sl_indices[1]]
            dz = abs(z1 - z0)
        else:
            dz = 5.0   # fallback

    # --------------------------------------------------
    # FLIP y + SIMPLE DIAGONAL AFFINE (MRINavigator.jl convention)
    # y flip : SCT expects y=0 at top (radiological convention)
    # affine : voxel sizes only, no patient orientation
    # --------------------------------------------------
    ref_volume_flipped = ref_volume[:, ::-1, :]
    affine = np.diag([-pixel_size_x, pixel_size_y, dz, 1.0]).astype(np.float32)

    nii = nib.Nifti1Image(ref_volume_flipped.T, affine)
    nib.save(nii, str(ref_path))

    print(f"\nReference volume saved : {ref_path}")
    print(f"  shape={ref_volume.shape}  "
          f"range=[{ref_volume.min():.4f}, {ref_volume.max():.4f}]")

    return ref_path


def run_sct_centerline(output_dir, nKy, nKx_recon, nSlice):
    """
    Run sct_get_centerline on the reference volume and save the centerline CSV.
    Called after build_reference_volume() to detect the spinal cord centerline.

    Parameters
    ----------
    output_dir : Path — directory containing ref_echo0.nii.gz

    Returns
    -------
    csv_path : Path — path to centerline CSV, or None if SCT failed
    """
    import subprocess

    output_dir  = Path(output_dir)
    ref_path    = output_dir / "ref_echo0.nii.gz"
    cl_nii_path = output_dir / "ref_echo0_centerline.nii.gz"
    csv_path    = output_dir / "ref_echo0_centerline.csv"
    csv_path_crop = output_dir / "ref_echo0_centerline_crop.csv"

    if csv_path.exists():
        print(f"Centerline CSV already present : {csv_path}")
        return csv_path

    if not ref_path.exists():
        print(f"ERROR : reference volume not found : {ref_path}")
        return None

    print(f"Running sct_get_centerline -c t2s ...")
    print(f"  Input  : {ref_path}")
    print(f"  Output : {cl_nii_path}")

    result = subprocess.run(
        [
            "sct_get_centerline",
            "-i", str(ref_path.resolve()),
            "-c", "t2s",
            "-o", str(cl_nii_path.resolve()),
        ],
        capture_output=True, text=True
    )

    if result.returncode != 0:
        print("SCT failed — check SCT installation and PATH")
        return None
    
    # Find SCT output CSV
    sct_csv = output_dir / "ref_echo0_centerline.csv"

    if not sct_csv.exists():
        # SCT sometimes appends _centerline to the output name
        sct_csv = output_dir / "ref_echo0_centerline_centerline.csv"

    if not sct_csv.exists():
        print(f"SCT CSV not found. Files in output_dir :")
        for f in output_dir.glob("*"):
            print(f"  {f.name}")
        return None

    # # Print raw SCT output for verification
    df = pd.read_csv(sct_csv, header=None, names=["x", "y", "z"])
    print(f"\nRaw SCT output ({len(df)} points) :")
    print(df.to_string())
    
    # --------------------------------------------------
    # COORDINATE TRANSFORMATION → original image space
    # SCT output is in cropped (150×150) y-flipped volume
    # 1. x_img = SCT_x + x_lo_crop
    # 2. y_img = nKy - (SCT_y + y_lo_crop)
    # --------------------------------------------------
    CROP_SIZE = 150
    CROP_HALF = CROP_SIZE // 2
    cx        = nKx_recon // 2
    cy        = nKy // 2
    x_lo_crop = cx - CROP_HALF
    y_lo_crop = cy - CROP_HALF

    # --------------------------------------------------
    # CSV 1 — cropped volume space (raw SCT output)
    # Coordinates in the 150×150 y-flipped cropped volume
    # Useful for debugging and visualizing on ref_echo0.nii.gz
    # --------------------------------------------------
    with open(csv_path_crop, "w") as f:
        f.write("x_crop,y_crop,z_anat\n")   # header for readability
        for _, row in df.head(nSlice).iterrows():
            z_sl = int(round(row["z"]))
            f.write(f"{row['x']:.4f},{row['y']:.4f},{z_sl}\n")

    print(f"\nCenterline in cropped image space CSV saved : {csv_path_crop}")
    print(f"  (coordinates in 150×150 y-flipped cropped volume)")

    # --------------------------------------------------
    # CSV 2 — original image space (transformed)
    # Coordinates in the full 384×384 image
    # Used by apply_nav_mask_from_centerline in the pipeline
    # --------------------------------------------------
    print(f"\n=== Coordinate transformation → original image space ===")
    print(f"  x_lo_crop={x_lo_crop}, y_lo_crop={y_lo_crop}, nKy={nKy}")

    with open(csv_path, "w") as f:
        for _, row in df.head(nSlice).iterrows():
            z_sl  = int(round(row["z"]))
            x_img = row["x"] + x_lo_crop
            y_img = nKy - (row["y"] + y_lo_crop)
            f.write(f"{x_img:.4f},{y_img:.4f},{z_sl}\n")
            print(f"  z={z_sl:2d} → crop({row['x']:.0f},{row['y']:.0f}) "
                  f"→ img({x_img:.1f},{y_img:.1f})")

    print(f"\nCenterline in full image space CSV saved : {csv_path}")

    return csv_path

def apply_nav_mask_from_centerline(nav_line, x_img, nKx_recon, width=20):
    """
    Apply spatial mask to navigator line around spinal cord centerline.
    Navigator is FFT'd to spatial domain, masked, and returned masked.

    nav_line  : (nKx,) complex — navigator kspace line
    x_img     : float — centerline x in image coordinates (0 to nKx_recon)
    nKx_recon : int — image size along x (after oversampling removal)
    width     : int — half-width of mask in navigator spatial domain
    """
    # FFT to spatial domain
    nav_spatial = np.fft.fftshift(np.fft.fft(nav_line))
    N = len(nav_spatial)

    # Scale x_img (image coords) to navigator spatial domain
    x_center = int(x_img / nKx_recon * N)

    # Apply rectangular mask
    mask = np.zeros_like(nav_spatial)
    lo = max(0, x_center - width)
    hi = min(N, x_center + width)
    mask[lo:hi] = 1

    return nav_spatial * mask

def process_raw(raw, mrdHeader):
    rep_index = raw.acquisitions[0].idx.repetition

    # Preprocessing acquisitions
    # First repetition will contain a noise acq. Extract it and keep it for all reps.
    if raw.noise is None:
        raw.extract_noise()
    
    # Remove useless acquisitions at the beginning of each rep
    raw.remove_phase_stabilization_references()

    # Build kspace and navigator data structure
    print("Building kspace array...")
    kspace, navigator, acs_mask = raw.build_kspace()
    kspace = kspace[[rep_index], ...]           # (1, 1, 1, 1, nEcho=4, nSlice=15, 1, nKy=384, nKx=768, nCoils=4)
    navigator = navigator[[rep_index], ...]     # (1, 1, 1, 1, 1,       nSlice=15, 1, nKy=384, nKx=768, nCoils=4)
    acs_mask = acs_mask[[rep_index], ...]       # same as kspace

    # Save it for tests
    #raw.save_kspace("ice_data.npz")
    
    # Load precomputed kspace
    #raw.load_kspace("ice_data.npz")

    # --------------------------------------------------
    # Build reference volume and save it under NifTi
    CENTERLINE_DIR = Path("/workspaces/siemens-fire/Icesimu_output/sct_centerline")

    ref_path = build_reference_volume(
        kspace     = kspace,
        acs_mask   = acs_mask,
        mrdHeader  = mrdHeader,
        output_dir = CENTERLINE_DIR,
        raw = raw
    )

    print(f"Reference volume ready : {ref_path}")
    print("SCT centerline detection will now be run on this volume.")

    # Run SCT centerline detection — results saved to CENTERLINE_DIR for inspection
    # csv_path = run_sct_centerline(output_dir=CENTERLINE_DIR)
    csv_path = run_sct_centerline(
        output_dir = CENTERLINE_DIR,
        nKy        = kspace.shape[7],
        nKx_recon  = mrdHeader.encoding[0].reconSpace.matrixSize.x,
        nSlice     = kspace.shape[5]
    )

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
    # NAVIGATOR PREPARATION — same as original
    # S : (rep, samples=nKx, lines=nKy, slices=nSlice, coils)
    # --------------------------------------------------
    # # Extract phase for each repetition (TEMP for current data with multiple reps)
    # # TODO: This is not robust. Find a general way to deal with the multiple dims
    # S = navigator[:, 0, 0, 0, 0, :, 0, :, :, :] # shape : (rep, slices, lines, samples, coils)
    # S = np.transpose(S, (0, 3, 2, 1, 4)) # shape : (rep, samples, lines, slices, coils)

    S = navigator[:, 0, 0, 0, 0, :, 0, :, :, :]
    S = np.transpose(S, (0, 3, 2, 1, 4))

    # --------------------------------------------------
    # CENTERLINE MASKING on navigator lines
    # S[rep, samples, ky, sl, coil] — samples axis = nav readout
    # apply_nav_mask_from_centerline works on (nKx,) 1D line
    # center_x is in image coordinates (0 to nKx_recon)
    # --------------------------------------------------
    if use_mask:
        print("Applying centerline mask to navigator...")
        nKx_full  = S.shape[1]   # full readout size (with oversampling)
        nSlice_S  = S.shape[3]
        nRep_S    = S.shape[0]
        nKy_S     = S.shape[2]
        nCoils_S  = S.shape[4]

        S_masked = S.copy()
        for sl in range(nSlice_S):
            center_x = center_x_per_slice.get(sl)
            if center_x is None:
                print(f"  sl={sl} : no centerline — using full line")
                continue
            for rep in range(nRep_S):
                for ky in range(nKy_S):
                    for co in range(nCoils_S):
                        nav_line = S[rep, :, ky, sl, co]   # (nKx_full,)
                        S_masked[rep, :, ky, sl, co] = apply_nav_mask_from_centerline(
                            nav_line, center_x, nKx_recon  = mrdHeader.encoding[0].reconSpace.matrixSize.x, width=70
                        )
        S = S_masked
        print("Centerline masking applied.")

    # --------------------------------------------------
    # PHASE EXTRACTION — unchanged from original
    # --------------------------------------------------
    print("Computing corrections...")
    phase_extractions = np.stack([phase_extraction(s, raw.noise.data.T) for s in S])
    nx = mrdHeader.encoding[0].encodedSpace.matrixSize.x
    # les TE sont dans le header.
    # Par contre, le temps d'écho du navigateur n'est pas dans le header
    # Il serait supposé être dans les user_int des acqs, mais quand on passe
    # par FIRE, il est absent.
    # dt est dans le header de l'acquisition.
    echo_times = np.array(mrdHeader.sequenceParameters.TE, dtype=np.float32) * 1e-3
    print ("echo_times=", echo_times)
    navigator_te = 24e-3
    dt = 5e-6

    # Apply navigator correction
    print("Applying corrections...")
    field_estimates = field_conversion(phase_extractions, navigator_te) # rad/s
    corrected = kspace_correction(kspace, field_estimates, nx, echo_times, dt)

    # Use GRAPPA to fill in missing kspace lines
    print("GRAPPA correction...")
    corrected = grappa_reconstruction(corrected, corrected * acs_mask.squeeze())

    # Preprocessing before sending it back to ICE
    print("Reconstruction...")
    corrected = raw_to_image(corrected)
    corrected = mag_images(corrected)

    # Save images for faster testing
    # np.save("images.npy", images)

    field_of_view = (ctypes.c_float(mrdHeader.encoding[0].reconSpace.fieldOfView_mm.x), 
                            ctypes.c_float(mrdHeader.encoding[0].reconSpace.fieldOfView_mm.y), 
                            ctypes.c_float(mrdHeader.encoding[0].reconSpace.fieldOfView_mm.z))

    header_map = {}
    for acq in raw.acquisitions:
        key = (acq.idx.contrast, acq.idx.slice)
        if key not in header_map:
            header_map[key] = acq.getHead()
    acq_headers = [header_map[(c, s)] for c in range(4) for s in range(15)]
    for i, h in enumerate(acq_headers):
        print(f"Image {i}: contrast={h.idx.contrast}, slice={h.idx.slice}")
    ismrmrd_images = convert_to_ismrmrd_images(corrected, acq_headers, field_of_view)
    
    print("Sending images...")
    return ismrmrd_images

def mag_images(images):
    return np.abs(images).astype(np.float64)

def convert_float64_to_int16(data):
    return np.around(data * (2**12 - 1)/data.max()).astype(np.int16)

def convert_to_ismrmrd_images(images, acq_headers, fov):
    images_out = []
    *_, y, x = images.shape
    data = convert_float64_to_int16(images)
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

def raw_to_image(raw, noise_data=None):
    # assumed shape : (repetitions, echoes, slices, y, x, coils)
    #                 (0          , 1     , 2,    , 3, 4, 5) 

    data = np.flip(raw, (3, 4)) # inverser les données en x et y for some reason
    data = reconstruct_image(data)
    data *= np.prod(data.shape) # FFT scaling, for consistency with ICE apparently

    # Walsh coil combination
    print("  Using Inati coil combination with prewhitening...")
    data = coil_combination_Inati(data, noise_data=noise_data)

    # Remove readout oversampling by cropping
    data = remove_oversampling(data, 4)

    return data

def remove_oversampling(data, readout_axis, oversampling_factor=2):
    nx = data.shape[readout_axis]
    recon_size = nx // oversampling_factor
    start = (nx - recon_size) // 2

    return data.take(np.arange(start, start+recon_size), axis=readout_axis)

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
    rep, echo, slc, y, x, nCoils = data.shape

    # Prewhitening matrix from noise scan
    if noise_data is not None:
        dmtx = coils.calculate_prewhitening(noise_data)
        print(f"  Prewhitening matrix computed : shape={np.asarray(dmtx).shape}")
    else:
        dmtx = None

    combined = np.zeros((rep, echo, slc, y, x), dtype=np.float32)

    for r in range(rep):
        for e in range(echo):
            for s in range(slc):
                # Extract coil images : (y, x, nCoils) → (nCoils, y, x)
                img_coils = np.moveaxis(data[r, e, s], -1, 0)

                # Prewhitening
                if dmtx is not None:
                    img_coils = coils.apply_prewhitening(img_coils, dmtx)

                # Inati returns (csm, combined_image) directly
                _, combined_complex = coils.calculate_csm_inati_iter(
                    img_coils,
                    smoothing = smoothing,
                    niter     = niter,
                    thresh    = 1e-3
                )

                # Magnitude
                combined[r, e, s] = np.abs(combined_complex)

    return combined


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
    print("navigator.shape", navigator.shape)   # (768, 384, 15, coils=(4,8,...))
    print("noise.shape=", noise.shape)      # (768, coils=(4,8,...))

    # subtract first navigator phase to remove static phase contributions
    delta_S = navigator * np.exp(-1j*np.angle(navigator[:, [0], :, :]))  # (samples=768, lines=384, slices=15, coils=(4,8,...))
    print("delta_S.shape=", delta_S.shape)

    w = np.abs(delta_S) / np.std(noise, axis=0)   # (768, 384, 15, coils=(4,8,...))  (TODO: check if should need to raise to power 2)
    print("w.shape=", w.shape)
    # RuntimeWarning here because of dividing by zero
    w_tilde = w/np.sum(w, axis=(0, 3), keepdims=True) # (768, 384, 15, coils=(4,8,...))
    print("w_tilde.shape=", w_tilde.shape)
    # Replace resulting NaNs by zero
    w_tilde[np.isnan(w_tilde)] = 0.0

    delta_S = np.sum(w_tilde * delta_S, axis=(0, 3)) # (384, 15)  (lines, slices)
    print("delta_S.shape=", delta_S.shape)

    delta_phi_mean = np.angle(np.mean(delta_S, axis=0)) # (15,)  (slices,)
    print("delta_phi_mean.shape=", delta_phi_mean.shape)

    delta_S_tilde = delta_S * np.exp(-1j * delta_phi_mean) # (384, 15)  (lines, slices)
    print("delta_S_tilde.shape=", delta_S_tilde.shape) 

    delta_phi = np.angle(delta_S_tilde) # (384, 15)  (lines, slices)
    print("delta_phi.shape=", delta_phi.shape)

    return delta_phi

def field_conversion(nav_phases, te_nav):
    print ("nav_phases.shape=", nav_phases.shape) # (1, 384, 15)  (1, lines, slices)
    return nav_phases/te_nav    # rad/s

def kspace_correction(raw_data, field_estimates, n_samples, echo_times, dt):
    print("raw_data.shape", raw_data.shape)         # (1, 1, 1, 1, 4, 15, 1, 384, 768, coils=(4,8,...)) (1, 1, 1, 1, echo, slices, 1, lines, samples, coils)
    print("field_estimates.shape", field_estimates.shape)    # (1, 384, 15)  (1, lines, slices)
    print("echo_times.shape", echo_times.shape)     # (4,)  (echo,)

    t = np.array([(j - n_samples/2)*dt for j in range(n_samples)], dtype=np.float32)
    print("t.shape", t.shape) # (768,) (samples,)
    t = np.repeat(t[:, np.newaxis], echo_times.shape[0], axis=1)
    print("t.shape", t.shape) # (768, 4) (samples, echo)

    t += echo_times # will probably crash
    print("t.shape", t.shape) # (768, 4) (samples, echo)

    demodulation = np.exp(-1j * np.einsum('rlp,je->replj', field_estimates, t))
    print("demodulation.shape", demodulation.shape) # (1, 4, 15, 384, 768) (1, echo, slices, lines, samples)

    # TODO: handle all the dimensions correctly instead of squeezing
    corrected = raw_data.squeeze()*demodulation[..., np.newaxis] # (1, 4, 15, 384, 768, coils=(4,8,...)) (1, echo, slices, lines, samples, coils)
    print("corrected.shape", corrected.shape)
    return corrected
 
def grappa_reconstruction(kspace, acs_lines, R=2, kernel_size=(5, 5)):
    
    *leading, y, x, c = kspace.shape

    kspace_r   = kspace.reshape(-1, y, x, c)
    acs_mask_r = acs_lines.reshape(-1, y, x, c)

    results = np.memmap(
        os.path.join(mkdtemp(), 'grappa.dat'),
        dtype=np.complex64, mode='w+',
        shape=kspace_r.shape
    )

    for i in range(kspace_r.shape[0]):
        # kspace_r[i]   : (nKy, nKx, nCoils) → reorder to (nCoils, nKy, nKx)
        k_coils   = np.moveaxis(kspace_r[i],   -1, 0)   # (nCoils, nKy, nKx)
        acs_coils = np.moveaxis(acs_mask_r[i], -1, 0)   # (nCoils, nKy, nKx)

        # Extract ACS line indices from mask
        # acs_coils[0] is boolean mask of shape (nKy, nKx)
        # ACS lines are those where the entire ky line is non-zero
        acs_line_mask = np.any(np.abs(acs_coils[0]) > 1e-8, axis=-1)
        acs_idx       = np.where(acs_line_mask)[0]   # (nACS,)

        if len(acs_idx) == 0:
            print(f"  WARNING : no ACS lines found for batch {i}")
            results[i] = kspace_r[i]
            continue

        # Apply GRAPPA — returns (nCoils, nKy, nKx)
        k_filled = apply_grappa(k_coils, acs_idx, R=R, kernel_size=kernel_size)

        # Reorder back to (nKy, nKx, nCoils)
        results[i] = np.moveaxis(k_filled, 0, -1)

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