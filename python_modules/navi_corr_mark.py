import ismrmrd
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
from pygrappa import grappa
from tempfile import mkdtemp
from sklearn import cluster as KMeans
from grappa import grappa
from SLR import SLR
from common import SiemensRAW

'''
This script is a modified version of navi_corr.py found on the navi-me-gre branch of siemens-fire. Parts of the process_raw function were changed in order to implement the a new reconstruction method.
The method bins navigators and seperates the lines of kspace into these bins with the goal to seperate lines of kspace based on when they were acquired in the breathing cycle. The original recon script can be found in this repo,
it's called Example_Recon.ipynb.
'''
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
    kspace = kspace[[rep_index], ...]
    navigator = navigator[[rep_index], ...]
    acs_mask = acs_mask[[rep_index], ...]

    # Whiten noise
    print("Whitening noise...")
    noise_data = raw.noise.data.T
    cov = np.cov(noise_data, rowvar=False)
    inv_sq_cov = np.linalg.cholesky(np.linalg.inv(cov)).T.conj()

    kspace = np.einsum('...c,cc->...c', kspace, inv_sq_cov)
    navigator = np.einsum('...c,cc->...c', navigator, inv_sq_cov)


    num_shots = navigator.shape[3]  # ny dimension
    nav_features = navigator[0, 0, 0, :, :, :].reshape(num_shots, -1)
    nav_features = np.abs(nav_features)
    
    # k-means clustering to assign each shot to a motion bin
    kmeans = KMeans(n_clusters=8, random_state=42, n_init=10).fit(nav_features)
    idx = kmeans.labels_

    _, neco, _, ny, nx, nc = kspace.shape
    nbins = 8
    R = 2
    
    dat = np.zeros((nx, ny, nbins, neco, nc), dtype='complex64')
    cnt = np.zeros((nx, ny, nbins, neco, nc)) 

    # Distribute the data lines for this single repetition
    for j in range(ny):
        # Check if this line is an active k-space line or part of the ACS auto-calibration region
        # We look at any echo/coil just to see if data was collected on line 'j'
        is_line_acquired = np.any(np.abs(kspace[0, 0, 0, j, :, 0]) > 0)
        is_acs_line = np.any(acs_mask[0, 0, 0, j, :, 0] > 0)
        
        # Grab the motion bin assigned to this line
        current_bin = idx[j]
        
        if is_line_acquired or is_acs_line:
            # Extract the line (shape: echoes, nx, coils) -> transpose to (nx, echoes, coils)
            line_data = kspace[0, :, 0, j, :, :].transpose((1, 0, 2))
            
            # Place it into the target matrix
            dat[:, j, current_bin, :, :] += line_data
            cnt[:, j, current_bin, :, :] += 1

    # Average out overlapping lines safely
    mask_nonzero = (cnt != 0)
    dat[mask_nonzero] = dat[mask_nonzero] / cnt[mask_nonzero]

    init = np.zeros((nc, nx, ny, neco), dtype=np.complex64)

    for echo in range(neco):
        # 1. Pull the accelerated kspace for this echo (ny, nx, nc) -> Transpose to (nc, nx, ny)
        # Note: Using your target script's transpose structure (2, 1, 0)
        k_input = kspace[0, echo, 0, :, :, :].transpose((2, 1, 0))
        
        # 2. Pull the full ACS auto-calibration region for this echo
        # We isolate where the acs_mask is True along the ny direction
        acs_indices = np.where(np.any(acs_mask[0, echo, 0, :, :, 0], axis=1))[0]
        ref_input = kspace[0, echo, 0, acs_indices, :, :].transpose((2, 1, 0))
        
        # 3. Call your custom grappa library function
        # Passing: input_data, reference_data, (lamda/reg), (kernel_x, kernel_y)
        init[:, :, :, echo] = grappa.grappa(k_input, ref_input, (1, 2), (3, 2))

    # Copy the GRAPPA reconstruction across the 8 bins (Just like your script)
    init = np.tile(init[:, :, :, None, :].transpose((1, 2, 3, 4, 0)), (1, 1, nbins, 1, 1))

    # SLR recon
    r = 150

    # set number of iterations
    niters = 100

    # slr kernel size
    kernel = (5,5)
    out = SLR.ADMM(dat, SLR.c_matrix, kernel, r, niters=niters, init=init)

    # Preprocessing before sending it back to ICE
    img_space = ifftdim(init.reshape((nx, ny, nbins, neco, nc)), dims=(0,1))
    img_transposed = img_space.transpose((0, 1, 3, 2, 4))
    img_reshaped = img_transposed.reshape((nx, ny, neco, -1))
    mag_combined = sos(img_reshaped, axis=-1)
    corrected = mag_combined.transpose((2, 1, 0))

    corrected = corrected[None, :, None, :, :]
    corrected = remove_oversampling(corrected, readout_axis=4)
    corrected = corrected.astype(np.float64)

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

def raw_to_image(raw):
    # assumed shape : (repetitions, echoes, slices, y, x, coils)
    #                 (0          , 1     , 2,    , 3, 4, 5) 

    data = np.flip(raw, (3, 4)) # inverser les données en x et y for some reason
    data = reconstruct_image(data)
    data *= np.prod(data.shape) # FFT scaling, for consistency with ICE apparently

    # RMS (temporary)
    data = coil_combination(data)

    # Remove readout oversampling by cropping
    data = remove_oversampling(data, 4)

    return data

def remove_oversampling(data, readout_axis, oversampling_factor=2):
    nx = data.shape[readout_axis]
    recon_size = nx // oversampling_factor
    start = (nx - recon_size) // 2

    return data.take(np.arange(start, start+recon_size), axis=readout_axis)

def coil_combination(data, coil_axis=-1):
    return np.sqrt(np.sum(np.square(data), axis=coil_axis))


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
    delta_S = navigator * np.exp(-1j*np.angle(navigator[:, [0], :, :])) 

    w = np.abs(delta_S) / np.std(noise, axis=0) # TODO: check if should need to raise to power 2
    # RuntimeWarning here because of dividing by zero
    w_tilde = w/np.sum(w, axis=(0, 3), keepdims=True)
    # Replace resulting NaNs by zero
    w_tilde[np.isnan(w_tilde)] = 0.0

    delta_S = np.sum(w_tilde * delta_S, axis=(0, 3))

    delta_phi_mean = np.angle(np.mean(delta_S, axis=0))

    delta_S_tilde = delta_S * np.exp(-1j * delta_phi_mean)

    delta_phi = np.angle(delta_S_tilde)

    return delta_phi

def field_conversion(nav_phases, te_nav):
    return nav_phases/te_nav

def kspace_correction(raw_data, field_estimates, n_samples, echo_times, dt):
    t = np.array([(j - n_samples/2)*dt for j in range(n_samples)], dtype=np.float32)
    t = np.repeat(t[:, np.newaxis], echo_times.shape[0], axis=1)

    t += echo_times # will probably crash

    demodulation = np.exp(-1j * np.einsum('rlp,je->replj', field_estimates, t))

    # TODO: handle all the dimensions correctly instead of squeezing
    return raw_data.squeeze()*demodulation[..., np.newaxis] # again check shape for broadcasting


# handles centric k-space with shifting, along arbitrary dimensions
def fftdim(x, dims=None):
    return np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(x), axes=dims, norm="ortho"))

def ifftdim(x, dims=None):
    return np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(x), axes=dims, norm="ortho"))

# root-sum-of-squares combination
def sos(x, axis=-1):
    return np.sqrt(np.sum(np.abs(x)**2, axis=axis))
