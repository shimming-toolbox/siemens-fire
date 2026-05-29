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

from common import SiemensRAW

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

    # Save it for tests
    #raw.save_kspace("ice_data.npz")
    
    # Load precomputed kspace
    #raw.load_kspace("ice_data.npz")

    # Extract phase for each repetition (TEMP for current data with multiple reps)
    # TODO: This is not robust. Find a general way to deal with the multiple dims
    S = navigator[:, 0, 0, 0, 0, :, 0, :, :, :] # shape : (rep, slices, lines, samples, coils)
    S = np.transpose(S, (0, 3, 2, 1, 4)) # shape : (rep, samples, lines, slices, coils)

    print("Computing corrections...")
    phase_extractions = np.stack([phase_extraction(s, raw.noise.data.T) for s in S])
    nx = mrdHeader.encoding[0].encodedSpace.matrixSize.x
    # les TE sont dans le header.
    # Par contre, le temps d'écho du navigateur n'est pas dans le header
    # Il serait supposé être dans les user_int des acqs, mais quand on passe
    # par FIRE, il est absent.
    # dt est dans le header de l'acquisition.
    echo_times = np.array(mrdHeader.sequenceParameters.TE, dtype=np.float32) * 1e-3
    navigator_te = 24e-3
    dt = 5e-6

    # Apply navigator correction
    print("Applying corrections...")
    field_estimates = field_conversion(phase_extractions, navigator_te)
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

    ismrmrd_images = convert_to_ismrmrd_images(corrected, raw.acquisitions[0].getHead(), field_of_view)
    print("Sending images...")
    return ismrmrd_images

def mag_images(images):
    return np.abs(images).astype(np.float64)

def convert_float64_to_int16(data):
    return np.around(data * (2**12 - 1)/data.max()).astype(np.int16)

def convert_to_ismrmrd_images(images, acq_header, fov):
    images_out = []
    
    *_, y, x = images.shape

    data = convert_float64_to_int16(images)

    for i, img in enumerate(data.reshape(-1, y, x)):
        ismrmrd_image = ismrmrd.Image.from_array(img, transpose=False)
        ismrmrd_image.setHead(mrdhelper.update_img_header_from_raw(ismrmrd_image.getHead(), acq_header))
        ismrmrd_image.field_of_view = fov
        ismrmrd_image.image_index = i

        # Set ISMRMRD Meta Attributes
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
    np.fft.ifft2(image, axes=axes, out=image)
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

def grappa_reconstruction(kspace, acs_lines):

    def _grappa(kspace, acs_lines):
        calib = np.trim_zeros(acs_lines)
        # transpose, because pygrappa assumes (x, y, c) and not (y, x, c) like we have
        calib = np.transpose(calib, (1, 0, 2))
        k = np.transpose(kspace, (1, 0, 2))
        # tranpose back to keep it consistent with the rest of our code
        return np.transpose(grappa(k, calib), (1, 0, 2))

    *leading, y, x, c = kspace.shape

    # (-1) in reshape means it's inferred from the other dims
    kspace_r = kspace.reshape(-1, y, x, c)
    acs_lines_r = acs_lines.reshape(-1, y, x, c)

    results = np.memmap(os.path.join(mkdtemp(), 'grappa.dat'), dtype=np.complex64, mode='w+', shape=kspace_r.shape)
    for i in range(kspace_r.shape[0]):
        results[i] = _grappa(kspace_r[i], acs_lines_r[i])
    results.flush()

    return results.reshape(*leading, y, x, c)
