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
    raw = SiemensRAW(mrdHeader)
    try:
        for item in connection:
            # ----------------------------------------------------------
            # Raw k-space data messages
            # ----------------------------------------------------------
            if isinstance(item, ismrmrd.Acquisition):
                raw.add_acq(item)

            elif item is None:
                break

            else:
                logging.error("Unsupported data type %s", type(item).__name__)

        # Preprocessing acquisitions
        raw.extract_noise()
        raw.remove_phase_stabilization_references()

        # Build kspace and navigator data structure
        raw.build_kspace()

        # Save it for tests
        raw.save_kspace("ice_data.npz")

        # Extract phase for each repetition (TEMP for current data with multiple reps)
        S = raw.navigator[:, 0, 0, 0, 0, :, 0, :, :, :] # shape : (rep, slices, lines, samples, coils)

        S = np.transpose(S, (0, 3, 2, 1, 4)) # shape : (rep, samples, lines, slices, coils)

        phase_extractions = np.stack([phase_extraction(s, raw.noise.data.T) for s in S])
        nx = mrdHeader.encoding[0].encodedSpace.matrixSize.x

        # À coder autrement qu'en dur comme ça. dt est dans le header
        # d'acq, les temps d'echo disparaissent dans FIRE par contre.
        # Cependant, ils sont dans le header MRD.
        echo_times = np.array([6900e-6, 10920e-6, 14940e-6, 18960e-6])
        navigator_te = 24e-3
        dt = 5e-6

        corrected_raw = kspace_correction(raw.kspace, field_conversion(phase_extractions, navigator_te), nx, echo_times, dt)

        corrected_images = grappa_reconstruction(corrected_raw, corrected_raw * raw.acs_mask.squeeze())

    except Exception as e:
        logging.error(traceback.format_exc())
        connection.send_logging(constants.MRD_LOGGING_ERROR, traceback.format_exc())

    finally:
        connection.send_close()

def phase_extraction(navigator, noise):
    """
    
    Arguments:
    navigator -- shape : (j, l, p, c) -> (samples, lines, slices, coils)
    noise     -- shape : (j, c)       -> (samples, coils)
    """
    delta_S = navigator * np.exp(-1j*np.angle(navigator[:, [0], :, :]))

    w = np.abs(delta_S) / np.std(noise, axis=0) # TODO: check if should need to raise to power 2
    w_tilde = w/np.sum(w, axis=(0, 3), keepdims=True)
    w_tilde[np.isnan(w_tilde)] = 0.0

    delta_S = np.sum(w_tilde * delta_S, axis=(0, 3))

    delta_phi_mean = np.angle(np.mean(delta_S, axis=0))

    delta_S_tilde = delta_S * np.exp(-1j * delta_phi_mean)

    delta_phi = np.angle(delta_S_tilde)

    return delta_phi

def field_conversion(nav_phases, te_nav):
    return nav_phases/te_nav

def kspace_correction(raw_data, field_estimates, n_samples, echo_times, dt):
    t = np.array([(j - n_samples/2)*dt for j in range(n_samples)])
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

    results = np.stack([_grappa(kspace_r[i], acs_lines_r[i]) for i in range(kspace_r.shape[0])])

    return results.reshape(*leading, y, x, c)
