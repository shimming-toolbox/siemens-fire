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

        phase_extractions = np.concat([phase_extraction(s, raw.noise.data.T) for s in S])
        np.save("phases.npy", phase_extractions)


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
