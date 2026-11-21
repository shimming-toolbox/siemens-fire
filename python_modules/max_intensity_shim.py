# Adapted from Kelvin Chow's python-ismrmrd-server repository
# Source: https://github.com/kspaceKelvin/python-ismrmrd-server

from datetime import datetime
from multiprocessing import connection
import ismrmrd
import os
import logging
import traceback
import json
import numpy as np
import xml.dom.minidom
import base64
import mrdhelper # Custom module for MRD helper functions found in the python-ismrmrd-server repository
import constants # Custom module for constants found in the python-ismrmrd-server repository
import subprocess
import nibabel as nib
import random
import shutil


# Folder for debug output files
debugFolder = "/tmp/share/debug"
shimFolder = f"{debugFolder}/max_intensity_shim"
mrd2niiFolder = f"{shimFolder}/mrd2nii"
dataFolder = "/tmp/share/saved_data"

def process(connection, config, mrdHeader):
    logging.info("Config: \n%s", config)

    global shimFolder, mrd2niiFolder
    shimFolder = f"{debugFolder}/max_intensity_shim" + datetime.now().strftime("_%Y-%m-%d-%H%M%S" + "_" + str(random.randint(0,100)))
    mrd2niiFolder = f"{shimFolder}/mrd2nii"

    # mrdHeader should be xml formatted MRD header, but may be a string
    # if it failed conversion earlier
    try:
        # Disabled due to incompatibility between PyXB and Python 3.8:
        # https://github.com/pabigot/pyxb/issues/123
        # # logging.info("MRD header: \n%s", mrdHeader.toxml('utf-8'))

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

    imgGroup = []
    dset = None

    try:
        for item in connection:
            if isinstance(item, ismrmrd.Image):
                imgGroup.append(item)

    except Exception as e:
        logging.error(traceback.format_exc())
        connection.send_logging(constants.MRD_LOGGING_ERROR, traceback.format_exc())

    finally:
        process_acquisition(imgGroup, connection, config, mrdHeader, dset)
        connection.send_close()

def process_acquisition(imgGroup, connection, config, mrdHeader, dset):
    if isinstance(config, str):
        config_filename = f"{config}.json"
        config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), config_filename)
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"{config_path} does not exist.")

        with open(config_path, "r") as f:
            config_dict = json.load(f)
    elif isinstance(config, dict):
        config_dict = config
    else:
        raise RuntimeError("Config should be a string or a dictionary")

    # Create folder, if necessary
    if not os.path.exists(debugFolder):
        os.makedirs(debugFolder)
        logging.debug("Created folder " + debugFolder + " for debug output files")

    if os.path.exists(shimFolder):
        shutil.rmtree(shimFolder)
        logging.debug("Removed existing folder " + shimFolder + " for max intensity shim output files")

    os.makedirs(shimFolder)
    logging.debug("Created folder " + shimFolder + " for max intensity shim output files")
    
    if not os.path.exists(mrd2niiFolder):
        os.makedirs(mrd2niiFolder)
        logging.debug("Created folder " + mrd2niiFolder + " for MRD max intensity shim output files")
    subprocess.run(['mrd2nii',
                    '-i', connection.mrdFilePath,
                    '-o', mrd2niiFolder],
                    check=True)

    fnames = os.listdir(mrd2niiFolder)

    fname_mag = ""
    for fname in fnames:
        if fname.endswith("_magnitude_echo-1.nii.gz"):
            fname_mag = os.path.join(mrd2niiFolder, fname)
            logging.info("Found magnitude image: %s", fname_mag)
            break
    if fname_mag == "":
        raise FileNotFoundError("Could not find magnitude image in converted MRD data")
    
    use_already_computed_mask = mrdhelper.get_json_config_param(config_dict, 'use_already_computed_mask', default=False, type='bool')
    # TEMP: TODO REMOVE
    # use_already_computed_mask = True
    
    if use_already_computed_mask:
        fname_mask = os.path.join(dataFolder, "mask.nii.gz")
        if not os.path.exists(fname_mask):
            raise FileNotFoundError("Could not find mask image: " + fname_mask)

    else:
        nii_mag = nib.load(fname_mag)
        if nii_mag.ndim != 4:
            raise RuntimeError("Magnitude image is not 4D as expected for max intensity shim")
        mag_mean = np.mean(nii_mag.get_fdata(), axis=(3))
        nii_mag_mean = nib.Nifti1Image(mag_mean, nii_mag.affine, nii_mag.header)
        nii_mag_mean.set_data_dtype(np.float32)
        fname_mag_mean = os.path.join(shimFolder, "magnitude_mean.nii.gz")
        nib.save(nii_mag_mean, fname_mag_mean)
        logging.info("Saved mean magnitude image over time to: %s", fname_mag_mean)

        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = "0"
        env["SCT_USE_GPU"] = "1"
        path_sct_binaries = '/opt/code/spinalcordtoolbox/bin'
        fname_mask = os.path.join(shimFolder, 'mask_sct_deepseg_seg.nii.gz')
        subprocess.run([os.path.join(path_sct_binaries, 'sct_deepseg'), 'spinalcord',
                        '-i', fname_mag_mean,
                        '-o', fname_mask],
                        env=env,
                        check=True)

    fname_output = os.path.join(shimFolder, "scanner_shim.txt")
    fname_applied_shims = os.path.join(dataFolder, "scanner_shim_last_acq.txt")
    try:
        ret = subprocess.run(['st_b0shim', 'max-intensity',
                        '-i', fname_mag,
                        '--coefs', fname_applied_shims,
                        '--mask', fname_mask,
                        '-o', shimFolder],
                        capture_output=True,
                        text=True,
                        check=True)
        logging.info("Command stdout: %s", ret.stdout)
        logging.info("Command stderr: %s", ret.stderr)
    except subprocess.CalledProcessError as e:
        logging.error("st_b0shim failed with return code %d", e.returncode)
        logging.error("Command stdout: %s", e.stdout)
        logging.error("Command stderr: %s", e.stderr)
        if os.path.exists(fname_output):
            os.remove(fname_output)
        return

    shutil.copyfile(fname_output, os.path.join(dataFolder, "scanner_shim.txt"))

    return
