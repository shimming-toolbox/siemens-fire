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
b0shimFolder = f"{debugFolder}/st_b0shim"
mrd2niiFolder = f"{b0shimFolder}/mrd2nii"
dataFolder = "/tmp/share/saved_data"

def process(connection, config, mrdHeader):
    logging.info("Config: \n%s", config)

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

    if os.path.exists(b0shimFolder):
        shutil.rmtree(b0shimFolder)
        logging.debug("Removed existing folder " + b0shimFolder + " for B0 shim output files")

    os.makedirs(b0shimFolder)
    logging.debug("Created folder " + b0shimFolder + " for B0 shim output files")
    
    if not os.path.exists(mrd2niiFolder):
        os.makedirs(mrd2niiFolder)
        logging.debug("Created folder " + mrd2niiFolder + " for MRD b0shim output files")

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

    fname_fmap = os.path.join(dataFolder, "fieldmap.nii.gz")
    if not os.path.exists(fname_fmap):
        raise FileNotFoundError("Could not find fieldmap image: " + fname_fmap)
    
    fname_mask = os.path.join(dataFolder, "mask.nii.gz")
    if not os.path.exists(fname_mask):
        raise FileNotFoundError("Could not find mask image: " + fname_mask)
    
    slices_option = mrdhelper.get_json_config_param(config_dict, 'shim_method', default='auto', type='str')
    if slices_option not in ['auto', 'slicewise', 'volume']:
        raise ValueError("Invalid shim method: " + slices_option)
    if slices_option == 'slicewise':
        slices_option = 'auto'

    subprocess.run(['st_b0shim', 'dynamic',
                    '--fmap', fname_fmap,
                    '--target', fname_mag,
                    '--mask', fname_mask,
                    '--slices', slices_option,
                    '--scanner-coil-order', mrdhelper.get_json_config_param(config_dict, 'scanner-coil-order', default='0,1', type='str'),
                    '--optimizer-method', mrdhelper.get_json_config_param(config_dict, 'optimizer-method', default='pseudo_inverse', type='str'),
                    '--optimizer-criteria', mrdhelper.get_json_config_param(config_dict, 'optimizer-criteria', default='rmse', type='str'),
                    '--regularization-factor', str(mrdhelper.get_json_config_param(config_dict, 'regularization-factor', default=0.1, type='float')),
                    '--weighting-signal-loss', str(mrdhelper.get_json_config_param(config_dict, 'weighting-signal-loss', default=10, type='float')),
                    '--mask-dilation-kernel-size', str(mrdhelper.get_json_config_param(config_dict, 'mask-dilation-kernel-size', default=5, type='int')),
                    '--output-file-format-scanner', 'slicewise-hrd',
                    '-o', b0shimFolder],
                    check=True)

    shutil.copyfile(os.path.join(b0shimFolder, "scanner_shim.txt"), os.path.join(dataFolder, "scanner_shim.txt"))

    return


def extract_nii_slice_ordering_to_chronological(json_data, n_slices):
    slice_order_nii_to_chrono = {}
    slice_timing = json_data.get("SliceTiming")
    if slice_timing is not None and n_slices == 1:
        return {0: 0}
    indexes_chr_to_nii = np.argsort(slice_timing)
    indexes_nii_to_chr = np.argsort(indexes_chr_to_nii)
    for i_slice in range(n_slices):
        slice_order_nii_to_chrono[i_slice] = indexes_nii_to_chr[i_slice]
    return slice_order_nii_to_chrono


def create_debug_save_file():
    # Create savedata folder, if necessary
    if not os.path.exists(debugFolder):
        os.makedirs(debugFolder)

    mrdFilePath = os.path.join(debugFolder, "MRD_output_" + datetime.now().strftime("%Y-%m-%d-%H%M%S" + "_" + str(random.randint(0,100)) + ".h5"))

    # Create HDF5 file to store incoming MRD data
    logging.info("Incoming data will be saved to: '%s' in group '%s'", mrdFilePath, "dataset")
    dset = ismrmrd.Dataset(mrdFilePath, "dataset")
    return dset
