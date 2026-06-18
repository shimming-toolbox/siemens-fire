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
finstshimFolder = f"{debugFolder}/finst_shim"
mrd2niiFolder = f"{finstshimFolder}/mrd2nii"
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

    if os.path.exists(finstshimFolder):
        shutil.rmtree(finstshimFolder)
        logging.debug("Removed existing folder " + finstshimFolder + " for finst shim output files")

    os.makedirs(finstshimFolder)
    logging.debug("Created folder " + finstshimFolder + " for finst shim output files")
    
    if not os.path.exists(mrd2niiFolder):
        os.makedirs(mrd2niiFolder)
        logging.debug("Created folder " + mrd2niiFolder + " for MRD finst shim output files")
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
            raise RuntimeError("Magnitude image is not 4D as expected for finst_shim")
        mag_mean = np.mean(nii_mag.get_fdata(), axis=(3))
        nii_mag_mean = nib.Nifti1Image(mag_mean, nii_mag.affine, nii_mag.header)
        fname_mag_mean = os.path.join(finstshimFolder, "magnitude_mean.nii.gz")
        nib.save(nii_mag_mean, fname_mag_mean)
        logging.info("Saved mean magnitude image over time to: %s", fname_mag_mean)

        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = "0"
        env["SCT_USE_GPU"] = "1"
        path_sct_binaries = '/opt/code/spinalcordtoolbox/bin'
        fname_mask = os.path.join(finstshimFolder, 'mask_sct_deepseg_seg.nii.gz')
        subprocess.run([os.path.join(path_sct_binaries, 'sct_deepseg'), 'spinalcord',
                        '-i', fname_mag_mean,
                        '-o', fname_mask],
                        env=env,
                        check=True)

    fname_best_idx = os.path.join(finstshimFolder, "best_shim_index.txt")
    subprocess.run(['st_b0shim', 'max-intensity',
                    '-i', fname_mag,
                    '--mask', fname_mask,
                    '-o', fname_best_idx],
                    check=True)

    fname_output = os.path.join(finstshimFolder, "scanner_shim.txt")
    fname_applied_shims = os.path.join(dataFolder, "scanner_shim_last_acq.txt")
    best_idx_to_scanner_shim(fname_applied_shims, fname_best_idx, fname_output)

    shutil.copyfile(fname_output, os.path.join(dataFolder, "scanner_shim.txt"))

    return


def best_idx_to_scanner_shim(fname_shim_file, fname_best_indexes, fname_output):

    with open(fname_best_indexes, 'r') as f:
        lines = f.readlines()

    n_slices = int(lines[0])
    best_indexes = [int(x) for x in lines[1].split(' ')]

    with open(fname_shim_file, 'r') as f:
        lines = f.readlines()

    coefs = []
    for i_slice in range(n_slices):
        best_line = lines[(best_indexes[i_slice] * n_slices) + 1]
        values = [float(x.strip()) for x in best_line.split('|')]
        coefs.append(values)

    # Create column names: "orderX_channelY"
    column_names = ['f0', 'Gx', 'Gy', 'Gz']
    coefs_array = np.array(coefs)

    # Compute column widths
    # 1. Get max formatted value length in each column
    formatted_values = []
    col_widths = []

    for col_idx in range(coefs_array.shape[1]):
        col_vals = coefs_array[:, col_idx]
        formatted_col = [f"{val:.6f}" for val in col_vals]
        formatted_values.append(formatted_col)
        max_val_len = max(len(s) for s in formatted_col)
        col_name_len = len(column_names[col_idx])
        col_width = max(max_val_len, col_name_len)
        col_widths.append(col_width)

    # Write to file manually
    with open(fname_output, 'w') as f:
        # Write header (centered titles)
        header_cells = [column_names[i].center(col_widths[i]) for i in range(len(column_names))]
        header = ' | '.join(header_cells)
        f.write(header + '\n')

        # Write each row of shim values (right-aligned)
        nb_rows = coefs_array.shape[0]
        for row_idx in range(nb_rows):
            row_cells = [
                formatted_values[col_idx][row_idx].rjust(col_widths[col_idx])
                for col_idx in range(len(column_names))
            ]
            row_str = ' | '.join(row_cells)
            f.write(row_str + '\n')
