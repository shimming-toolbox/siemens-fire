# Adapted from Kelvin Chow's python-ismrmrd-server repository
# Source: https://github.com/kspaceKelvin/python-ismrmrd-server

from datetime import datetime
from multiprocessing import connection
import h5py
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
from mrd2nii.mrd2nii_main import get_main_dir


# Folder for debug output files
debugFolder = "/tmp/share/debug"
fmapFolder = f"{debugFolder}/st_fieldmap"
mrd2niiFolder = f"{fmapFolder}/mrd2nii"
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
        if dset is None:
            dset = create_debug_save_file()
        images = process_acquisition(imgGroup, connection, config, mrdHeader, dset)
        for image in images:
            connection.send_image(image)

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

    if os.path.exists(fmapFolder):
        shutil.rmtree(fmapFolder)
        logging.debug("Removed existing folder " + fmapFolder + " for fieldmap output files")

    os.makedirs(fmapFolder)
    logging.debug("Created folder " + fmapFolder + " for fieldmap output files")
    
    if not os.path.exists(mrd2niiFolder):
        os.makedirs(mrd2niiFolder)
        logging.debug("Created folder " + mrd2niiFolder + " for fieldmap output files")

    use_mask = mrdhelper.get_json_config_param(config_dict, 'use_mask', default=True, type='bool')
    fname_mask = os.path.join(dataFolder, "mask.nii.gz")
    if use_mask:
        if not os.path.exists(fname_mask):
            raise FileNotFoundError("Could not find mask image: " + fname_mask)

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
        raise FileNotFoundError("Could not find magnitude image with suffix '_magnitude_echo-1.nii.gz' in " + mrd2niiFolder)
    
    fnames_phases = []
    for fname in fnames:
        if fname.endswith("_phase_echo-1.nii.gz"):
            fnames_phases.append(os.path.join(mrd2niiFolder, fname))
            logging.info("Found phase image: %s", os.path.join(mrd2niiFolder, fname))
        if fname.endswith("_phase_echo-2.nii.gz"):
            fnames_phases.append(os.path.join(mrd2niiFolder, fname))
            logging.info("Found phase image: %s", os.path.join(mrd2niiFolder, fname))
    fnames_phases = sorted(fnames_phases)

    if len(fnames_phases) < 2:
        raise FileNotFoundError("Could not find at least two phase images with suffix '_phase_echo-1.nii.gz' and '_phase_echo-2.nii.gz' in " + mrd2niiFolder)

    fname_mask_saved = os.path.join(fmapFolder, "saved_mask.nii.gz")
    fname_fmap = os.path.join(fmapFolder, "fieldmap.nii.gz")
    env = os.environ.copy()
    env["FSLOUTPUTTYPE"] = "NIFTI_GZ"
    cmd = ['st_prepare_fieldmap',
           *fnames_phases,
           '--mag', fname_mag,
           '--unwrapper', mrdhelper.get_json_config_param(config_dict, 'unwrapper', default='prelude'),
           '--gaussian-filter', mrdhelper.get_json_config_param(config_dict, 'gaussian-filter', default='true'),
           '--sigma', mrdhelper.get_json_config_param(config_dict, 'sigma', default='1'),
            # '--2d', mrdhelper.get_json_config_param(config_dict, '2d', default='false'),
            '--savemask', fname_mask_saved,
            '-o', fname_fmap]

    if use_mask:
        logging.info("Using mask: %s", fname_mask)
        dilate_mask = mrdhelper.get_json_config_param(config_dict, 'dilate_mask', default=True, type='bool')
        if dilate_mask:
            dilate_mask_kernel_size = mrdhelper.get_json_config_param(config_dict, 'dilate_mask_kernel_size', default='3', type='str')
            fname_mask_dilated = os.path.join(fmapFolder, "mask_dilated.nii.gz")
            cmd_dilate = ['st_mask', 'modify-binary-mask', '-i', fname_mask, '-o', fname_mask_dilated, '--operation', 'dilate', '--size', dilate_mask_kernel_size, '--shape', 'sphere']
            subprocess.run(cmd_dilate,
                        env=env,
                        capture_output=True,
                        text=True,
                        check=True)
            cmd += ['--mask', fname_mask_dilated]
        else:
            cmd += ['--mask', fname_mask]
    else:
        logging.info("Using thresholding instead of mask")
        thr = mrdhelper.get_json_config_param(config_dict, 'threshold', default='0.1')
        cmd += ['--threshold', thr]

    subprocess.run(cmd,
                   env=env,
                   check=True)

    nii_fmap = nib.load(fname_fmap)
    data_fmap = nii_fmap.get_fdata().astype(np.int16)

    # Scale to int16 max
    # dyn_range = data_fmap.max() - data_fmap.min()
    # if dyn_range == 0:
    #    dyn_range = 1
    # data_fmap = (data_fmap - data_fmap.min()) / dyn_range
    # data_fmap *= 32767
    # data_fmap = data_fmap.astype(np.int16)

    fname_fmap_json = fname_fmap.replace('.nii.gz', '.json')
    with open(fname_fmap_json, 'r') as f:
        sidecar = json.load(f)

    # Copy to main data directory
    shutil.copyfile(fname_fmap_json, os.path.join(dataFolder, "fieldmap.json"))
    shutil.copyfile(fname_fmap, os.path.join(dataFolder, "fieldmap.nii.gz"))

    head = [img.getHead()                                  for img in imgGroup]
    meta = [ismrmrd.Meta.deserialize(img.attribute_string) for img in imgGroup]
    
    mag1_head = []
    mag1_meta = []
    magGroup = []

    # Extract first echo
    for i, h in enumerate(head):
        if h.image_type == ismrmrd.IMTYPE_MAGNITUDE and h.contrast == 0:
            mag1_head.append(h)
            mag1_meta.append(meta[i])
            magGroup.append(imgGroup[i])

    # Display MetaAttributes for first image
    logging.debug("MetaAttributes[0]: %s", ismrmrd.Meta.serialize(mag1_meta[0]))

    # Optional serialization of ICE MiniHeader
    if 'IceMiniHead' in mag1_meta[0]:
        logging.debug("IceMiniHead[0]: %s", base64.b64decode(mag1_meta[0]['IceMiniHead']).decode('utf-8'))

    currentSeries = 0
    
    # Find slice encoding direction
    dim_info = nii_fmap.header.get_dim_info()

    if None in dim_info:
        nii = nib.load(fname_mag)
        if np.allclose(nii.affine, nii_fmap.affine) and nii.shape == nii_fmap.shape:
            logging.warning("Input image and output mask have the same orientation and shape, assuming slice encoding direction is the same")
            dim_info = nii.header.get_dim_info()

    dim_of_freq_phase_slice_enc_directions = np.argsort(dim_info)

    if len(nii_fmap.shape) == 3:
        nb_slices = nii_fmap.shape[dim_of_freq_phase_slice_enc_directions[2]]
    elif len(nii_fmap.shape) == 4:
        raise NotImplementedError("4D output masks are not supported.")
        nb_vols = nii_fmap.shape[-1]
        nb_ch = None
        nb_slices = nii_fmap.shape[-2]

    # Output fmap is in shape [x y z]
    # Note: The MRD Image class stores data as [ch z y x]
    logging.debug("Output fmap is size %s" % (nii_fmap.shape,))

    if mrdHeader.encoding[0].encodedSpace.matrixSize.z != 1:
        # 3d
        slice_order_nii_to_chrono = {i: nb_slices - i - 1 for i in range(nb_slices)}
    else:
        slice_order_nii_to_chrono = extract_nii_slice_ordering_to_chronological(sidecar, nb_slices)

    # Create a list of MRD Image instances to return
    imagesOut = [None] * nb_slices
    for nii_slice_index in range(nb_slices):
        mrd_slice_index = slice_order_nii_to_chrono[nii_slice_index]
        
        # Select the slice
        slice_dim_in_nifti_coords = dim_of_freq_phase_slice_enc_directions[2]
        if slice_dim_in_nifti_coords == 0:
            tmp = data_fmap[nii_slice_index, :, :]
        elif slice_dim_in_nifti_coords == 1:
            tmp = data_fmap[:, nii_slice_index, :]
        elif slice_dim_in_nifti_coords ==2:
            tmp = data_fmap[:, :, nii_slice_index]
        else:
            raise NotImplementedError("Slice index is not 0, 1 or 2")
        
        # Rotate and flip axes for MRD format
        if slice_dim_in_nifti_coords == 2:
            out = np.flip(np.rot90(tmp, k=-1))[np.newaxis, np.newaxis, ...]
        elif slice_dim_in_nifti_coords == 1:
            out = np.flip(np.rot90(tmp, k=-1))[np.newaxis, np.newaxis, ...]
        elif slice_dim_in_nifti_coords == 0:
            out = np.flip(np.rot90(tmp, k=-1), axis=0)[np.newaxis, np.newaxis, ...]
        else:
            raise NotImplementedError("Slice index is not 0, 1 or 2")

        # Create new MRD instance for the mask
        imagesOut[mrd_slice_index] = ismrmrd.Image.from_array(out, transpose=False)

        # Create a copy of the original fixed header and update the data_type
        # (we changed it to int16 from all other types)
        oldHeader = mag1_head[mrd_slice_index]
        oldHeader.data_type = imagesOut[mrd_slice_index].data_type

        # Increment series number when flag detected (i.e. follow ICE logic for splitting series)
        if mrdhelper.get_meta_value(mag1_meta[mrd_slice_index], 'IceMiniHead') is not None:
            if mrdhelper.extract_minihead_bool_param(base64.b64decode(mag1_meta[mrd_slice_index]['IceMiniHead']).decode('utf-8'), 'BIsSeriesEnd') is True:
                currentSeries += 1

        imagesOut[mrd_slice_index].setHead(oldHeader)

        # Create a copy of the original ISMRMRD Meta attributes and update
        tmpMeta = mag1_meta[mrd_slice_index]
        tmpMeta['DataRole']                       = 'Image'
        tmpMeta['ImageProcessingHistory']         = ['PYTHON', 'FIELDMAP']
        tmpMeta['SequenceDescriptionAdditional']  = 'FIRE'
        tmpMeta['Keep_image_geometry']            = 1

        # Add image orientation directions to MetaAttributes if not already present
        if tmpMeta.get('ImageRowDir') is None:
            tmpMeta['ImageRowDir'] = ["{:.18f}".format(oldHeader.read_dir[0]), "{:.18f}".format(oldHeader.read_dir[1]), "{:.18f}".format(oldHeader.read_dir[2])]

        if tmpMeta.get('ImageColumnDir') is None:
            tmpMeta['ImageColumnDir'] = ["{:.18f}".format(oldHeader.phase_dir[0]), "{:.18f}".format(oldHeader.phase_dir[1]), "{:.18f}".format(oldHeader.phase_dir[2])]

        metaXml = tmpMeta.serialize()
        logging.debug("Image MetaAttributes: %s", xml.dom.minidom.parseString(metaXml).toprettyxml())
        logging.debug("Image data has %d elements", imagesOut[mrd_slice_index].data.size)

        imagesOut[mrd_slice_index].attribute_string = metaXml

        # Debug output
        if "xml" not in dset.list():
            dset.write_xml_header(mrdHeader.toXML())
        dset.append_image("image_%d" % imagesOut[mrd_slice_index].image_series_index, imagesOut[mrd_slice_index])

    dset.close()

    # Send a copy of original (unmodified) images back too
    if mrdhelper.get_json_config_param(config_dict, 'sendOriginal', default=False, type='bool') == True:
        logging.info('Sending a copy of original unmodified images due to sendOriginal set to True')
        # In reverse order so that they'll be in correct order as we insert them to the front of the list
        for image in reversed(imgGroup):
            # Create a copy to not modify the original inputs
            tmpImg = image

            # Change the series_index to have a different series
            tmpImg.image_series_index = 99

            # Ensure Keep_image_geometry is set to not reverse image orientation
            tmpMeta = ismrmrd.Meta.deserialize(tmpImg.attribute_string)
            tmpMeta['Keep_image_geometry'] = 1
            tmpImg.attribute_string = tmpMeta.serialize()

            imagesOut.insert(0, tmpImg)
    
    return imagesOut


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
    dset._file.require_group("dataset")
    return dset
