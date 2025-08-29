# Adapted from Kelvin Chow's python-ismrmrd-server repository
# Source: https://github.com/kspaceKelvin/python-ismrmrd-server

from datetime import datetime
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
from mrd2nii.mrd2nii_main import mrd2nii_volume


# Folder for debug output files
debugFolder = "/tmp/share/debug"
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

    # Continuously parse incoming data parsed from MRD messages
    currentSeries = 0
    imgGroup = []
    waveformGroup = []
    dset = None
    try:
        for item in connection:
            # ----------------------------------------------------------
            # Image data messages
            # ----------------------------------------------------------
            if isinstance(item, ismrmrd.Image):
                # When this criteria is met, run process_group() on the accumulated
                # data, which returns images that are sent back to the client.
                # e.g. when the series number changes:
                if item.image_series_index != currentSeries:
                    logging.info("Processing a group of images because series index changed to %d", item.image_series_index)
                    currentSeries = item.image_series_index
                    if dset is None:
                        dset = create_debug_save_file()
                    image = process_image(imgGroup, connection, config, mrdHeader, dset)
                    connection.send_image(image)
                    imgGroup = []
                # Only process magnitude images -- send phase images back without modification (fallback for images with unknown type)
                if (item.image_type is ismrmrd.IMTYPE_MAGNITUDE) or (item.image_type == 0):
                    imgGroup.append(item)
                else:
                    tmpMeta = ismrmrd.Meta.deserialize(item.attribute_string)
                    tmpMeta['Keep_image_geometry']    = 1
                    item.attribute_string = tmpMeta.serialize()

                    connection.send_image(item)
                    continue

            # ----------------------------------------------------------
            # Waveform data messages
            # ----------------------------------------------------------
            elif isinstance(item, ismrmrd.Waveform):
                waveformGroup.append(item)
            elif item is None:
                break
            else:
                logging.error("Unsupported data type %s", type(item).__name__)

        # Extract raw ECG waveform data. Basic sorting to make sure that data 
        # is time-ordered, but no additional checking for missing data.
        # ecgData has shape (5 x timepoints)
        if len(waveformGroup) > 0:
            waveformGroup.sort(key = lambda item: item.time_stamp)
            ecgData = [item.data for item in waveformGroup if item.waveform_id == 0]
            if len(ecgData) > 0:
                ecgData = np.concatenate(ecgData,1)

        # Process any remaining groups of image data.  This can 
        # happen if the trigger condition for these groups are not met.
        # This is also a fallback for handling image data, as the last
        # image in a series is typically not separately flagged.
        if len(imgGroup) > 0:
            if dset is None:
                dset = create_debug_save_file()
            logging.info("Processing a group of images (untriggered)")
            image = process_image(imgGroup, connection, config, mrdHeader, dset)
            connection.send_image(image)
            imgGroup = []

    except Exception as e:
        logging.error(traceback.format_exc())
        connection.send_logging(constants.MRD_LOGGING_ERROR, traceback.format_exc())

    finally:
        connection.send_close()
        if dset is not None:
            dset.close()

def process_image(imgGroup, connection, config, mrdHeader, dset):
    config_filename = f"{config}.json"
    config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), config_filename)
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"{config_path} does not exist.")

    with open(config_path, "r") as f:
        config_dict = json.load(f)

    if len(imgGroup) == 0:
        return []

    logging.info('-----------------------------------------------')
    logging.info(f'process_image called with {len(imgGroup)} images')
    logging.info('-----------------------------------------------')

    # Create folder, if necessary
    if not os.path.exists(debugFolder):
        os.makedirs(debugFolder)
        logging.debug("Created folder " + debugFolder + " for debug output files")

    logging.debug("Processing data with %d images of type %s", len(imgGroup), ismrmrd.get_dtype_from_data_type(imgGroup[0].data_type))

    head = [img.getHead()                                  for img in imgGroup]
    meta = [ismrmrd.Meta.deserialize(img.attribute_string) for img in imgGroup]
    
    # Display MetaAttributes for first image
    logging.debug("MetaAttributes[0]: %s", ismrmrd.Meta.serialize(meta[0]))

    # Optional serialization of ICE MiniHeader
    if 'IceMiniHead' in meta[0]:
        logging.debug("IceMiniHead[0]: %s", base64.b64decode(meta[0]['IceMiniHead']).decode('utf-8'))

    # Copy the MRD file to a new location for conversion
    mrd2nii_folder = debugFolder + "/mrd2nii_conversion"
    if not os.path.exists(mrd2nii_folder):
        os.makedirs(mrd2nii_folder)

    # Convert the MRD images to NIfTI format
    nii, sidecar = mrd2nii_volume(mrdHeader, imgGroup)
    fname_input_nii = os.path.join(mrd2nii_folder, "img.nii.gz")
    nib.save(nii, fname_input_nii)

    # Create the mask with the threshold method
    if mrdhelper.get_json_config_param(config_dict, 'method') == 'threshold':
        fname_output_mask_nii = os.path.join(debugFolder, 'mask_threshold.nii.gz')
        subprocess.run(['st_mask', 'threshold',
                        '-i', fname_input_nii,
                        '--scaled-thr',
                        '--thr', mrdhelper.get_json_config_param(config_dict, 'threshold_thr'),
                        '-o', fname_output_mask_nii],
                        check=True)
        output_mask_nii = nib.load(fname_output_mask_nii)
        output_mask = output_mask_nii.get_fdata()

    # Create the mask with the SC segmentation method
    elif mrdhelper.get_json_config_param(config_dict, 'method') == 'sct_deepseg':
        repetition = None
        for hrd in head:
            if repetition is None:
                repetition = hrd.repetition
            elif repetition != hrd.repetition:
                raise RuntimeError("sct_deepseg method only supports input images with the same repetition number")

        if head[0].repetition == 0:
            fname_seg_mask_nii = os.path.join(debugFolder, 'mask_sct_deepseg_seg.nii.gz')
            subprocess.run(['sct_deepseg', 'spinalcord',
                            '-i', fname_input_nii,
                            '-o', fname_seg_mask_nii],
                            check=True)
            fname_output_mask_nii = os.path.join(debugFolder, 'mask_sct_deepseg.nii.gz')
            subprocess.run(['sct_create_mask',
                            '-i', fname_input_nii,
                            '-p', f'centerline,{fname_seg_mask_nii}',
                            '-size', mrdhelper.get_json_config_param(config_dict, 'sct_deepseg_mask_size'),
                            '-o', fname_output_mask_nii],
                            check=True)
            output_mask_nii = nib.load(fname_output_mask_nii)
        else:
            # Only process the first repetition
            output_mask_nii = nib.load(fname_input_nii)
        
        output_mask = output_mask_nii.get_fdata()

    # Create the mask with the brain segmentation method
    elif mrdhelper.get_json_config_param(config_dict, 'method') == 'bet':
        raise NotImplementedError(f"Method {mrdhelper.get_json_config_param(config_dict, 'method')} is not implemented yet")
        # TODO: Implement BET method

    else :
        raise RuntimeError(f"Method {mrdhelper.get_json_config_param(config_dict, 'method')} is not available. Options are : \
                           'threshold', 'sct_deepseg' and 'bet'.")

    currentSeries = 0

    if len(output_mask.shape) == 3:
        nb_z = output_mask.shape[-1]  
    elif len(output_mask.shape) == 4:
        raise NotImplementedError("4D output masks are not supported.")
        nb_vols = output_mask.shape[-1]
        nb_ch = None
        nb_z = output_mask.shape[-2]

    # Output mask is in shape [x y z]
    # Note: The MRD Image class stores data as [ch z y x]
    logging.debug("Output mask is size %s" % (output_mask.shape,))

    slice_order_nii_to_chrono = extract_nii_slice_ordering_to_chronological(sidecar, nb_z)

    # Create a list of MRD Image instances to return
    imagesOut = [None] * nb_z
    for nii_slice_index in range(nb_z):
        mrd_slice_index = slice_order_nii_to_chrono[nii_slice_index]
        # Rotate each mask by 90 degrees clock-wise in the xy plane
        out = np.flip(np.rot90(output_mask[:, :, nii_slice_index], k=-1))[np.newaxis, np.newaxis, ...]
        # out = output_mask[:, :, mrd_slice_index][np.newaxis, np.newaxis, ...]

        # Create new MRD instance for the mask
        imagesOut[mrd_slice_index] = ismrmrd.Image.from_array(out, transpose=False)

        # Create a copy of the original fixed header and update the data_type
        # (we changed it to int16 from all other types)
        oldHeader = head[mrd_slice_index]
        oldHeader.data_type = imagesOut[mrd_slice_index].data_type

        # Increment series number when flag detected (i.e. follow ICE logic for splitting series)
        if mrdhelper.get_meta_value(meta[mrd_slice_index], 'IceMiniHead') is not None:
            if mrdhelper.extract_minihead_bool_param(base64.b64decode(meta[mrd_slice_index]['IceMiniHead']).decode('utf-8'), 'BIsSeriesEnd') is True:
                currentSeries += 1

        imagesOut[mrd_slice_index].setHead(oldHeader)

        # Create a copy of the original ISMRMRD Meta attributes and update
        tmpMeta = meta[mrd_slice_index]
        tmpMeta['DataRole']                       = 'Image'
        tmpMeta['ImageProcessingHistory']         = ['PYTHON', 'MASK']
        tmpMeta['SequenceDescriptionAdditional']  = 'FIRE'
        tmpMeta['Keep_image_geometry']            = 1

        if mrdhelper.get_json_config_param(config_dict, 'method') == 'threshold':
            tmpMeta['ImageProcessingHistory'].append('THRESHOLD')
        if mrdhelper.get_json_config_param(config_dict, 'method') == 'sct_deepseg':
            tmpMeta['ImageProcessingHistory'].append('SC SEGMENTATION')
        if mrdhelper.get_json_config_param(config_dict, 'method') == 'bet':
            tmpMeta['ImageProcessingHistory'].append('BRAIN SEGMENTATION')

        # Add image orientation directions to MetaAttributes if not already present
        if tmpMeta.get('ImageRowDir') is None:
            tmpMeta['ImageRowDir'] = ["{:.18f}".format(oldHeader.read_dir[0]), "{:.18f}".format(oldHeader.read_dir[1]), "{:.18f}".format(oldHeader.read_dir[2])]

        if tmpMeta.get('ImageColumnDir') is None:
            tmpMeta['ImageColumnDir'] = ["{:.18f}".format(oldHeader.phase_dir[0]), "{:.18f}".format(oldHeader.phase_dir[1]), "{:.18f}".format(oldHeader.phase_dir[2])]

        metaXml = tmpMeta.serialize()
        logging.debug("Image MetaAttributes: %s", xml.dom.minidom.parseString(metaXml).toprettyxml())
        logging.debug("Image data has %d elements", imagesOut[mrd_slice_index].data.size)

        imagesOut[mrd_slice_index].attribute_string = metaXml

        # Copy to main data directory
        shutil.copyfile(fname_output_mask_nii, os.path.join(dataFolder, "mask.nii.gz"))

        # Debug output
        dset.append_image("image_%d" % imagesOut[mrd_slice_index].image_series_index, imagesOut[mrd_slice_index])
        if "xml" not in dset.list():
            dset.write_xml_header(mrdHeader.toXML())
        

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
    return dset
