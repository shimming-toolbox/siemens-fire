import ismrmrd
import os
import logging
import traceback
import json
import numpy as np
import numpy.fft as fft
import xml.dom.minidom
import base64
import ctypes
import mrdhelper # Custom module for MRD helper functions found in the python-ismrmrd-server repository
import constants # Custom module for constants found in the python-ismrmrd-server repository
import subprocess
import nibabel as nib
import shutil
from time import perf_counter
import glob

# Folder for debug output files
debugFolder = "/tmp/share/debug"

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
    acqGroup = []
    imgGroup = []
    waveformGroup = []
    try:
        for item in connection:
            # ----------------------------------------------------------
            # Raw k-space data messages
            # ----------------------------------------------------------
            if isinstance(item, ismrmrd.Acquisition):
                # Accumulate all imaging readouts in a group
                if (not item.is_flag_set(ismrmrd.ACQ_IS_NOISE_MEASUREMENT) and
                    not item.is_flag_set(ismrmrd.ACQ_IS_PARALLEL_CALIBRATION) and
                    not item.is_flag_set(ismrmrd.ACQ_IS_PHASECORR_DATA) and
                    not item.is_flag_set(ismrmrd.ACQ_IS_NAVIGATION_DATA)):
                    acqGroup.append(item)

                # When this criteria is met, run process_raw() on the accumulated
                # data, which returns images that are sent back to the client.
                if item.is_flag_set(ismrmrd.ACQ_LAST_IN_SLICE):
                    logging.info("Processing a group of k-space data")
                    image = process_raw(acqGroup, connection, config, mrdHeader)
                    connection.send_image(image)
                    acqGroup = []

            # ----------------------------------------------------------
            # Image data messages
            # ----------------------------------------------------------
            elif isinstance(item, ismrmrd.Image):
                # When this criteria is met, run process_group() on the accumulated
                # data, which returns images that are sent back to the client.
                # e.g. when the series number changes:
                if item.image_series_index != currentSeries:
                    logging.info("Processing a group of images because series index changed to %d", item.image_series_index)
                    currentSeries = item.image_series_index
                    image = process_image(imgGroup, connection, config, mrdHeader)
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

        # Process any remaining groups of raw or image data.  This can 
        # happen if the trigger condition for these groups are not met.
        # This is also a fallback for handling image data, as the last
        # image in a series is typically not separately flagged.
        if len(acqGroup) > 0:
            logging.info("Processing a group of k-space data (untriggered)")
            image = process_raw(acqGroup, connection, config, mrdHeader)
            connection.send_image(image)
            acqGroup = []

        if len(imgGroup) > 0:
            logging.info("Processing a group of images (untriggered)")
            image = process_image(imgGroup, connection, config, mrdHeader)
            connection.send_image(image)
            imgGroup = []

    except Exception as e:
        logging.error(traceback.format_exc())
        connection.send_logging(constants.MRD_LOGGING_ERROR, traceback.format_exc())

    finally:
        connection.send_close()


def process_raw(acqGroup, connection, config, mrdHeader):
    if len(acqGroup) == 0:
        return []
    
    logging.info(f'-----------------------------------------------')
    logging.info(f'     process_raw called with {len(acqGroup)} readouts')
    logging.info(f'-----------------------------------------------')

    # Start timer
    tic = perf_counter()

    # Create folder, if necessary
    if not os.path.exists(debugFolder):
        os.makedirs(debugFolder)
        logging.debug("Created folder " + debugFolder + " for debug output files")

    # Format data into single [cha PE RO phs] array
    lin = [acquisition.idx.kspace_encode_step_1 for acquisition in acqGroup]
    phs = [acquisition.idx.phase                for acquisition in acqGroup]

    # Use the zero-padded matrix size
    data = np.zeros((acqGroup[0].data.shape[0], 
                     mrdHeader.encoding[0].encodedSpace.matrixSize.y, 
                     mrdHeader.encoding[0].encodedSpace.matrixSize.x, 
                     max(phs)+1), 
                    acqGroup[0].data.dtype)

    rawHead = [None]*(max(phs)+1)

    for acq, lin, phs in zip(acqGroup, lin, phs):
        if (lin < data.shape[1]) and (phs < data.shape[3]):
            # TODO: Account for asymmetric echo in a better way
            data[:,lin,-acq.data.shape[1]:,phs] = acq.data

            # center line of k-space is encoded in user[5]
            if (rawHead[phs] is None) or (np.abs(acq.getHead().idx.kspace_encode_step_1 - acq.getHead().idx.user[5]) < np.abs(rawHead[phs].idx.kspace_encode_step_1 - rawHead[phs].idx.user[5])):
                rawHead[phs] = acq.getHead()

    # Flip matrix in RO/PE to be consistent with ICE
    data = np.flip(data, (1, 2))

    logging.debug("Raw data is size %s" % (data.shape,))
    np.save(debugFolder + "/" + "raw.npy", data)

    # Fourier Transform
    data = fft.fftshift( data, axes=(1, 2))
    data = fft.ifft2(    data, axes=(1, 2))
    data = fft.ifftshift(data, axes=(1, 2))
    data *= np.prod(data.shape) # FFT scaling for consistency with ICE

    # Sum of squares coil combination
    # Data will be [PE RO phs]
    data = np.abs(data)
    data = np.square(data)
    data = np.sum(data, axis=0)
    data = np.sqrt(data)

    logging.debug("Image data is size %s" % (data.shape,))
    np.save(debugFolder + "/" + "img.npy", data)

    # Remove readout oversampling
    if mrdHeader.encoding[0].reconSpace.matrixSize.x != 0:
        offset = int((data.shape[1] - mrdHeader.encoding[0].reconSpace.matrixSize.x)/2)
        data = data[:,offset:offset+mrdHeader.encoding[0].reconSpace.matrixSize.x]

    # Remove phase oversampling
    if mrdHeader.encoding[0].reconSpace.matrixSize.y != 0:
        offset = int((data.shape[0] - mrdHeader.encoding[0].reconSpace.matrixSize.y)/2)
        data = data[offset:offset+mrdHeader.encoding[0].reconSpace.matrixSize.y,:]

    logging.debug("Image without oversampling is size %s" % (data.shape,))
    np.save(debugFolder + "/" + "imgCrop.npy", data)

    # Measure processing time
    toc = perf_counter()
    strProcessTime = "Total processing time: %.2f ms" % ((toc-tic)*1000.0)
    logging.info(strProcessTime)

    # Send this as a text message back to the client
    connection.send_logging(constants.MRD_LOGGING_INFO, strProcessTime)

    # Format as ISMRMRD image data
    imagesOut = []
    for phs in range(data.shape[2]):
        # Create new MRD instance for the processed image
        # data has shape [PE RO phs], i.e. [y x].
        # from_array() should be called with 'transpose=False' to avoid warnings, and when called
        # with this option, can take input as: [cha z y x], [z y x], or [y x]
        tmpImg = ismrmrd.Image.from_array(data[...,phs], transpose=False)

        # Set the header information
        tmpImg.setHead(mrdhelper.update_img_header_from_raw(tmpImg.getHead(), rawHead[phs]))
        tmpImg.field_of_view = (ctypes.c_float(mrdHeader.encoding[0].reconSpace.fieldOfView_mm.x), 
                                ctypes.c_float(mrdHeader.encoding[0].reconSpace.fieldOfView_mm.y), 
                                ctypes.c_float(mrdHeader.encoding[0].reconSpace.fieldOfView_mm.z))
        tmpImg.image_index = phs

        # Set ISMRMRD Meta Attributes
        tmpMeta = ismrmrd.Meta()
        tmpMeta['DataRole']               = 'Image'
        tmpMeta['ImageProcessingHistory'] = ['FIRE', 'PYTHON']
        tmpMeta['Keep_image_geometry']    = 1

        xml = tmpMeta.serialize()
        logging.debug("Image MetaAttributes: %s", xml)
        tmpImg.attribute_string = xml
        imagesOut.append(tmpImg)

    # Call process_image() to invert image contrast
    imagesOut = process_image(imagesOut, connection, config, mrdHeader)

    return imagesOut


def process_image(imgGroup, connection, config, mrdHeader):
    config_filename = f"{config}.json"
    config_path = None
    search_dir = "/opt/code/python-ismrmrd-server/"
    for file in glob.glob(os.path.join(search_dir, "*.json")):
        if os.path.basename(file) == config_filename:
            config_path = file
            break

    if config_path is None:
        raise FileNotFoundError(f"Could not find {config_filename} in {search_dir}.")

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
    fname_input_mrd = connection.mrdFilePath
    fname_copied_mrd = os.path.join(mrd2nii_folder, "imgOrig.h5")
    shutil.copy(fname_input_mrd, fname_copied_mrd)

    # Convert the MRD images to NIfTI format
    subprocess.run(['mrd2nii', '-i', mrd2nii_folder, '-o', mrd2nii_folder], check=True)

    # Rename the output NIfTI file
    fname_input_nii = os.path.join(mrd2nii_folder, "imgOrig.nii.gz")
    nii_files = [f for f in os.listdir(mrd2nii_folder) if f.endswith('.nii.gz')]
    if not nii_files:
        raise FileNotFoundError("No .nii.gz file found in {}".format(mrd2nii_folder))
    first_nii = os.path.join(mrd2nii_folder, nii_files[0])
    os.rename(first_nii, fname_input_nii)

    fname_output_mask_nii = debugFolder + '/mask_threshold.nii.gz'

    # Create the mask with the threshold method
    if mrdhelper.get_json_config_param(config_dict, 'method') == 'threshold':
        fname_output_mask_nii = debugFolder + '/mask_threshold.nii.gz'
        subprocess.run(['st_mask', 'threshold',
                        '-i', fname_input_nii,
                        '--thr', mrdhelper.get_json_config_param(config_dict, 'value'),
                        '-o', fname_output_mask_nii],
                        check=True)
        output_mask_nii = nib.load(fname_output_mask_nii)
        output_mask = output_mask_nii.get_fdata()

    # Create the mask with the SC segmentation method
    elif mrdhelper.get_json_config_param(config_dict, 'method') == 'sct_deepseg':
        raise NotImplementedError(f"Method {mrdhelper.get_json_config_param(config_dict, 'method')} is not implemented yet")
        # TODO: Implement SCT DeepSeg method

    # Create the mask with the brain segmentation method
    elif mrdhelper.get_json_config_param(config_dict, 'method') == 'bet':
        raise NotImplementedError(f"Method {mrdhelper.get_json_config_param(config_dict, 'method')} is not implemented yet")
        # TODO: Implement BET method

    else :
        raise RuntimeError(f"Method {mrdhelper.get_json_config_param(config_dict, 'method')} is not available. Options are : \
                           'threshold', 'sct_deepseg' and 'bet'.")

    currentSeries = 0

    # Output mask is in shape [x y z ch]
    # Note: The MRD Image class stores data as [ch z y x]
    # Extract mask data into a 5D array of size [img=ch*z 1 1 y x]
    nb_ch = output_mask.shape[-1]
    nb_z = output_mask.shape[-2]

    logging.debug("Output mask is size %s" % (output_mask.shape,))
    out = np.stack([output_mask[:, :, i_z, i_ch] for i_ch in range(nb_ch) for i_z in range(nb_z)], axis=0) # shape [img x y]
    out = out[:, np.newaxis, np.newaxis, :, :] # shape [img 1 1 x y]
    # Rotate each mask by 90 degrees clock-wise in the xy plane
    out = np.rot90(out, k=3, axes=(-1, -2)) # shape [img 1 1 y x]

    # Create a list of MRD Image instances to return
    nb_img = nb_ch * nb_z
    imagesOut = [None] * nb_img
    for iImg in range(nb_img):
        # Create new MRD instance for the mask
        imagesOut[iImg] = ismrmrd.Image.from_array(out[iImg, ...], transpose=False)

        # Create a copy of the original fixed header and update the data_type
        # (we changed it to int16 from all other types)
        oldHeader = head[iImg]
        oldHeader.data_type = imagesOut[iImg].data_type

        # Increment series number when flag detected (i.e. follow ICE logic for splitting series)
        if mrdhelper.get_meta_value(meta[iImg], 'IceMiniHead') is not None:
            if mrdhelper.extract_minihead_bool_param(base64.b64decode(meta[iImg]['IceMiniHead']).decode('utf-8'), 'BIsSeriesEnd') is True:
                currentSeries += 1

        imagesOut[iImg].setHead(oldHeader)

        # Create a copy of the original ISMRMRD Meta attributes and update
        tmpMeta = meta[iImg]
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
        logging.debug("Image data has %d elements", imagesOut[iImg].data.size)

        imagesOut[iImg].attribute_string = metaXml

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
