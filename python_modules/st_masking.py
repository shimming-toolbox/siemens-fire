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
                if dset is None:
                    dset = create_debug_save_file()

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
        if dset is not None:
            dset.close()
        connection.send_close()
                

def process_image(imgGroup, connection, config, mrdHeader, dset):
    # todo: TEMP:
    #config['parameters']['method'] = 'sct_deepseg'
    #config['parameters']['sct_deepseg_create_circular_mask'] = 'True'
    #config['parameters']['sct_deepseg_mask_size'] = '20'
    #config['parameters']['method'] = 'sct_propseg'
    #config['parameters']['sct_propseg_contrast'] = 't2s'
    #config['parameters']['sct_propseg_include_csf'] = 'True'

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
    fname_input_nii = os.path.join(dataFolder, "anat.nii.gz")
    if nii.ndim == 4:
        nii = nib.Nifti1Image(np.mean(nii.get_fdata(), axis=3), nii.affine, header=nii.header)
    nib.save(nii, fname_input_nii)

    method = mrdhelper.get_json_config_param(config_dict, 'method')
    logging.info(f"Creating mask using method: {method}")
    # Create the mask with the threshold method
    if method == 'threshold':
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
    elif method == 'sct_deepseg':
        repetition = None
        for hrd in head:
            if repetition is None:
                repetition = hrd.repetition
            elif repetition != hrd.repetition:
                raise RuntimeError(f"{method} method only supports input images with the same repetition number")

        fname_output_mask_nii = os.path.join(debugFolder, 'mask_sct_deepseg.nii.gz')
        if head[0].repetition == 0:
            env = os.environ.copy()
            env["CUDA_VISIBLE_DEVICES"] = "0"
            env["SCT_USE_GPU"] = "1"
            path_sct_binaries = '/opt/code/spinalcordtoolbox/bin'
            fname_seg_mask_nii = os.path.join(debugFolder, 'mask_sct_deepseg_seg.nii.gz')
            subprocess.run([os.path.join(path_sct_binaries, 'sct_deepseg'), 'spinalcord',
                            '-i', fname_input_nii,
                            '-o', fname_seg_mask_nii],
                           env=env,
                           check=True)
            # Optionally create a circular mask around the centerline
            if mrdhelper.get_json_config_param(config_dict, 'sct_deepseg_create_circular_mask', default=True, type='bool'):
                subprocess.run([os.path.join(path_sct_binaries, 'sct_create_mask'),
                                '-i', fname_input_nii,
                                '-p', f'centerline,{fname_seg_mask_nii}',
                                '-size', mrdhelper.get_json_config_param(config_dict, 'sct_deepseg_mask_size'),
                                '-o', fname_output_mask_nii],
                                check=True)
            else:
                fname_output_mask_nii = fname_seg_mask_nii

            output_mask_nii = nib.load(fname_output_mask_nii)
        else:
            # Only process the first repetition
            output_mask_nii = nib.load(fname_output_mask_nii)
        
        output_mask = output_mask_nii.get_fdata()

    elif method == 'sct_propseg':
        repetition = None
        for hrd in head:
            if repetition is None:
                repetition = hrd.repetition
            elif repetition != hrd.repetition:
                raise RuntimeError(f"{method} method only supports input images with the same repetition number")

        fname_output_mask_nii = os.path.join(debugFolder, 'mask_sct_propseg.nii.gz')
        if head[0].repetition == 0:
            env = os.environ.copy()
            env["CUDA_VISIBLE_DEVICES"] = "0"
            env["SCT_USE_GPU"] = "1"
            path_sct_binaries = '/opt/code/spinalcordtoolbox/bin'
            fname_seg_mask_nii = os.path.join(debugFolder, 'mask_sct_propseg_seg.nii.gz')
            cmd = [os.path.join(path_sct_binaries, 'sct_propseg'),
                            '-c', mrdhelper.get_json_config_param(config_dict, 'sct_propseg_contrast'),
                            '-i', fname_input_nii,
                            '-o', fname_seg_mask_nii]
            include_csf = mrdhelper.get_json_config_param(config_dict, 'sct_propseg_include_csf', default=False, type='bool')
            if include_csf:
                cmd.append('-CSF')
                basename = os.path.basename(fname_input_nii)
                fname_seg_csf_nii = os.path.join(debugFolder, f'{basename.replace(".nii.gz", "_CSF_seg.nii.gz")}')
            subprocess.run(cmd,
                           env=env,
                           check=True)
            
            # Optionally create a circular mask around the centerline
            if include_csf:
                subprocess.run([os.path.join(path_sct_binaries, 'sct_maths'),
                                '-i', fname_seg_mask_nii,
                                '-add', fname_seg_csf_nii,
                                '-o', fname_output_mask_nii],
                                check=True)
            else:
                fname_output_mask_nii = fname_seg_mask_nii

            output_mask_nii = nib.load(fname_output_mask_nii)
        else:
            # Only process the first repetition
            output_mask_nii = nib.load(fname_output_mask_nii)
        
        output_mask = output_mask_nii.get_fdata()

    # Create the mask with the brain segmentation method
    elif method == 'bet':
        repetition = None
        for hrd in head:
            if repetition is None:
                repetition = hrd.repetition
            elif repetition != hrd.repetition:
                raise RuntimeError("bet method only supports input images with the same repetition number")
        
        fname_output_mask_nii = os.path.join(debugFolder, 'mask_bet.nii.gz')
        if head[0].repetition == 0:
            subprocess.run(['/root/shimming-toolbox/python/bin/bet2',
                            fname_input_nii,
                            os.path.join(debugFolder, 'tmp'),
                            '-f', mrdhelper.get_json_config_param(config_dict, 'bet_f', default='0.5'),
                            '-g', mrdhelper.get_json_config_param(config_dict, 'bet_g', default='0'),
                            '-m',
                            '-n'],
                            check=True)
            os.rename(os.path.join(debugFolder, 'tmp_mask.nii.gz'), fname_output_mask_nii)
            output_mask_nii = nib.load(fname_output_mask_nii)
        else:
            # Only process the first repetition
            output_mask_nii = nib.load(fname_output_mask_nii)

        output_mask = output_mask_nii.get_fdata()

    else :
        raise RuntimeError(f"Method {method} is not available. Options are : \
                           'threshold', 'sct_deepseg' 'sct_propseg'and 'bet'.")

    currentSeries = 0

    # Find slice encoding direction
    dim_info = output_mask_nii.header.get_dim_info()

    if None in dim_info:
        if np.allclose(nii.affine, output_mask_nii.affine, rtol=1e-5, atol=1e-5) and nii.shape == output_mask_nii.shape:
            logging.warning("Input image and output mask have the same orientation and shape, assuming slice encoding direction is the same")
            dim_info = nii.header.get_dim_info()

    dim_of_freq_phase_slice_enc_directions = np.argsort(dim_info)

    if len(output_mask.shape) == 3:
        nb_slices = output_mask.shape[dim_of_freq_phase_slice_enc_directions[2]]
    elif len(output_mask.shape) == 4:
        raise NotImplementedError("4D output masks are not supported.")
        nb_vols = output_mask.shape[-1]
        nb_ch = None
        nb_slices = output_mask.shape[-2]

    # Output mask is in shape [x y z]
    # Note: The MRD Image class stores data as [ch z y x]
    logging.debug("Output mask is size %s" % (output_mask.shape,))

    if mrdHeader.encoding[0].encodedSpace.matrixSize.z != 1:
        # 3d
        slice_order_nii_to_chrono = {i: i for i in range(nb_slices)}
    else:
        slice_order_nii_to_chrono = extract_nii_slice_ordering_to_chronological(sidecar, nb_slices)

    # Create a list of MRD Image instances to return
    imagesOut = [None] * nb_slices
    for nii_slice_index in range(nb_slices):
        mrd_slice_index = slice_order_nii_to_chrono[nii_slice_index]
        
        # Select the slice
        slice_dim_in_nifti_coords = dim_of_freq_phase_slice_enc_directions[2]
        if slice_dim_in_nifti_coords == 0:
            tmp = output_mask[nii_slice_index, :, :]
        elif slice_dim_in_nifti_coords == 1:
            tmp = output_mask[:, nii_slice_index, :]
        elif slice_dim_in_nifti_coords ==2:
            tmp = output_mask[:, :, nii_slice_index]
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
        
        # Scale to int16 max
        if out.max() > 1:
            raise RuntimeError("Output mask has values greater than 1. Scaling to max int16 will result in an overflow.")
        if out.min() < -1:
            raise RuntimeError("Output mask has values less than 1. Scaling to max int16 will result in an overflow.")
        out *= 32767 
        out = out.astype(np.int16)
        data = np.array(imgGroup[mrd_slice_index].data).astype(np.float64)

        # Scale to int16
        bits_stored = 12
        if (mrdhelper.get_userParameterLong_value(mrdHeader, "BitsStored") is not None):
            bits_stored = mrdhelper.get_userParameterLong_value(mrdHeader, "BitsStored")
        max_val = (2 ** bits_stored) - 1
        data *= max_val / data.max()
        data = np.around(data)
        data = data.astype(np.int16)
        out[out == 0] = data[out == 0]

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
        if "xml" not in dset.list():
            dset.write_xml_header(mrdHeader.toXML())
        dset.append_image("image_%d" % imagesOut[mrd_slice_index].image_series_index, imagesOut[mrd_slice_index])
        

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
