#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import nibabel as nib
from nibabel.processing import resample_from_to as nib_resample_from_to
import logging
from scipy.ndimage import binary_dilation, binary_erosion, binary_opening, generate_binary_structure, iterate_structure
from joblib import Parallel, delayed
import os


def load_mask(nii_input, nii_mask, path_output):
    # Load the mask
    if nii_mask is None:
        raise ValueError("Input mask is None")
    else:
        # Masks must be 3d
        if len(nii_mask.shape) != 3:
            raise ValueError("Input mask must be 3d")
        # If the mask is of a different shape, resample it.
        elif not np.all(nii_mask.shape == nii_input.shape) or not np.all(nii_mask.affine == nii_input.affine):
            nii_mask = resample_mask(nii_mask, nii_input)
            nib.save(nii_mask, os.path.join(path_output, 'resampled_mask.nii.gz'))
        else:
            logging.info("No need to resample the mask")

    return nii_mask


def resample_mask(nii_mask_from, nii_target, from_slices=None, dilation_kernel='None', dilation_size=None):
    """
    Select the appropriate slices from ``nii_mask_from`` using ``from_slices`` and resample onto ``nii_target``

    Args:
        nii_mask_from (nib.Nifti1Image): Mask to resample from. False or 0 signifies not included.
        nii_target (nib.Nifti1Image): Target image to resample onto.
        from_slices (tuple): Tuple containing the slices to select from nii_mask_from. None selects all the slices.

    Returns:
        nib.Nifti1Image: Mask resampled with nii_target.shape and nii_target.affine.
    """
    mask_from = nii_mask_from.get_fdata()

    if from_slices is None:
        from_slices = tuple(range(mask_from.shape[2]))

    # Initialize a sliced mask and select the slices from from_slices
    sliced_mask = np.full_like(mask_from, fill_value=0)
    sliced_mask[:, :, from_slices] = mask_from[:, :, from_slices]

    # Create nibabel object of sliced mask
    nii_mask = nib.Nifti1Image(sliced_mask.astype(float), nii_mask_from.affine, header=nii_mask_from.header)

    # Resample the sliced mask onto nii_target
    logging.info(nii_mask.shape)
    logging.info(nii_target.shape)
    nii_mask_target = nib_resample_from_to(nii_mask,
                                           nii_target,
                                           order=0,
                                           mode='grid-constant',
                                           cval=0,
                                           out_class=nib.Nifti1Image)

    if dilation_kernel in [None, 'None']:
        return nii_mask_target

    # Dilate the mask to add more pixels in particular directions
    mask_dilated = modify_binary_mask(nii_mask_target.get_fdata(), dilation_kernel, dilation_size, 'dilate')

    nii_full_mask_target = resample_from_to(nii_mask_from, nii_target, order=0, mode='grid-constant', cval=0)

    # Make sure the mask is within the original ROI
    mask_dilated_in_roi = np.logical_and(mask_dilated, nii_full_mask_target.get_fdata())
    nii_mask_dilated = nib.Nifti1Image(mask_dilated_in_roi, nii_mask_target.affine, header=nii_mask_target.header)

    return nii_mask_dilated


def resample_from_to(nii_from_img, nii_to_vox_map, order=2, mode='nearest', cval=0., out_class=nib.Nifti1Image):
    """ Wrapper to nibabel's ``resample_from_to`` function. Resample image `from_img` to mapped voxel space
    `to_vox_map`. The wrapper adds support for 2D input data (adds a singleton) and for 4D time series.
    For more info, refer to nibabel.processing.resample_from_to.

    Args:
        nii_from_img (nibabel.Nifti1Image): Nibabel object with 2D, 3D or 4D array. The 4d case will be treated as a
                                            timeseries.
        nii_to_vox_map (nibabel.Nifti1Image): Nibabel object with
        order (int): Refer to nibabel.processing.resample_from_to
        mode (str): Refer to nibabel.processing.resample_from_to
        cval (scalar): Refer to nibabel.processing.resample_from_to
        out_class: Refer to nibabel.processing.resample_from_to

    Returns:
        nibabel.Nifti1Image: Return a Nibabel object with the resampled data. The 4d case will have an extra dimension
                             for the different time points.

    """

    from_img = nii_from_img.get_fdata()
    if from_img.ndim == 2:
        nii_from_img_3d = nib.Nifti1Image(np.expand_dims(from_img, -1), nii_from_img.affine, header=nii_from_img.header)
        nii_resampled = nib_resample_from_to(nii_from_img_3d, nii_to_vox_map, order=order, mode=mode, cval=cval,
                                             out_class=out_class)

    elif from_img.ndim == 3:
        nii_resampled = nib_resample_from_to(nii_from_img, nii_to_vox_map, order=order, mode=mode, cval=cval,
                                             out_class=out_class)

    elif from_img.ndim == 4:
        nt = from_img.shape[3]
        results = Parallel(-1, backend='loky')(
            delayed(_resample_4d)(it, nii_from_img, nii_to_vox_map, order, mode, cval, out_class)
            for it in range(nt))
        resampled_4d = np.array([results[it] for it in range(nt)]).transpose(1, 2, 3, 0)
        nii_resampled = nib.Nifti1Image(resampled_4d, nii_to_vox_map.affine, header=nii_to_vox_map.header)

    else:
        raise NotImplementedError("Dimensions of input can only be 2D, 3D or 4D")

    return nii_resampled


def _resample_4d(i, nii_from_img, nii_to_vox_map, order, mode, cval, out_class):
    nii_from_img_3d = nib.Nifti1Image(nii_from_img.get_fdata()[..., i], nii_from_img.affine, header=nii_from_img.header)
    resampled_image = nib_resample_from_to(nii_from_img_3d, nii_to_vox_map, order=order, mode=mode,
                                           cval=cval, out_class=out_class).get_fdata()
    return resampled_image


def modify_binary_mask(mask, shape='sphere', size=3, operation='dilate'):
    """
    Dilates or erodes a binary mask according to different shapes and kernel size

    Args:
        mask (numpy.ndarray): 3d array containing the binary mask.
        shape (str): 3d kernel to perform the dilation. Allowed shapes are: 'sphere', 'cross', 'line', 'cube', 'None'.
                     'line' uses 3 line kernels to extend in each directions by "(size - 1) / 2" only if that direction
                     is smaller than (size - 1) / 2
        size (int): Length of a side of the 3d kernel. Must be odd.
        operation (str): Operation to perform. Allowed operations are: 'dilate', 'erode'.

    Returns:
        numpy.ndarray: Dilated/eroded mask.
    """
    mask_operations = {'dilate': binary_dilation, 'erode': binary_erosion}

    if size % 2 == 0 or size < 3:
        raise ValueError("Size must be odd and greater or equal to 3")

    if operation not in mask_operations:
        raise ValueError(f"Operation <{operation}> not supported. Supported operations are: {list(mask_operations.keys())}")

    # Find the middle pixel, will always work since we check size is odd
    mid_pixel = int((size - 1) / 2)

    if shape == 'sphere':
        # Define kernel to perform the dilation
        struct_sphere_size1 = generate_binary_structure(3, 1)
        struct = iterate_structure(struct_sphere_size1, mid_pixel)

        # Dilate
        mask_dilated = mask_operations[operation](mask, structure=struct)

    elif shape == 'None':
        mask_dilated = mask

    else:
        raise ValueError("Use of non supported algorithm for dilating the mask")

    return mask_dilated


def parse_slices(json_data):
    """
    Parse the BIDS sidecar associated with the input nifti file.

    Args:
        fname_nifti (str): Full path to a NIfTI file
    Returns:
        list: 1D list containing tuples of dim3 slices to shim. (dim1, dim2, dim3)
    """

    # Make sure tag SliceTiming exists
    if 'SliceTiming' in json_data:
        slice_timing = json_data['SliceTiming']
    else:
        raise RuntimeError("No tag SliceTiming to automatically parse slice data")

    # If SliceEncodingDirection exists and is negative, SliceTiming is reversed
    if 'SliceEncodingDirection' in json_data:
        if json_data['SliceEncodingDirection'][-1] == '-':
            logging.debug("SliceEncodeDirection is negative, SliceTiming parsed backwards")
            slice_timing.reverse()

    # Return the indexes of the sorted slice_timing
    slice_timing = np.array(slice_timing)
    list_slices = np.argsort(slice_timing)
    slices = []
    # Construct the list of tuples
    while len(list_slices) > 0:
        # Find if the first index has the same timing as other indexes
        # shim_group = tuple(list_slices[list_slices == list_slices[0]])
        shim_group = tuple(np.where(slice_timing == slice_timing[list_slices[0]])[0].astype(np.int32).tolist())
        # Add this as a tuple
        slices.append(shim_group)

        # Since the list_slices is sorted by slice_timing, the only similar values will be at the beginning
        n_to_remove = len(shim_group)
        list_slices = list_slices[n_to_remove:]

    return slices