import os
import shutil
import subprocess
import json
from ismrmrd import Dataset
import h5py
import nibabel as nib
import numpy as np
import glob

from python_modules import __PATH_TESTING_DATA__, __PATH_REPO__, __TMP_SHARE_DEBUG__, __TMP_SHARE_SAVEDDATA__
from . import DEBUG


def test_st_masking_sct_deepseg_20mm_circular(tmpdir):
    fname_logfile = os.path.join(__PATH_TESTING_DATA__, "log", "test_st_masking_sct_deepseg_20mm_circular.log")
    if os.path.exists(fname_logfile):
        os.remove(fname_logfile)

    config = {
        "version": "1.0",
        "parameters": {
            "method": "sct_deepseg",
            "threshold_thr": "0.1",
            "sct_deepseg_create_circular_mask": "True",
            "sct_deepseg_mask_size": "20",
            "sct_propseg_contrast": "t2s",
            "sct_propseg_include_csf": "True",
            "sendOriginal": "False",
            "bet_f": "0.5",
            "bet_g": "0.0"
        }
    }

    fname_input = os.path.join(__PATH_TESTING_DATA__, "localizer_for_segmentation_2025-10-23-155042_92_invivo.h5")
    dset = Dataset(fname_input)
    dsetConfigAdditional = dset._dataset.require_dataset('configAdditional',shape=(1,), dtype=h5py.special_dtype(vlen=bytes))
    dsetConfigAdditional[0] = bytes(json.dumps(config), 'utf-8')
    dset.close()

    ret = subprocess.run(
        [os.path.join("/root", "shimming-toolbox", "python", "bin", "python"),
         os.path.join("/opt", "code", "python-ismrmrd-server", "client.py"),
         fname_input,
         "-p", "9020",
         "-c", "st_masking",
         "-o", tmpdir / "output.mrd",
         "-l", fname_logfile,
         ], 
         check=True,
    )

    assert os.path.exists(tmpdir / "output.mrd")
    assert ret.returncode == 0

    fname_anat = os.path.join(__TMP_SHARE_SAVEDDATA__, "anat.nii.gz")
    fname_mask = os.path.join(__TMP_SHARE_SAVEDDATA__, "mask.nii.gz")
    assert os.path.exists(fname_mask)
    assert os.path.exists(fname_anat)
    assert nib.load(fname_mask).shape == nib.load(fname_anat).shape
    assert np.allclose(nib.load(fname_mask).affine, nib.load(fname_anat).affine)
    assert nib.load(fname_anat).get_fdata().max() > 0
    assert nib.load(fname_mask).get_fdata()[121,126,46] == 1
    assert nib.load(fname_mask).get_fdata()[121,171,38] == 0
    
    if DEBUG:
        path_debug = os.path.join(__PATH_REPO__, "data", "debug")
        if os.path.exists(path_debug):
            shutil.rmtree(path_debug)
        shutil.copytree(__TMP_SHARE_SAVEDDATA__, path_debug, dirs_exist_ok=True)


def test_st_masking_sct_propseg_csf(tmpdir):
    fname_logfile = os.path.join(__PATH_TESTING_DATA__, "log", "test_st_masking_sct_propseg_csf.log")
    if os.path.exists(fname_logfile):
        os.remove(fname_logfile)

    config = {
        "version": "1.0",
        "parameters": {
            "method": "sct_propseg",
            "threshold_thr": "0.1",
            "sct_deepseg_create_circular_mask": "True",
            "sct_deepseg_mask_size": "20",
            "sct_propseg_contrast": "t2s",
            "sct_propseg_include_csf": "True",
            "sendOriginal": "False",
            "bet_f": "0.5",
            "bet_g": "0.0"
        }
    }

    fname_input = os.path.join(__PATH_TESTING_DATA__, "localizer_for_segmentation_2025-10-23-155042_92_invivo.h5")
    dset = Dataset(fname_input)
    dsetConfigAdditional = dset._dataset.require_dataset('configAdditional',shape=(1,), dtype=h5py.special_dtype(vlen=bytes))
    dsetConfigAdditional[0] = bytes(json.dumps(config), 'utf-8')
    dset.close()

    ret = subprocess.run(
        [os.path.join("/root", "shimming-toolbox", "python", "bin", "python"),
         os.path.join("/opt", "code", "python-ismrmrd-server", "client.py"),
         fname_input,
         "-p", "9020",
         "-c", "st_masking",
         "-o", tmpdir / "output.mrd",
         "-l", fname_logfile,
         ], 
         check=True,
    )

    assert os.path.exists(tmpdir / "output.mrd")
    assert ret.returncode == 0

    fname_anat = os.path.join(__TMP_SHARE_SAVEDDATA__, "anat.nii.gz")
    fname_mask = os.path.join(__TMP_SHARE_SAVEDDATA__, "mask.nii.gz")
    assert os.path.exists(fname_mask)
    assert os.path.exists(fname_anat)
    assert nib.load(fname_mask).shape == nib.load(fname_anat).shape
    assert np.allclose(nib.load(fname_mask).affine, nib.load(fname_anat).affine)
    assert nib.load(fname_anat).get_fdata().max() > 0
    assert nib.load(fname_mask).get_fdata()[121,126,46] == 1
    assert nib.load(fname_mask).get_fdata()[121,171,38] == 0
    
    if DEBUG:
        path_debug = os.path.join(__PATH_REPO__, "data", "debug")
        if os.path.exists(path_debug):
            shutil.rmtree(path_debug)
        shutil.copytree(__TMP_SHARE_SAVEDDATA__, path_debug, dirs_exist_ok=True)


def test_st_masking_threshold(tmpdir):
    fname_logfile = os.path.join(__PATH_TESTING_DATA__, "log", "test_st_masking_threshold.log")
    if os.path.exists(fname_logfile):
        os.remove(fname_logfile)

    config = {
        "version": "1.0",
        "parameters": {
            "method": "threshold",
            "threshold_thr": "0.1",
            "sct_deepseg_create_circular_mask": "True",
            "sct_deepseg_mask_size": "20",
            "sct_propseg_contrast": "t2s",
            "sct_propseg_include_csf": "True",
            "sendOriginal": "False",
            "bet_f": "0.5",
            "bet_g": "0.0"
        }
    }

    fname_input = os.path.join(__PATH_TESTING_DATA__, "localizer_for_segmentation_2025-10-23-155042_92_invivo.h5")
    dset = Dataset(fname_input)
    dsetConfigAdditional = dset._dataset.require_dataset('configAdditional',shape=(1,), dtype=h5py.special_dtype(vlen=bytes))
    dsetConfigAdditional[0] = bytes(json.dumps(config), 'utf-8')
    dset.close()

    ret = subprocess.run(
        [os.path.join("/root", "shimming-toolbox", "python", "bin", "python"),
         os.path.join("/opt", "code", "python-ismrmrd-server", "client.py"),
         fname_input,
         "-p", "9020",
         "-c", "st_masking",
         "-o", tmpdir / "output.mrd",
         "-l", fname_logfile,
         ], 
         check=True,
    )

    assert os.path.exists(tmpdir / "output.mrd")
    assert ret.returncode == 0

    fname_anat = os.path.join(__TMP_SHARE_SAVEDDATA__, "anat.nii.gz")
    fname_mask = os.path.join(__TMP_SHARE_SAVEDDATA__, "mask.nii.gz")
    assert os.path.exists(fname_mask)
    assert os.path.exists(fname_anat)
    assert nib.load(fname_mask).shape == nib.load(fname_anat).shape
    assert np.allclose(nib.load(fname_mask).affine, nib.load(fname_anat).affine)
    assert nib.load(fname_anat).get_fdata().max() > 0
    assert nib.load(fname_mask).get_fdata()[121,126,46] == 1
    assert nib.load(fname_mask).get_fdata()[121,171,38] == 1

    path_debug = os.path.join(__PATH_REPO__, "data", "debug")
    if os.path.exists(path_debug):
        shutil.rmtree(path_debug)
    shutil.copytree(__TMP_SHARE_SAVEDDATA__, path_debug, dirs_exist_ok=True)

    # We run mrd2nii on the debug folder to convert the output.mrd to nifti for inspection
    subprocess.run(
        ["mrd2nii",
         "-i", __TMP_SHARE_DEBUG__,
         "-o", path_debug],
        check=True
    )

    data_output = nib.load(glob.glob(os.path.join(path_debug, "*.nii.gz"))[0]).get_fdata()
    data_mask_output = np.zeros_like(data_output)
    data_mask_output[data_output > 30000] = 1
    assert np.allclose(data_mask_output, nib.load(fname_mask).get_fdata())


def test_st_masking_bet(tmpdir):
    fname_logfile = os.path.join(__PATH_TESTING_DATA__, "log", "test_st_masking_bet.log")
    if os.path.exists(fname_logfile):
        os.remove(fname_logfile)

    config = {
        "version": "1.0",
        "parameters": {
            "method": "bet",
            "threshold_thr": "0.1",
            "sct_deepseg_create_circular_mask": "True",
            "sct_deepseg_mask_size": "20",
            "sct_propseg_contrast": "t2s",
            "sct_propseg_include_csf": "True",
            "sendOriginal": "False",
            "bet_f": "0.5",
            "bet_g": "0.0"
        }
    }

    fname_input = os.path.join(__PATH_TESTING_DATA__, "localizer_for_segmentation_2025-10-23-155042_92_invivo.h5")
    dset = Dataset(fname_input)
    dsetConfigAdditional = dset._dataset.require_dataset('configAdditional',shape=(1,), dtype=h5py.special_dtype(vlen=bytes))
    dsetConfigAdditional[0] = bytes(json.dumps(config), 'utf-8')
    dset.close()

    ret = subprocess.run(
        [os.path.join("/root", "shimming-toolbox", "python", "bin", "python"),
         os.path.join("/opt", "code", "python-ismrmrd-server", "client.py"),
         fname_input,
         "-p", "9020",
         "-c", "st_masking",
         "-o", tmpdir / "output.mrd",
         "-l", fname_logfile,
         ], 
         check=True,
    )

    assert os.path.exists(tmpdir / "output.mrd")
    assert ret.returncode == 0

    fname_anat = os.path.join(__TMP_SHARE_SAVEDDATA__, "anat.nii.gz")
    fname_mask = os.path.join(__TMP_SHARE_SAVEDDATA__, "mask.nii.gz")
    assert os.path.exists(fname_mask)
    assert os.path.exists(fname_anat)
    assert nib.load(fname_mask).shape == nib.load(fname_anat).shape
    assert np.allclose(nib.load(fname_mask).affine, nib.load(fname_anat).affine)
    assert nib.load(fname_anat).get_fdata().max() > 0
    assert nib.load(fname_mask).get_fdata()[121,126,46] == 1
    assert nib.load(fname_mask).get_fdata()[121,171,38] == 1
    
    if DEBUG:
        path_debug = os.path.join(__PATH_REPO__, "data", "debug")
        if os.path.exists(path_debug):
            shutil.rmtree(path_debug)
        shutil.copytree(__TMP_SHARE_SAVEDDATA__, path_debug, dirs_exist_ok=True)


def test_st_masking_threshold_3d_acq(tmpdir):
    fname_logfile = os.path.join(__PATH_TESTING_DATA__, "log", "test_st_masking_threshold_3d_acq.log")
    if os.path.exists(fname_logfile):
        os.remove(fname_logfile)

    config = {
        "version": "1.0",
        "parameters": {
            "method": "threshold",
            "threshold_thr": "0.1",
            "sct_deepseg_create_circular_mask": "True",
            "sct_deepseg_mask_size": "20",
            "sct_propseg_contrast": "t2s",
            "sct_propseg_include_csf": "True",
            "sendOriginal": "False",
            "bet_f": "0.5",
            "bet_g": "0.0"
        }
    }

    fname_input = os.path.join(__PATH_TESTING_DATA__, "T1w_2025-09-16-171846_81.h5")
    dset = Dataset(fname_input)
    dsetConfigAdditional = dset._dataset.require_dataset('configAdditional',shape=(1,), dtype=h5py.special_dtype(vlen=bytes))
    dsetConfigAdditional[0] = bytes(json.dumps(config), 'utf-8')
    dset.close()

    ret = subprocess.run(
        [os.path.join("/root", "shimming-toolbox", "python", "bin", "python"),
         os.path.join("/opt", "code", "python-ismrmrd-server", "client.py"),
         fname_input,
         "-p", "9020",
         "-c", "st_masking",
         "-o", tmpdir / "output.mrd",
         "-l", fname_logfile,
         ], 
         check=True,
    )

    assert os.path.exists(tmpdir / "output.mrd")
    assert ret.returncode == 0

    fname_anat = os.path.join(__TMP_SHARE_SAVEDDATA__, "anat.nii.gz")
    fname_mask = os.path.join(__TMP_SHARE_SAVEDDATA__, "mask.nii.gz")
    assert os.path.exists(fname_mask)
    assert os.path.exists(fname_anat)
    assert nib.load(fname_mask).shape == nib.load(fname_anat).shape
    assert np.allclose(nib.load(fname_mask).affine, nib.load(fname_anat).affine)

    path_debug = os.path.join(__PATH_REPO__, "data", "debug")
    if os.path.exists(path_debug):
        shutil.rmtree(path_debug)
    shutil.copytree(__TMP_SHARE_SAVEDDATA__, path_debug, dirs_exist_ok=True)

    # We run mrd2nii on the debug folder to convert the output.mrd to nifti for inspection
    subprocess.run(
        ["mrd2nii",
         "-i", __TMP_SHARE_DEBUG__,
         "-o", path_debug],
        check=True
    )

    data_output = nib.load(glob.glob(os.path.join(path_debug, "*.nii.gz"))[0]).get_fdata()
    data_mask_output = np.zeros_like(data_output)
    data_mask_output[data_output > 30000] = 1
    assert np.allclose(data_mask_output, nib.load(fname_mask).get_fdata())