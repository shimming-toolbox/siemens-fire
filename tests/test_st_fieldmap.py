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



def test_st_fieldmap_use_mask(tmpdir):
    # Copy mask into saved_data folder
    fname_mask = os.path.join(__PATH_TESTING_DATA__, "gre_mask.nii.gz")
    shutil.copy(fname_mask, os.path.join(__TMP_SHARE_SAVEDDATA__, "mask.nii.gz"))

    fname_logfile = os.path.join(__PATH_TESTING_DATA__, "test_st_fieldmap_use_mask.log")
    if os.path.exists(fname_logfile):
        os.remove(fname_logfile)

    config = {
        "version": "1.0",
        "parameters": {
            "unwrapper": "prelude",
            "gaussian-filter": "True",
            "sigma": "1",
            "use_mask": "True",
            "dilate_mask": "False",
            "dilate_mask_kernel_size": "3"
        }
    }

    fname_input = os.path.join(__PATH_TESTING_DATA__, "gre_tra_rot_2025-10-07-223140_89.h5")
    dset = Dataset(fname_input)
    dsetConfigAdditional = dset._dataset.require_dataset('configAdditional',shape=(1,), dtype=h5py.special_dtype(vlen=bytes))
    dsetConfigAdditional[0] = bytes(json.dumps(config), 'utf-8')
    dset.close()

    ret = subprocess.run(
        [os.path.join("/root", "shimming-toolbox", "python", "bin", "python"),
         os.path.join("/opt", "code", "python-ismrmrd-server", "client.py"),
         fname_input,
         "-p", "9020",
         "-c", "st_fieldmap",
         "-o", tmpdir / "output.mrd",
         "-l", fname_logfile,
         ], 
         check=True,
    )

    assert os.path.exists(tmpdir / "output.mrd")
    assert ret.returncode == 0

    fname_fieldmap = os.path.join(__TMP_SHARE_SAVEDDATA__, "fieldmap.nii.gz")
    assert os.path.exists(fname_fieldmap)

    fname_saved_mask = os.path.join(__TMP_SHARE_DEBUG__, "st_fieldmap", "saved_mask.nii.gz")
    assert os.path.exists(fname_saved_mask)

    nii_saved_mask = nib.load(fname_saved_mask)
    nii_mask = nib.load(fname_mask)
    assert nii_saved_mask.shape == nii_mask.shape
    assert np.allclose(nii_saved_mask.get_fdata(), nii_mask.get_fdata())
    
    fname_expected_fmap = os.path.join(__PATH_TESTING_DATA__, "test_st_fieldmap_use_mask", "fieldmap.nii.gz")
    fname_expected_json = os.path.join(__PATH_TESTING_DATA__, "test_st_fieldmap_use_mask", "fieldmap.json")
    assert os.path.exists(fname_expected_fmap)
    assert os.path.exists(fname_expected_json)
    nii_expected_fmap = nib.load(fname_expected_fmap)
    fname_computed_fmap = os.path.join(__TMP_SHARE_SAVEDDATA__, "fieldmap.nii.gz")
    nii_computed_fmap = nib.load(fname_computed_fmap)
    assert nii_computed_fmap.shape == nii_expected_fmap.shape
    assert np.allclose(nii_computed_fmap.get_fdata(), nii_expected_fmap.get_fdata())
    subprocess.run(['mrd2nii',
                    '-i', glob.glob(os.path.join(__TMP_SHARE_DEBUG__, "MRD_output_*.h5"))[0],
                    '-o', tmpdir],
                    check=True)
    fname_nii_output = glob.glob(os.path.join(tmpdir, "*.nii.gz"))[0]
    nii_mrd_output = nib.load(fname_nii_output)
    assert nii_mrd_output.shape == nii_expected_fmap.shape
    assert np.allclose(nii_mrd_output.get_fdata(), nii_expected_fmap.get_fdata(), atol=1.0)

    if DEBUG:
        path_debug = os.path.join(__PATH_REPO__, "data", "debug")
        if os.path.exists(path_debug):
            shutil.rmtree(path_debug)
        shutil.copytree(__TMP_SHARE_SAVEDDATA__, path_debug, dirs_exist_ok=True)
        shutil.copy(fname_saved_mask, os.path.join(path_debug, "saved_mask.nii.gz"))
        shutil.copy(fname_nii_output, os.path.join(path_debug, "mrd_output.nii.gz"))


def test_st_fieldmap_use_mask_dilate(tmpdir):
    # Copy mask into saved_data folder
    fname_mask = os.path.join(__PATH_TESTING_DATA__, "gre_mask.nii.gz")
    shutil.copy(fname_mask, os.path.join(__TMP_SHARE_SAVEDDATA__, "mask.nii.gz"))

    fname_logfile = os.path.join(__PATH_TESTING_DATA__, "test_st_fieldmap_use_mask_dilate.log")
    if os.path.exists(fname_logfile):
        os.remove(fname_logfile)

    config = {
        "version": "1.0",
        "parameters": {
            "unwrapper": "prelude",
            "gaussian-filter": "True",
            "sigma": "1",
            "use_mask": "True",
            "dilate_mask": "True",
            "dilate_mask_kernel_size": "3"
        }
    }

    fname_input = os.path.join(__PATH_TESTING_DATA__, "gre_tra_rot_2025-10-07-223140_89.h5")
    dset = Dataset(fname_input)
    dsetConfigAdditional = dset._dataset.require_dataset('configAdditional',shape=(1,), dtype=h5py.special_dtype(vlen=bytes))
    dsetConfigAdditional[0] = bytes(json.dumps(config), 'utf-8')
    dset.close()

    ret = subprocess.run(
        [os.path.join("/root", "shimming-toolbox", "python", "bin", "python"),
         os.path.join("/opt", "code", "python-ismrmrd-server", "client.py"),
         fname_input,
         "-p", "9020",
         "-c", "st_fieldmap",
         "-o", tmpdir / "output.mrd",
         "-l", fname_logfile,
         ], 
         check=True,
    )

    assert os.path.exists(tmpdir / "output.mrd")
    assert ret.returncode == 0

    fname_fieldmap = os.path.join(__TMP_SHARE_SAVEDDATA__, "fieldmap.nii.gz")
    assert os.path.exists(fname_fieldmap)

    fname_saved_mask = os.path.join(__TMP_SHARE_DEBUG__, "st_fieldmap", "saved_mask.nii.gz")
    assert os.path.exists(fname_saved_mask)

    nii_saved_mask = nib.load(fname_saved_mask)
    fname_expected_mask = os.path.join(__PATH_TESTING_DATA__, "test_st_fieldmap_use_mask_dilate", "mask_dilated.nii.gz")
    nii_mask_expected = nib.load(fname_expected_mask)
    assert nii_saved_mask.shape == nii_mask_expected.shape
    assert np.allclose(nii_saved_mask.get_fdata(), nii_mask_expected.get_fdata())

    fname_expected_fmap = os.path.join(__PATH_TESTING_DATA__, "test_st_fieldmap_use_mask_dilate", "fieldmap.nii.gz")
    fname_expected_json = os.path.join(__PATH_TESTING_DATA__, "test_st_fieldmap_use_mask_dilate", "fieldmap.json")
    assert os.path.exists(fname_expected_fmap)
    assert os.path.exists(fname_expected_json)
    nii_expected_fmap = nib.load(fname_expected_fmap)
    fname_computed_fmap = os.path.join(__TMP_SHARE_SAVEDDATA__, "fieldmap.nii.gz")
    nii_computed_fmap = nib.load(fname_computed_fmap)
    assert nii_computed_fmap.shape == nii_expected_fmap.shape
    assert np.allclose(nii_computed_fmap.get_fdata(), nii_expected_fmap.get_fdata())
    subprocess.run(['mrd2nii',
                    '-i', glob.glob(os.path.join(__TMP_SHARE_DEBUG__, "MRD_output_*.h5"))[0],
                    '-o', tmpdir],
                    check=True)
    fname_nii_output = glob.glob(os.path.join(tmpdir, "*.nii.gz"))[0]
    nii_mrd_output = nib.load(fname_nii_output)
    assert nii_mrd_output.shape == nii_expected_fmap.shape
    assert np.allclose(nii_mrd_output.get_fdata(), nii_expected_fmap.get_fdata(), atol=1.0)

    if DEBUG:
        path_debug = os.path.join(__PATH_REPO__, "data", "debug")
        if os.path.exists(path_debug):
            shutil.rmtree(path_debug)
        shutil.copytree(__TMP_SHARE_SAVEDDATA__, path_debug, dirs_exist_ok=True)
        shutil.copy(fname_saved_mask, os.path.join(path_debug, "saved_mask.nii.gz"))
        shutil.copy(fname_nii_output, os.path.join(path_debug, "mrd_output.nii.gz"))
