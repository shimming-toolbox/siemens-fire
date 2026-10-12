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



def test_st_b0shim_gre_slicewise(tmpdir):
    # Copy mask into saved_data folder
    fname_mask = os.path.join(__PATH_TESTING_DATA__, "gre_mask.nii.gz")
    shutil.copy(fname_mask, os.path.join(__TMP_SHARE_SAVEDDATA__, "mask.nii.gz"))
    fname_fieldmap = os.path.join(__PATH_TESTING_DATA__, "fieldmap.nii.gz")
    shutil.copy(fname_fieldmap, os.path.join(__TMP_SHARE_SAVEDDATA__, "fieldmap.nii.gz"))
    shutil.copy(fname_fieldmap.replace(".nii.gz", ".json"), os.path.join(__TMP_SHARE_SAVEDDATA__, "fieldmap.json"))

    fname_logfile = os.path.join(__PATH_TESTING_DATA__, "log", "test_st_b0shim_gre_slicewise.log")
    if os.path.exists(fname_logfile):
        os.remove(fname_logfile)

    config = {
        "version": "1.0",
        "parameters": {
            "shim_method": "slicewise",
            "scanner-coil-order": "0,1",
            "optimizer-method": "pseudo_inverse",
            "optimizer-criteria": "rmse",
            "regularization-factor": "0.1",
            "weighting-signal-loss": "10",
            "mask-dilation-kernel-size": "5",
            "off-channels": "1,2"
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
         "-c", "st_b0shim",
         "-o", tmpdir / "output.mrd",
         "-l", fname_logfile,
         ], 
         check=True,
    )

    assert os.path.exists(tmpdir / "output.mrd")
    assert ret.returncode == 0

    fname_shim_file = os.path.join(__TMP_SHARE_SAVEDDATA__, "scanner_shim.txt")
    assert os.path.exists(fname_shim_file)

    nii_target = nib.load(glob.glob(os.path.join(__TMP_SHARE_DEBUG__, "st_b0shim", "mrd2nii", "*nii.gz"))[0])

    with open(fname_shim_file, 'r') as f:
        lines = f.readlines()
        assert len(lines) == nii_target.shape[2] + 1
        assert lines[1] == " 17.942297 | 0.000000 | 0.000000 | -0.000602\n"


def test_st_b0shim_gre_volume(tmpdir):
    # Copy mask into saved_data folder
    fname_mask = os.path.join(__PATH_TESTING_DATA__, "gre_mask.nii.gz")
    shutil.copy(fname_mask, os.path.join(__TMP_SHARE_SAVEDDATA__, "mask.nii.gz"))
    fname_fieldmap = os.path.join(__PATH_TESTING_DATA__, "fieldmap.nii.gz")
    shutil.copy(fname_fieldmap, os.path.join(__TMP_SHARE_SAVEDDATA__, "fieldmap.nii.gz"))
    shutil.copy(fname_fieldmap.replace(".nii.gz", ".json"), os.path.join(__TMP_SHARE_SAVEDDATA__, "fieldmap.json"))

    fname_logfile = os.path.join(__PATH_TESTING_DATA__, "log", "test_st_b0shim_gre_volume.log")
    if os.path.exists(fname_logfile):
        os.remove(fname_logfile)

    config = {
        "version": "1.0",
        "parameters": {
            "shim_method": "volume",
            "scanner-coil-order": "0,1,2",
            "optimizer-method": "pseudo_inverse",
            "optimizer-criteria": "rmse",
            "regularization-factor": "0.1",
            "weighting-signal-loss": "10",
            "mask-dilation-kernel-size": "5",
            "off-channels": "1,2"
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
         "-c", "st_b0shim",
         "-o", tmpdir / "output.mrd",
         "-l", fname_logfile,
         ], 
         check=True,
    )

    assert os.path.exists(tmpdir / "output.mrd")
    assert ret.returncode == 0

    fname_shim_file = os.path.join(__TMP_SHARE_SAVEDDATA__, "scanner_shim.txt")
    assert os.path.exists(fname_shim_file)
    fname_adjvalidate_commands = os.path.join(__TMP_SHARE_SAVEDDATA__, "commands_to_apply_coefs_at_the_scanner.txt")
    assert os.path.exists(fname_adjvalidate_commands)

    with open(fname_shim_file, 'r') as f:
        lines = f.readlines()
        assert len(lines) == 1
        assert lines[0] == "13.685491, 0.000000, 0.000000, 1.459667, -7.421840, 0.250997, -0.545989, 0.141669, -0.003329, \n"


def test_st_b0shim_gre_slicewise_gradient_cs(tmpdir):
    # Copy mask into saved_data folder
    fname_mask = os.path.join(__PATH_TESTING_DATA__, "gre_mask.nii.gz")
    shutil.copy(fname_mask, os.path.join(__TMP_SHARE_SAVEDDATA__, "mask.nii.gz"))
    fname_fieldmap = os.path.join(__PATH_TESTING_DATA__, "fieldmap.nii.gz")
    shutil.copy(fname_fieldmap, os.path.join(__TMP_SHARE_SAVEDDATA__, "fieldmap.nii.gz"))
    shutil.copy(fname_fieldmap.replace(".nii.gz", ".json"), os.path.join(__TMP_SHARE_SAVEDDATA__, "fieldmap.json"))

    fname_logfile = os.path.join(__PATH_TESTING_DATA__, "log", "test_st_b0shim_gre_slicewise_gradient_cs.log")
    if os.path.exists(fname_logfile):
        os.remove(fname_logfile)

    config = {
        "version": "1.0",
        "parameters": {
            "shim_method": "slicewise",
            "scanner-coil-order": "0,1",
            "opt-cs": "gradient-cs",
            "optimizer-method": "lin_lsq",
            "optimizer-criteria": "rmse",
            "regularization-factor": "0.1",
            "weighting-signal-loss": "10",
            "mask-dilation-kernel-size": "5",
            "off-channels": "1,2"
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
         "-c", "st_b0shim",
         "-o", tmpdir / "output.mrd",
         "-l", fname_logfile,
         ], 
         check=True,
    )

    assert os.path.exists(tmpdir / "output.mrd")
    assert ret.returncode == 0

    fname_shim_file = os.path.join(__TMP_SHARE_SAVEDDATA__, "scanner_shim.txt")
    assert os.path.exists(fname_shim_file)

    nii_target = nib.load(glob.glob(os.path.join(__TMP_SHARE_DEBUG__, "st_b0shim", "mrd2nii", "*nii.gz"))[0])

    with open(fname_shim_file, 'r') as f:
        lines = f.readlines()
        assert len(lines) == nii_target.shape[2] + 1
        assert lines[1] == " 16.351514 | -0.000126 | -0.000128 | -0.000713\n"


def test_st_b0shim_gre_slicewise_convert(tmpdir):
    # Copy mask into saved_data folder
    fname_mask = os.path.join(__PATH_TESTING_DATA__, "gre_mask.nii.gz")
    shutil.copy(fname_mask, os.path.join(__TMP_SHARE_SAVEDDATA__, "mask.nii.gz"))
    fname_fieldmap = os.path.join(__PATH_TESTING_DATA__, "fieldmap.nii.gz")
    shutil.copy(fname_fieldmap, os.path.join(__TMP_SHARE_SAVEDDATA__, "fieldmap.nii.gz"))
    shutil.copy(fname_fieldmap.replace(".nii.gz", ".json"), os.path.join(__TMP_SHARE_SAVEDDATA__, "fieldmap.json"))

    shutil.copy(os.path.join(__PATH_TESTING_DATA__, "coefs_gradientcs.txt"), os.path.join(__TMP_SHARE_SAVEDDATA__, "coefs_gradientcs.txt"))

    fname_logfile = os.path.join(__PATH_TESTING_DATA__, "log", "test_st_b0shim_gre_slicewise_convert.log")
    if os.path.exists(fname_logfile):
        os.remove(fname_logfile)

    config = {
        "version": "1.0",
        "parameters": {
            "shim_method": "slicewise",
            "scanner-coil-order": "0,1",
            "optimizer-method": "pseudo_inverse",
            "optimizer-criteria": "rmse",
            "regularization-factor": "0.1",
            "weighting-signal-loss": "10",
            "mask-dilation-kernel-size": "5",
            "off-channels": "1,2"
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
         "-c", "st_b0shim",
         "-o", tmpdir / "output.mrd",
         "-l", fname_logfile,
         ], 
         check=True,
    )

    assert os.path.exists(tmpdir / "output.mrd")
    assert ret.returncode == 0

    fname_converted_file_nof0 = os.path.join(__TMP_SHARE_SAVEDDATA__, "coefs_shimcs.txt")
    assert os.path.exists(fname_converted_file_nof0)
    fname_converted_file_f0 = os.path.join(__TMP_SHARE_SAVEDDATA__, "coefs_shimcs_f0.txt")
    assert os.path.exists(fname_converted_file_f0)

    with open(fname_converted_file_nof0, 'r') as f:
        lines = f.readlines()
        assert lines[1] == "0.000000 | 0.035912 | 0.036466 | 0.203668\n"

    with open(fname_converted_file_f0, 'r') as f:
        lines = f.readlines()
        assert lines[1] == "-879.359762 | 0.035912 | 0.036466 | 0.203668\n"
