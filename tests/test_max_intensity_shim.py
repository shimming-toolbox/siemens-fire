import os
import shutil
import subprocess
import json
from ismrmrd import Dataset
import h5py

from python_modules import __PATH_TESTING_DATA__, __PATH_REPO__, __TMP_SHARE_DEBUG__, __TMP_SHARE_SAVEDDATA__


def test_max_intensity_shim(tmpdir):

    # Copy mask into saved_data folder
    fname_mask = os.path.join(__PATH_TESTING_DATA__, "gre_mask.nii.gz")
    shutil.copy(fname_mask, os.path.join(__TMP_SHARE_SAVEDDATA__, "mask.nii.gz"))

    # Create dummy txt file containing shim values to be used for the test
    fname_shim_file = os.path.join(__TMP_SHARE_SAVEDDATA__, "scanner_shim_last_acq.txt")
    with open(fname_shim_file, "w") as f:
        f.write("f0 | Gx | Gy | Gz\n")
        for i_vol in range(100):
            for i_slice in range(20):
                idx = (i_vol * 1000) + (i_slice * 10)
                f.write(f"{idx} | {idx + 1} | {idx + 2} | {idx + 3}\n")

    fname_logfile = os.path.join(__PATH_TESTING_DATA__, "log", "test_max_intensity_shim.log")
    if os.path.exists(fname_logfile):
        os.remove(fname_logfile)

    config = {
        "version": "1.0",
        "parameters": {
        "use_already_computed_mask": "True",
        }
    }

    fname_input = os.path.join(__PATH_TESTING_DATA__, "ep2d_bold_shimming_nomad_opt_surr_fz_2025-10-23-163657_98.h5")
    dset = Dataset(fname_input)
    dsetConfigAdditional = dset._dataset.require_dataset('configAdditional',shape=(1,), dtype=h5py.special_dtype(vlen=bytes))
    dsetConfigAdditional[0] = bytes(json.dumps(config), 'utf-8')
    dset.close()

    ret = subprocess.run(
        [os.path.join("/root", "shimming-toolbox", "python", "bin", "python"),
         os.path.join("/opt", "code", "python-ismrmrd-server", "client.py"),
         fname_input,
         "-p", "9020",
         "-c", "max_intensity_shim",
         "-o", tmpdir / "output.mrd",
         "-l", fname_logfile,
         ], 
         check=True,
    )

    assert os.path.exists(tmpdir / "output.mrd")
    assert ret.returncode == 0

    fname_shim_file = os.path.join(__TMP_SHARE_SAVEDDATA__, "scanner_shim.txt")
    assert os.path.exists(fname_shim_file)
    with open(fname_shim_file, "r") as f:
        lines = f.readlines()
        assert lines[0].strip() == "93000.000000 | 93001.000000 | 93002.000000 | 93003.000000"
