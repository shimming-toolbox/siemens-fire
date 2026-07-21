import os
import shutil
import subprocess
import json
from ismrmrd import Dataset
import h5py

from python_modules import __PATH_TESTING_DATA__, __PATH_REPO__, __TMP_SHARE_DEBUG__, __TMP_SHARE_SAVEDDATA__


def test_nomad(tmpdir):
    # Debug does not work well. There seems to be an interaction between how NOMAD instances are killed and the debugger.

    # Copy mask into saved_data folder
    fname_mask = os.path.join(__PATH_TESTING_DATA__, "gre_mask.nii.gz")
    shutil.copy(fname_mask, os.path.join(__TMP_SHARE_SAVEDDATA__, "mask.nii.gz"))

    fname_logfile = os.path.join(__PATH_TESTING_DATA__, "test_nomad.log")
    if os.path.exists(fname_logfile):
        os.remove(fname_logfile)

    config = {
        "version": "1.0",
        "parameters": {
        "channels_to_shim": "z",
        "objective": "Sig int",
        "use_surrogate": "False",
        "use_f0_offset_from_gradients": "True",
        "f0_bounds": [-100, 100],
        "xyz_bounds": [-0.2, 0.2]
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
         "-c", "nomad_f0xyz_opt",
         "-o", tmpdir / "output.mrd",
         "-l", fname_logfile,
         ], 
         check=True,
    )

    assert os.path.exists(tmpdir / "output.mrd")
    assert ret.returncode == 0

    fname_shim_file = os.path.join(__TMP_SHARE_SAVEDDATA__, "scanner_shim.txt")
    assert os.path.exists(fname_shim_file)
