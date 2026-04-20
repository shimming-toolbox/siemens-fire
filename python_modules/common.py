import ismrmrd
import numpy as np
import mrd_flags

EXCLUSION_FLAGS = [
    "ISMRMRD_ACQ_IS_NOISE_MEASUREMENT",
#    "ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION",
#    "ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING",
#    "ISMRMRD_ACQ_IS_PHASE_STABILIZATION",
    "ISMRMRD_ACQ_IS_PHASE_STABILIZATION_REFERENCE",
]

KSPACE_LAYOUT = [
    "repetition",
    "phase",
    "set",
    "segment",
    "contrast", # echo
    "slice",
    "kspace_encoding_step_2",
    "kspace_encoding_step_1",
]

KSPACE_LAYOUT_IDX = [
    "repetition",
    "phase",
    "set",
    "segment",
    "contrast", # echo
    "slice",
    "kspace_encode_step_2",
    "kspace_encode_step_1",
]

class SiemensRAW:
    def __init__(self, filename: str) -> None:
        self.dset = ismrmrd.Dataset(filename, "dataset")
        self.header = ismrmrd.xsd.CreateFromDocument(self.dset.read_xml_header())
        self.n_acq = self.dset.number_of_acquisitions()

        self.kspace = None
        self.navigator = None

    def reconstruct_images():
        pass

    def build_kspace(self) -> None:
        kspace = np.zeros(self._get_kspace_dims(), dtype=np.complex128)
        navigator = np.zeros_like(kspace)

        used_idx = set()
        for i in range(self.n_acq):
            acq = self.dset.read_acquisition(i)

            if not self._check_for_flags(acq):
                continue

            # Handle reverse acquisitions (data.shape = (n_coils, samples))
            data = acq.data if not self._is_acq_reverse(acq) else acq.data[:, ::-1]
            # new shape: (samples, n_coils)
            data = data.T

            idx = (getattr(acq.idx, d) for d in KSPACE_LAYOUT_IDX)
            if idx in used_idx:
                raise IndexError(f"Index already filled: {idx}")
            used_idx.add(idx)

            if acq.is_flag_set(ismrmrd.ACQ_IS_PHASE_STABILIZATION):
                navigator[*idx, :, :] = data
            else:
                kspace[*idx, :, :] = data

        self.kspace = kspace
        self.navigator = navigator

    def load_kspace(self, filename: str) -> None:
        data = np.load(filename)
        self.kspace = data['kspace']
        self.navigator = data['navigator']

    def save_kspace(self, filename: str) -> None:
        np.savez(filename, kspace=self.kspace, navigator=self.navigator)

    def _get_kspace_dims(self):
        encoding_limits = self.header.encoding[0].encodingLimits
        get_dim = lambda x: getattr(encoding_limits, x).maximum + 1 # minimum must be zero
        nx = self.header.encoding[0].encodedSpace.matrixSize.x
        n_coils = self.header.acquisitionSystemInformation.receiverChannels

        first_dims = [get_dim(d) for d in KSPACE_LAYOUT]

        return tuple(first_dims + [nx, n_coils])

    def _check_for_flags(self, acq: ismrmrd.Acquisition) -> bool:
        for f in EXCLUSION_FLAGS:
            if acq.is_flag_set(mrd_flags.flags[f]):
                return False
        return True
    
    def _is_acq_reverse(self, acq: ismrmrd.Acquisition):
        return acq.is_flag_set(ismrmrd.ACQ_IS_REVERSE)

        