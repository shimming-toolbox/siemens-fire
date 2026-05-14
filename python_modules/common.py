import ismrmrd
import numpy as np
import mrd_flags
import itertools

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
    def __init__(self, mrd_header) -> None:
        #self.dset = ismrmrd.Dataset(filename, "dataset")
        self.header = mrd_header
        #self.n_acq = self.dset.number_of_acquisitions()

        self.acquisitions: list[ismrmrd.Acquisition] = []
        self.acs_mask = None # assuming they are the same for every raw kspace
        self.kspace = None
        self.navigator = None
        self.noise = None

    def reconstruct_images():
        pass

    def add_acq(self, acq: ismrmrd.Acquisition) -> None:
        self.acquisitions.append(acq)

    def extract_noise(self) -> None:
        # Assume noise is first acquisition
        if self.acquisitions[0].is_flag_set(ismrmrd.ACQ_IS_NOISE_MEASUREMENT):
            self.noise = self.acquisitions[0]
            self.acquisitions = self.acquisitions[1:]
        else:
            raise ValueError("Noise acquisition not found")

    def remove_phase_stabilization_references(self) -> None:
        # At beginning of each repetition, you have n_slices * (n_echo + 1)
        # of ISMRMRD_ACQ_IS_PHASE_STABILIZATION_REFERENCE acquisitions. The +1
        # is for the phase stabilization acq, labeled as echo=0.

        k_space_encode_1_steps = set()
        for acq in self.acquisitions:
            k_space_encode_1_steps.add(acq.getHead().idx.kspace_encode_step_1)

        ny = len(k_space_encode_1_steps)

        n_reps, _, _, _, n_echoes, n_slices, _, _, _, _ = self._get_kspace_dims()
        n_echoes += 1 # fifth echo for phase stabilization

        indices_to_remove = range(0)
        for n in range(n_reps):
            a = n * (n_echoes*n_slices*ny) + n * (n_echoes*n_slices)
            b = a + n_echoes*n_slices
            indices_to_remove = itertools.chain(indices_to_remove, range(a, b))
        indices_to_remove = set(indices_to_remove) # for constant lookup

        self.acquisitions = [acq for i, acq in enumerate(self.acquisitions) if i not in indices_to_remove]

    def build_kspace(self) -> None:
        kspace = np.zeros(self._get_kspace_dims(), dtype=np.complex128)
        navigator = np.zeros_like(kspace)
        acs_mask = np.zeros_like(kspace)

        used_idx = set()
        for i in range(len(self.acquisitions)):
            acq = self.acquisitions[i]

            if not self._check_for_flags(acq):
                continue

            # Handle reverse acquisitions (data.shape = (n_coils, samples))
            data = acq.data if not self._is_acq_reverse(acq) else acq.data[:, ::-1]
            # new shape: (samples, n_coils)
            data = data.T

            idx = tuple(getattr(acq.idx, d) for d in KSPACE_LAYOUT_IDX)
            if idx in used_idx and not (acq.idx.contrast == 0 and acq.scan_counter % 5 == 1):
                raise IndexError(f"Index already filled: {idx}")
            used_idx.add(idx)

            if acq.idx.contrast == 0 and acq.scan_counter % 5 == 1: # TODO: double check for every repetition
                navigator[*idx, :, :] = data
            else:
                kspace[*idx, :, :] = data
                if acq.is_flag_set(ismrmrd.ACQ_IS_PARALLEL_CALIBRATION) or acq.is_flag_set(ismrmrd.ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING):
                    acs_mask[*idx, :, :] = True

        self.kspace = kspace
        self.navigator = navigator
        self.acs_mask = acs_mask

    def load_kspace(self, filename: str) -> None:
        data = np.load(filename)
        self.kspace = data['kspace']
        self.navigator = data['navigator']
        self.acs_mask = data['acs_mask']

    def save_kspace(self, filename: str) -> None:
        np.savez(filename, kspace=self.kspace, navigator=self.navigator, acs_mask=self.acs_mask)

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

        