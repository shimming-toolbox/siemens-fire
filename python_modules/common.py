import ismrmrd
import numpy as np
import itertools
from ismrmrd import constants
from collections import defaultdict

EXCLUSION_FLAGS = [
    constants.ACQ_IS_NOISE_MEASUREMENT,
    constants.ACQ_IS_PHASE_STABILIZATION_REFERENCE
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
        n_echo = mrd_header.encoding[0].encodingLimits.contrast.maximum + 1

        self.acquisitions: list[ismrmrd.Acquisition] = []
        self.acs_mask = None # assuming they are the same for every raw kspace
        self.kspace = None
        self.navigator = None
        self.noise = None

        self.navigator_detector = NavigatorDetector(n_echo)

    def reconstruct_images():
        pass

    def add_acq(self, acq: ismrmrd.Acquisition) -> None:
        self.acquisitions.append(acq)

    def reset_acq(self) -> None:
        self.acquisitions = []

    def extract_noise(self) -> None:
        # Assume noise is first acquisition
        if self.acquisitions[0].is_flag_set(ismrmrd.ACQ_IS_NOISE_MEASUREMENT):
            self.noise = self.acquisitions[0]
            self.acquisitions = self.acquisitions[1:]
        else:
            raise ValueError("Noise acquisition not found")
    
    def set_noise_acq(self, noise_acq) -> None:
        self.noise = noise_acq

    def remove_phase_stabilization_references(self) -> None:
        # At beginning of each repetition, you have n_slices * (n_echo + 1)
        # of ISMRMRD_ACQ_IS_PHASE_STABILIZATION_REFERENCE acquisitions. The +1
        # is for the phase stabilization acq, labeled as echo=0.

        _, _, _, _, n_echoes, n_slices, _, _, _, _ = self._get_kspace_dims()
        n_echoes += 1 # fifth echo for phase stabilization

        # TODO: verify this for each repetition
        a = 0
        b = a + n_echoes*n_slices
        indices_to_remove = range(a, b)
        indices_to_remove = set(indices_to_remove) # for constant lookup

        self.acquisitions = [acq for i, acq in enumerate(self.acquisitions) if i not in indices_to_remove]

    def build_kspace(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        dims = self._get_kspace_dims()

        kspace = np.zeros(dims, dtype=np.complex64)
        # Navigator has same shape, except for the echoes
        navigator = np.zeros(dims[:4] + (1,) + dims[5:], dtype=np.complex64)
        # Only need one contrast dimension for navigator
        #navigator = navigator[:, :, :, :, [0], :, :, :, :, :]
        acs_mask = np.zeros_like(kspace, dtype=np.bool)

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

            if self.navigator_detector.is_nav(acq):
                navigator[*idx, :, :] = data
            else:
                kspace[*idx, :, :] = data
                if acq.is_flag_set(ismrmrd.ACQ_IS_PARALLEL_CALIBRATION) or acq.is_flag_set(ismrmrd.ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING):
                    acs_mask[*idx, :, :] = True

        #self.kspace = kspace
        #self.navigator = navigator
        #self.acs_mask = acs_mask

        return kspace, navigator, acs_mask

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
            if acq.is_flag_set(f):   # f is now the integer value directly, not a string key
                return False
        return True
    
    def _is_acq_reverse(self, acq: ismrmrd.Acquisition):
        return acq.is_flag_set(ismrmrd.ACQ_IS_REVERSE)
    
# Navigator line detection for the following sequence RF-E1-E2-E3-E4-NAV
class NavigatorDetector:
    """
    Detects navigator acquisitions based on their position within the shot cycle.

    The sequence structure is fixed : RF - E1 - E2 - E3 - E4 - NAV
    meaning for every group of (nEcho + 1) acquisitions sharing the same
    (slice, repetition), the last one is the navigator.

    This position-based detection is used because the navigator acquisition
    does not carry a dedicated ISMRMRD flag in this sequence. It is identified
    purely by its index within the shot cycle.

    The counter is maintained per (slice, repetition) key because acquisitions
    from different slices are interleaved in the dataset and must be tracked independently.
    """
    def __init__(self, n_echo):
        self.n_echo = n_echo # Number of echoes per shot: defines the cycle length (n_echo + 1)
        self.counter = defaultdict(int) # Per(slice, rep) acquisition counter: tracks position within the shot cycle
        self.nav_per_slice = defaultdict(int)

    def is_nav(self, acq):
        """
        Determine whether an acquisition is a navigator based on its position
        in the shot cycle for its (slice, repetition) group.

        Shot cycle for n_echo=4 :
            counter mod 5 == 0 → echo 1  (imaging)
            counter mod 5 == 1 → echo 2  (imaging)
            counter mod 5 == 2 → echo 3  (imaging)
            counter mod 5 == 3 → echo 4  (imaging)
            counter mod 5 == 4 → NAV     (navigator)  ← (counter - 1) % 5 == 4 == n_echo

        Parameters
        ----------
        acq (ismrmrd.Acquisition): Acquisition object to classify.

        Returns
        -------
        (bool): True if this acquisition is the navigator of its shot, False if imaging echo.
        """
        # Group acquisitions by (slice, repetition): each group has its own cycle counter
        sl = acq.idx.slice
        rep = acq.idx.repetition
        key = (sl, rep)

        # Increment counter BEFORE the modulo check, then use (counter - 1) so that
        # the first acquisition of a group (counter=1) maps to position 0 in the cycle
        self.counter[key] += 1

        # Navigator is at position n_echo in the 0-based cycle (last of each n_echo+1 group)
        if (self.counter[key] - 1) % (self.n_echo + 1) == self.n_echo:
            self.nav_per_slice[sl] += 1
            return True

        return False

        