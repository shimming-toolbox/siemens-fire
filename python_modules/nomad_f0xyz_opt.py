import copy
import ctypes
from dataclasses import dataclass
from functools import partial
import ismrmrd
from joblib import Parallel, delayed
import json
import logging
from multiprocessing import Process, Queue, JoinableQueue
import nibabel as nib
from nibabel.processing import resample_from_to as nib_resample_from_to
import numpy as np
import os
import PyNomad
from scipy.ndimage import binary_dilation, binary_erosion, generate_binary_structure, iterate_structure, center_of_mass
import signal
import shutil
from threading import Thread
import traceback
import mrdhelper # Custom module for MRD helper functions found in the python-ismrmrd-server repository

from skimage.metrics import structural_similarity as ssim
from skimage.metrics import normalized_mutual_information as nmi
from scipy.ndimage import center_of_mass

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable

from mrd2nii.mrd2nii_main import mrd2nii_volume, mrd2nii_stack
from shimmingtoolbox.coils.coil import ScannerCoil


# Folder for debug output files
debugFolder = None
dataFolder = os.path.abspath("/tmp/share/saved_data")
fname_all_currents = None
fname_currents_per_volume = None
fname_solution = None


NB_CHANNELS_TO_SEND_TO_SEQ = 4
f0_ub = 100
f0_lb = -100
gradx_ub = 0.2
gradx_lb = -0.2
grady_ub = 0.2
grady_lb = -0.2
gradz_ub = 0.2
gradz_lb = -0.2
sum_cstr = np.inf

# Signal interrupt for bb_block
is_sigint = False

# Surrogate global variables
is_surrogate_ready = False
surr_mask = None  # 4d fmap coord (X, Y, Z, ishim)
time_ordering_to_nifti_slices = None

# Obj with sig int
is_obj_with_sigint = False

# Obj with mutual information
is_obj_with_mi = False
# mi_metric = 'ssim'
mi_metric = 'nmi'
if mi_metric == 'nmi':
    HYPERPARAM_MI = 150.0
elif mi_metric == 'ssim':
    HYPERPARAM_MI = 80.0
else:
    raise ValueError(f"mi_metric {mi_metric} not recognized")


# obj with fmap
is_obj_with_fmap = False
is_obj_with_fmap_ready = None
obj_mask = None
obj_time_ordering_to_nifti_slices = None
list_max_rmse = []
MAX_RMSE_SCALAR = 1.2

# Slice location global variables
is_slice_loc_ready = False
slice_loc = None
target_iso = None


@dataclass
class Currents:
    currents: list
    slice: int
    f: float = None
    repetition: int = None


class MyFeedbackData(ctypes.Structure):
    _pack_ = 1                            # Align on 1 byte boundaries
    _fields_ = [
        ('slice_nb', ctypes.c_int32),                #          
        ('currents',  ctypes.c_float * NB_CHANNELS_TO_SEND_TO_SEQ) #          NB_CHANNELS_TO_SEND_TO_SEQ bytes
    ]                                               #       =  NB_CHANNELS_TO_SEND_TO_SEQ + 1 bytes total


class ReturnThread(Thread):
    def __init__(self, group=None, target=None, name=None, args=(), kwargs=None, daemon=None):
        super().__init__(group=group, target=target, name=name, args=args, kwargs=kwargs, daemon=daemon)
        self._return = None

    def run(self):
        if self._target is not None:
            self._return = self._target(*self._args, **self._kwargs)

    def join(self):
        super().join()
        return self._return


def find_and_read_anat_within_data():
    # For now, take the latest created mask in dataFolder
    anat_files = [f for f in os.listdir(dataFolder) if "anat" in f and f.endswith('.nii.gz')]
    if not anat_files:
        raise RuntimeError(f"No anat files found. Is it in {dataFolder}?")

    latest_anat_file = max(anat_files, key=lambda f: os.path.getctime(os.path.join(dataFolder, f)))
    logging.info(f"Using latest anat file: {latest_anat_file}")
    return nib.load(os.path.join(dataFolder, latest_anat_file))


def find_and_read_mask_within_data():
    # For now, take the latest created mask in dataFolder
    mask_files = [f for f in os.listdir(dataFolder) if "mask" in f and f.endswith('.nii.gz')]
    if not mask_files:
        raise RuntimeError(f"No mask files found. Is it in {dataFolder}?")

    latest_mask_file = max(mask_files, key=lambda f: os.path.getctime(os.path.join(dataFolder, f)))
    logging.info(f"Using latest mask file: {latest_mask_file}")
    return nib.load(os.path.join(dataFolder, latest_mask_file))


def find_and_read_fmap_within_data():
    # For now, take the latest created fieldmap in dataFolder
    fieldmap_files = [f for f in os.listdir(dataFolder) if "fieldmap" in f and f.endswith('.nii.gz')]
    if not fieldmap_files:
        raise RuntimeError(f"No fieldmap files found. Is it in {dataFolder}?")

    latest_fieldmap_file = max(fieldmap_files, key=lambda f: os.path.getctime(os.path.join(dataFolder, f)))
    logging.info(f"Using latest fieldmap file: {latest_fieldmap_file}")
    nii = nib.load(os.path.join(dataFolder, latest_fieldmap_file))

    json_path = os.path.join(dataFolder, latest_fieldmap_file.replace('.nii.gz', '.json'))
    with open(json_path, 'r') as f:
        json_data = json.load(f)

    return nii, json_data


def process_print_queue(print_queue):

    messages = []
    # Process the print queue
    while True:
        message = print_queue.get()
        if message is None:
            print_queue.task_done()
            break
        messages.append(message)
        print_queue.task_done()

    return messages

def process_worker_results(result_queue):
    currents = []
    while True:
        current = result_queue.get()
        if current is None:
            result_queue.task_done()
            break
        currents.append(current)
        result_queue.task_done()

    return currents

def send_coefs_to_scanner(send_to_scanner_queue, print_queue, connection):
    print_queue.put("Starting thread to send coefs to Scanner (in thread)")
    n_currents = 0
    
    while True:
        currents = send_to_scanner_queue.get()
        
        # None is sent if the thread is stopped
        if currents is not None:
            n_currents += 1
            print_queue.put(f"Sending currents to scanner: {currents.currents}, slice: {currents.slice}")
            feedback_data = MyFeedbackData()
            feedback_data.slice_nb = currents.slice
            for i in range(NB_CHANNELS_TO_SEND_TO_SEQ):
                feedback_data.currents[i] = currents.currents[i]

            try:
                connection.send_feedback("nomad_current", feedback_data)
            except BrokenPipeError as e:
                print_queue.put(f"BrokenPipeError when sending currents to the scanner: {e}")

        else:
            print_queue.put(f"Number of currents sent to the scanner: {n_currents}")
            break

def stop(send_to_scanner_thread, nomad_instance_stopped, workers, result_queue, result_thread,
         print_queue, print_queue_thread, nb_slices, nb_repetitions, chronological_to_nifti_ordering, fatsat, send_to_scanner_queue, f_queues, channels_to_shim):
    
    def process_results(currents, nb_slices, nb_repetitions, chronological_to_nifti_ordering):

        logging.info("Processing results")
        opt_currents = np.zeros((nb_slices, NB_CHANNELS_TO_SEND_TO_SEQ))
        max_sig_int = np.zeros((nb_slices,))
        logging.info(f"Number of processed currents received in result_queue: {len(currents)}")

        all_currents = np.zeros((nb_slices, nb_repetitions + 1, NB_CHANNELS_TO_SEND_TO_SEQ + 1))

        for current in currents:
            if current.f > max_sig_int[current.slice]:
                max_sig_int[current.slice] = current.f
                opt_currents[current.slice, :] = current.currents
            all_currents[current.slice, current.repetition, :-1] = current.currents
            all_currents[current.slice, current.repetition, -1] = current.f

        currents_per_volume = np.empty(shape=(nb_slices, nb_repetitions + 1, NB_CHANNELS_TO_SEND_TO_SEQ + 1))
        currents_per_volume[:] = np.nan

        for current in currents:
            a_cur = copy.deepcopy(current.currents)
            a_cur.append(current.f)
            nifti_slice = chronological_to_nifti_ordering[current.slice][0]
            currents_per_volume[nifti_slice, current.repetition, :] = a_cur

        logging.info("Processed currents")
        np.save(fname_all_currents, all_currents)
        np.save(fname_currents_per_volume, currents_per_volume)
        return opt_currents
    
    logging.info("Stopping NomadOpt")
    if send_to_scanner_thread is not None:
        logging.info("Closing Scanner queue")
        send_to_scanner_queue.put(None)
        logging.info("Stopping Scanner Process")
        send_to_scanner_thread.join()
    else:
        logging.info("Scanner thread was not started")

    logging.info("Stopping all queues")
    for i, q in enumerate(f_queues):
        if not nomad_instance_stopped[i]:
            q.put(None)

    extra_nones = 5
    for i, q in enumerate(f_queues):
        for _ in range(extra_nones):
            q.put(None)

    logging.info("Stopping worker thread")
    for i, worker in enumerate(workers):
        logging.info(f"Stopped worker {i}")
        worker.join()

    for i, q in enumerate(f_queues):
        cnt = 0
        while not q.empty():
            logging.info(f"{q.get()=}")
            cnt += 1
        logging.info(f"Stopped queue {i} took {extra_nones - cnt} extra None(s)")

    result_queue.put(None)
    result_currents = result_thread.join()
    opt_currents = process_results(result_currents, nb_slices, nb_repetitions, chronological_to_nifti_ordering)

    def write_text_file(data, nb_slices, fatsat):
        logging.info("Writing currents to text file")

        def parse_chrono_to_nifti_2_nifti_to_chrono(chrono_to_nifti):
            def find_max_slice(chrono_to_nifti):
                max_slice = -1
                for shim_group in chrono_to_nifti:
                    if max(shim_group) > max_slice:
                        max_slice = max(shim_group)
                return max_slice

            n_slices = find_max_slice(chrono_to_nifti) + 1
            nifti_to_chrono = []
            for i in range(n_slices):
                for i_sg, shim_group in enumerate(chrono_to_nifti):
                    if i in shim_group:
                        nifti_to_chrono.append(i_sg)
                        break
            return nifti_to_chrono
        
        # Output file is nifti ordered, not chronologically ordered so we need to convert
        nifti_to_chrono = parse_chrono_to_nifti_2_nifti_to_chrono(chronological_to_nifti_ordering)

        with open(fname_solution, 'w') as f:
            for i_slice_nifti_idx in range(nb_slices):
                a_slice_chrono_idx = nifti_to_chrono[i_slice_nifti_idx]
                if fatsat:
                    for i_channel in range(NB_CHANNELS_TO_SEND_TO_SEQ):
                        f.write(f"0.0, ")
                    f.write("\n")
                for i_channel in range(NB_CHANNELS_TO_SEND_TO_SEQ):
                    f.write(f"{data[a_slice_chrono_idx, i_channel]}, ")
                f.write("\n")
        logging.info(f"Currents written to {fname_solution}")

    # Write the currents to a text file
    write_text_file(opt_currents, nb_slices, fatsat)
    shutil.copyfile(fname_solution, os.path.join(dataFolder, "scanner_shim.txt"))

    # Create interactive plot of the objective function and the currents
    opt_channels = [channel for channel in channels_to_shim]
    create_obj_interactive_plot(fname_currents_per_volume, path_output=os.path.join(debugFolder, 'obj_plots'), opt_channels=opt_channels)

    logging.info("Processing print queue")
    print_queue.put(None)
    messages = print_queue_thread.join()
    for message in messages:
        logging.info(message)
    print_queue.join()
    logging.info(opt_currents)


def process(connection, config, metadata):

    logging.info("Optimizing")
    is_processing = True
    
    logging.info(f"Config: {config}\n")
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

    use_surrogate = mrdhelper.get_json_config_param(config_dict, "use_surrogate", default=False, type="bool")
    channels_to_shim = mrdhelper.get_json_config_param(config_dict, "channels_to_shim", default="fxyz")
    objective = mrdhelper.get_json_config_param(config_dict, "objective", default="Sig int")
    use_f0_offset_from_gradients = mrdhelper.get_json_config_param(config_dict, "use_f0_offset_from_gradients", default=True, type="bool")
    # Todo TEMP: Remove this
    #objective = "Sig int + mi"
    #channels_to_shim = 'fz'
    #use_surrogate = False
    #use_f0_offset_from_gradients = True

    global debugFolder, fname_solution, is_obj_with_sigint, is_obj_with_mi, is_obj_with_fmap, fname_all_currents, fname_currents_per_volume
    is_obj_with_sigint = False
    is_obj_with_mi = False
    is_obj_with_fmap = False
    debugFolder = os.path.abspath(f"/tmp/share/debug/nomad_files_{channels_to_shim}")
    if use_surrogate:
        debugFolder += "_surrogate"
    
    found_objective = False
    if "SIG INT" in objective.upper():
        debugFolder += "_si"
        found_objective = True
        is_obj_with_sigint = True
    if "MI" in objective.upper():
        debugFolder += "_mi"
        found_objective = True
        is_obj_with_mi = True
        logging.info("Using MI objective function with metric: %s", mi_metric)
    if "FMAP" in objective.upper():
        debugFolder += "_fmap"
        found_objective = True
        is_obj_with_fmap = True
    
    if not found_objective:
        raise ValueError(f"Objective '{objective}' not recognized. Available objectives are: 'sig int', 'mi', 'fmap', or combinations of these separated by +")
    if not is_obj_with_sigint and not is_obj_with_mi:
        raise ValueError("At least one of 'sig int', 'mi' objective function should be selected")
    
    if use_f0_offset_from_gradients:
        debugFolder += "_f0_offset"
    if os.path.exists(debugFolder):
        shutil.rmtree(debugFolder)
    os.makedirs(debugFolder)
    fname_solution = os.path.join(debugFolder, "scanner_shim.txt")
    fname_all_currents = os.path.join(debugFolder, f"all_currents.npy")
    fname_currents_per_volume = os.path.join(debugFolder, f"currents_per_volume.npy")

    # Metadata should be MRD formatted header, but may be a string
    try:
        # Disabled due to incompatibility between PyXB and Python 3.8:
        # https://github.com/pabigot/pyxb/issues/123
        # # logging.info("Metadata: \n%s", metadata.toxml('utf-8'))

        logging.info("Incoming dataset contains %d encodings", len(metadata.encoding))
        logging.info(
            "First encoding is of type '%s', with a matrix size of (%s x %s x %s) and a field of view of (%s x %s x %s)mm^3",
            metadata.encoding[0].trajectory,
            metadata.encoding[0].encodedSpace.matrixSize.x,
            metadata.encoding[0].encodedSpace.matrixSize.y,
            metadata.encoding[0].encodedSpace.matrixSize.z,
            metadata.encoding[0].encodedSpace.fieldOfView_mm.x,
            metadata.encoding[0].encodedSpace.fieldOfView_mm.y,
            metadata.encoding[0].encodedSpace.fieldOfView_mm.z)

        # Extract ordering
        mrd_idx_to_order_idx = {}
        order_idx_to_mrd_idx = {}
        nb_slices = 0

        for param in metadata.userParameters.userParameterLong:
            if param.name.startswith("RelativeSliceNumber_"):
                i_slice = int(param.name.split("_")[-1]) - 1
                order_idx_to_mrd_idx[int(param.value)] = i_slice
                mrd_idx_to_order_idx[i_slice] = int(param.value)
                nb_slices += 1
        logging.info(f"mrd_idx_to_order_idx: {mrd_idx_to_order_idx}")
        logging.info(f"order_idx_to_mrd_idx: {order_idx_to_mrd_idx}")

    except:
        raise ValueError("Improperly formatted metadata: \n%s", metadata)

    logging.info("Initialize nomad_f0xyz_opt")
    global is_sigint
    is_sigint = False
    global is_surrogate_ready
    is_surrogate_ready = False
    global is_obj_with_fmap_ready
    is_obj_with_fmap_ready = False
    global is_slice_loc_ready
    is_slice_loc_ready = False

    # Constants
    fatsat = False

    # Options
    max_evals = 2000
    logging.info(f"Objective: {objective}, use_surrogate: {use_surrogate}, channels_to_shim: {channels_to_shim}")

    nii_orig_mask = find_and_read_mask_within_data()
    nii_mask = None
    
    nomad_instance_stopped = [True for _ in range(nb_slices)]
    f_queues = []
    send_to_scanner_queue = Queue()
    send_to_scanner_thread = None
    result_queue = JoinableQueue()
    result_thread = None
    print_queue = JoinableQueue()
    print_queue_thread = None
    workers = []
    chronological_to_nifti_ordering = None

    nii_fmap = None
    nii_coil = None
    nii_orig_coil_profiles = None
    nii_coil_profiles = None

    # MI variables
    nii_anat_epi_space = None
    square_mask_coords = None

    if use_surrogate or is_obj_with_fmap:
        nii_fmap, json_data = find_and_read_fmap_within_data()
        isocenter = get_isocenter(json_data)
        st_scanner_constraints = {"name": "scanner",
                                  "coef_sum_max": sum_cstr,
                                  "coef_channel_minmax": {"0": [[f0_lb, f0_ub]],
                                                          "1": [[gradx_lb, gradx_ub],
                                                                [grady_lb, grady_ub],
                                                                [gradz_lb, gradz_ub]]}}
        scanner_coil = ScannerCoil(nii_fmap.shape[:3], nii_fmap.affine, constraints=st_scanner_constraints, orders=(0, 1),
                                   manufacturer="Siemens", isocenter=isocenter)
        nii_coil =nib.Nifti1Image(scanner_coil.profile, nii_fmap.affine, header=nii_fmap.header)
        nii_orig_coil_profiles = nii_coil
        nii_coil_profiles = nii_orig_coil_profiles
        # resample_from_to(nii_orig_coil_profiles, nii_fmap,  order=1, mode='grid-constant',  cval=0)

    if use_surrogate:
        logging.info("Using surrogate model")
        surr_queues = [Queue() for _ in range(nb_slices)]
    else:
        surr_queues = None

    if is_obj_with_fmap:
        logging.info("Using fmap obj function")
        obj_fmap_queues = [Queue() for _ in range(nb_slices)]
    else:
        obj_fmap_queues = None
    
    if use_f0_offset_from_gradients:
        slice_loc_queues = [Queue() for _ in range(nb_slices)]
    else:
        slice_loc_queues = None

    print_queue_thread = ReturnThread(target=process_print_queue, args=(print_queue,))
    print_queue_thread.start()
    result_thread = ReturnThread(target=process_worker_results, args=(result_queue,))
    result_thread.start()

    logging.info("Initializing f queues")
    # Initialize the queues that communicate between main thread and the optimization threads
    for i in range(nb_slices):
        f_queues.append(Queue())

    # Start the thread that communicates with the scanner
    logging.info("Starting thread to send coefs to Scanner")
    send_to_scanner_thread = Thread(target=send_coefs_to_scanner, args=(send_to_scanner_queue, print_queue, connection))
    send_to_scanner_thread.start()

    # Parse channels_to_shim
    if len(channels_to_shim) == 0:
        raise ValueError("channels_to_shim should not be empty")
    if len(channels_to_shim) > 4:
        raise ValueError("channels_to_shim should contain at most 4 characters")
    if len(channels_to_shim.replace("f", "").replace("x", "").replace("y", "").replace("z", "")) != 0:
        raise ValueError("channels_to_shim should only contain 'f', 'x', 'y' and 'z'")
    nb_channels_to_optimize = len(channels_to_shim)

    x0 = [0 for _ in range(nb_channels_to_optimize)]
    # Prisma has bounds of 2300 uT/m, the input needs to be in mT/m (we don't use min max values available since it's too large)
    lb = []
    ub = []
    if 'f' in channels_to_shim:
        lb.append(f0_lb)
        ub.append(f0_ub)
    if 'x' in channels_to_shim:
        lb.append(gradx_lb)
        ub.append(gradx_ub)
    if 'y' in channels_to_shim:
        lb.append(grady_lb)
        ub.append(grady_ub)
    if 'z' in channels_to_shim:
        lb.append(gradz_lb)
        ub.append(gradz_ub)

    # start the optimizer threads
    for i in range(nb_slices):
        logging.info(f"Starting nomad instance for slice {i}")
        nomad_instance_stopped[i] = False
        worker = Process(target=nomad_instance, args=(
            i, x0, lb, ub,
            nb_channels_to_optimize,
            max_evals,
            sum_cstr,
            f_queues[i],
            result_queue,
            send_to_scanner_queue,
            print_queue,
            use_surrogate,
            surr_queues[i] if surr_queues is not None else None,
            nii_fmap,
            nii_coil_profiles,
            is_obj_with_fmap,
            obj_fmap_queues[i] if obj_fmap_queues is not None else None,
            channels_to_shim,
            use_f0_offset_from_gradients,
            slice_loc_queues[i] if slice_loc_queues is not None else None)
        )
        worker.start()
        workers.append(worker)

    img_processed = 0
    # Continuously parse incoming data parsed from MRD messages
    currentSeries = 0
    nb_repetitions = 0
    # main thread will receive data from the scanner and give it to the right queues
    try:
        im_missing_to_complete_one_volume = nb_slices
        volume_images = []
        is_mask_converted = False

        for item in connection:
            if np.all(nomad_instance_stopped):
                raise TerminationException("All Nomad instances are stopped, exiting")
            if not is_processing:
                raise TerminationException("Stopped processing, exiting")
            # ----------------------------------------------------------
            # Raw k-space data messages
            # ----------------------------------------------------------
            if isinstance(item, ismrmrd.Acquisition):
                if item.is_flag_set(ismrmrd.ACQ_LAST_IN_SLICE):
                    print("Raw data sent")

            # ----------------------------------------------------------
            # Image data messages
            # ----------------------------------------------------------
            elif isinstance(item, ismrmrd.Image):
                if item.image_series_index != currentSeries:
                    logging.debug("Series of images sent")
                    currentSeries = item.image_series_index
                header = item.getHead()

                if (item.image_type is ismrmrd.IMTYPE_MAGNITUDE) or (item.image_type == 0):
                    logging.info(f"Image received: Slice: {header.slice}, Repetition: {header.repetition}")
                    # Make sure we receive the images in the right order
                    if header.repetition > (nb_repetitions + 1):
                        raise RuntimeError("Receiving images in wrong order")
                    nb_repetitions = header.repetition

                    # If we have all the images for a single volume, convert to NIfTI,
                    # otherwise, continue to the next item
                    if not is_mask_converted:
                        # Fill "volume_images" for a single volume
                        if im_missing_to_complete_one_volume > 0:
                            if header.repetition == 0 and header.phase == 0 and header.average == 0 and header.contrast == 0 and header.set == 0:
                                volume_images.append(item)
                                im_missing_to_complete_one_volume -= 1
                            else:
                                raise RuntimeError("Receiving images in wrong order")

                        if im_missing_to_complete_one_volume <= 0:
                            # Convert target from mrd 2 NIfTI
                            nii, json_data = mrd2nii_volume(metadata, volume_images)
                            nib.save(nii, f"{debugFolder}/nii_target.nii.gz")
                            with open(f"{debugFolder}/nii_target.json", 'w') as outfile:
                                json.dump(json_data, outfile, indent=2)

                            chronological_to_nifti_ordering = parse_slices(json_data)

                            # Mask
                            nii_mask = load_mask(nii, nii_orig_mask, path_output=debugFolder)
                            is_mask_converted = True

                            # Slice locations
                            if use_f0_offset_from_gradients:
                                slice_com_per_slice = []
                                iso = get_isocenter(json_data)
                                # Get slice locations
                                for i_slice in range(nb_slices):
                                    nifti_slice = chronological_to_nifti_ordering[i_slice][0]
                                    slice_com_per_slice.append(get_slice_mask_center_of_mass(nii_mask, nifti_slice))
                                    logging.debug(f"Slice center of mass for slice {i_slice} (nifti slice {nifti_slice}): {slice_com_per_slice[i_slice]}, isocenter: {iso}")
                                    slice_loc_queues[i_slice].put(slice_com_per_slice[i_slice])
                                    slice_loc_queues[i_slice].put(iso)
                                
                                np.save(os.path.join(debugFolder, "slice_com_per_slice.npy"), np.array(slice_com_per_slice))
                                np.save(os.path.join(debugFolder, "iso.npy"), np.array(iso))

                            # We need to convert the mask in the same geometry as the fmap to use the surrogate
                            if is_obj_with_fmap or use_surrogate:
                                # Mask needs to be resampled to the target
                                slices = [(i,) for i in range(nb_slices)]
                                results_mask = Parallel(-1, backend='loky')(
                                    delayed(resample_mask)(nii_mask, nii_fmap, slices[i],
                                                            dilation_kernel='sphere',
                                                            dilation_size=3)
                                    for i in range(nb_slices))
                                mask_per_slice = np.array([results_mask[it].get_fdata() for it in range(nb_slices)]).transpose(1, 2, 3, 0)
                                nib.save(nii_coil_profiles, f"{debugFolder}/nii_coil_profiles.nii.gz")
                                nib.save(nii_fmap, f"{debugFolder}/nii_fmap.nii.gz")
                                nii_mask_sliced = nib.Nifti1Image(mask_per_slice, nii_fmap.affine, header=nii_fmap.header)
                                nib.save(nii_mask_sliced, f"{debugFolder}/nii_sliced_mask.nii.gz")
                                time_ordering = parse_slices(json_data)
                                
                                for i_slice in range(nb_slices):
                                    # list_max_rmse needs to be initialized. It's time ordered
                                    nifti_slice = time_ordering[i_slice][0]
                                    res = nii_fmap.get_fdata()[mask_per_slice[..., nifti_slice] > 0]
                                    b0_rmse_coef = np.sqrt(np.mean(res ** 2))
                                    if b0_rmse_coef == 0:
                                        logging.warning("b0_rmse_coef is 0 for slice %d, setting it to 50", i_slice)
                                        b0_rmse_coef = 50
                                    list_max_rmse.append(MAX_RMSE_SCALAR * b0_rmse_coef)

                                if use_surrogate:
                                    # Send the mask when we have it so that the surrogate opt can start processing
                                    for i in range(nb_slices):
                                        surr_queues[i].put(mask_per_slice)
                                        surr_queues[i].put(time_ordering)
                                        surr_queues[i].put(list_max_rmse)
                                        if use_f0_offset_from_gradients:
                                            surr_queues[i].put(slice_com_per_slice[i])
                                            surr_queues[i].put(iso)
                                if is_obj_with_fmap:
                                    # Send the mask to the objective function
                                    for i in range(nb_slices):
                                        obj_fmap_queues[i].put(mask_per_slice)
                                        obj_fmap_queues[i].put(time_ordering)
                                        obj_fmap_queues[i].put(list_max_rmse)
                            
                            if is_obj_with_mi:
                                # Load anatomical image in EPI space
                                nii_anat_epi_space = load_anat_in_epi_space(nii)

                                # Pre-compute square mask coords
                                square_mask_coords = extract_rect_mask_coords(nii_mask.get_fdata())

                            for a_item in volume_images:
                                ret_code, nomad_instance_stopped = process_image(a_item, mrd_idx_to_order_idx, metadata, nomad_instance_stopped, is_processing, f_queues, result_queue, nii_mask,
                                                                                 is_obj_with_sigint, is_obj_with_mi, nii_anat_epi_space, square_mask_coords)
                                if ret_code:
                                    img_processed += 1
                        continue

                    ret_code, nomad_instance_stopped = process_image(item, mrd_idx_to_order_idx, metadata, nomad_instance_stopped, is_processing, f_queues, result_queue, nii_mask,
                                                                     is_obj_with_sigint, is_obj_with_mi, nii_anat_epi_space, square_mask_coords)

                    if ret_code:
                        img_processed += 1

            # ----------------------------------------------------------
            # Waveform data messages
            # ----------------------------------------------------------
            elif isinstance(item, ismrmrd.Waveform):
                # print("Waveform sent")
                pass

            elif item is None:
                break

            else:
                print("Unsupported data type %s", type(item).__name__)

    except Exception as e:
        logging.error(traceback.format_exc())
        connection.send_logging("ERROR   ", traceback.format_exc())

    finally:
        logging.info(f"Img processed: {img_processed}")
        stop(send_to_scanner_thread,
             nomad_instance_stopped,
             workers,
             result_queue,
             result_thread,
             print_queue,
             print_queue_thread,
             nb_slices,
             nb_repetitions,
             chronological_to_nifti_ordering,
             fatsat,
             send_to_scanner_queue,
             f_queues,
             channels_to_shim)
        connection.send_close()

def process_image(item, mrd_idx_to_order_idx, metadata, nomad_instance_stopped, is_processing, f_queues, result_queue, nii_mask,
                  is_obj_with_sigint, is_obj_with_mi, nii_anat_epi_space, square_mask_coords):
    # The slices are already in order of acquisition time
    slice_mrd = item.getHead().slice
    slice_nii = mrd_idx_to_order_idx[slice_mrd]

    if nomad_instance_stopped[slice_mrd]:
        logging.info("Nomad instance %d is stopped, not processing image", slice_mrd)
        return 0, nomad_instance_stopped

    if not is_processing:
        logging.info("Not processing image since 'is_processing' is false")
        return 0, nomad_instance_stopped

    currents_applied = True if item.user_int[5] == 1 else False
    # Todo: TEMP
    #currents_applied = True

    logging.info(f"B0 offset: {item.user_int[6] - (2 ** 16) if item.user_int[6] > 2 ** 15 else item.user_int[6]}")

    if currents_applied:
        mask = nii_mask.get_fdata()[..., slice_nii]

        if np.sum(np.abs(mask)) == 0:
            logging.info(f"Mask is empty, not processing image for slice: {slice_mrd}")
            # If nothing to shim in that slice, stop the instance
            f_queues[slice_mrd].put(None)
            # "optimal" currents should be 0s
            cur = Currents([0 for _ in range(NB_CHANNELS_TO_SEND_TO_SEQ)], slice_mrd, -np.inf, item.getHead().repetition)
            result_queue.put(cur)
            nomad_instance_stopped[slice_mrd] = True
            return 1, nomad_instance_stopped

        # Use mrd2nii to convert the received slice to NIfTI format
        nii_tmp = mrd2nii_stack(metadata, item, include_slice_gap=True)
        if item.getHead().repetition == 20:
            nib.save(nii_tmp, f"{debugFolder}/slice{slice_nii}_rep{item.getHead().repetition}.nii.gz")
            nib.Nifti1Image(mask, nii_tmp.affine, header=nii_tmp.header)
            nib.save(nib.Nifti1Image(mask, nii_tmp.affine, header=nii_tmp.header), f"{debugFolder}/mask_slice{slice_nii}.nii.gz")
        
        obj = calculate_obj_function(is_obj_with_sigint, is_obj_with_mi, nii_tmp, mask, nii_anat_epi_space, square_mask_coords, slice_nii)
        f_queues[slice_mrd].put((obj, item.getHead().repetition))
        logging.info(f"Slice {slice_mrd}, repetition {item.getHead().repetition}: objective: {obj}")
    else:
        logging.info(f"Image received for slice: {slice_mrd}, but no currents were applied")
        return 0, nomad_instance_stopped

    return 1, nomad_instance_stopped


def calculate_obj_function(use_sig_int, use_mi, nii_slice, mask, nii_anat_epi_space, square_mask_coords, slice_nii):
    if use_sig_int and not use_mi:
        obj = np.average(nii_slice.get_fdata(), weights=mask)
    elif not use_sig_int and use_mi:
        mutual_info = get_mututal_information_obj(nii_slice.get_fdata(), nii_anat_epi_space.get_fdata()[..., slice_nii], square_mask_coords[slice_nii], metric=mi_metric)
        obj = HYPERPARAM_MI * mutual_info
    elif use_sig_int and use_mi:
        avg_sig_int = np.average(nii_slice.get_fdata(), weights=mask)
        mutual_info = get_mututal_information_obj(nii_slice.get_fdata(), nii_anat_epi_space.get_fdata()[..., slice_nii], square_mask_coords[slice_nii], metric=mi_metric)
        obj = avg_sig_int + (HYPERPARAM_MI * mutual_info)
    else:
        raise RuntimeError("Unreachable. At least one of 'sig int' or 'mi' objective function should be selected")

    return obj


def nomad_instance(i_slice, x0, lb, ub, nb_channels, max_evals, sum_cstr, f_queue, result_queue, send_to_scanner_queue,
                   print_queue, use_surrogate, surr_queue, nii_fmap, nii_coil_profiles, use_fmap_obj, obj_fmap_queue, channels_to_shim, use_f0_offset_from_gradients,
                   slice_loc_queue):
    fname_cache = os.path.join(debugFolder, f"cache_slice{i_slice}.txt")
    fname_stats = os.path.join(debugFolder, f"stats_slice{i_slice}.txt")
    granularity = "("
    if "f" in channels_to_shim:
        granularity += "1 "
    if "x" in channels_to_shim:
        granularity += "0.001 "
    if "y" in channels_to_shim:
        granularity += "0.001 "
    if "z" in channels_to_shim:
        granularity += "0.001 "
    granularity = granularity[:-1] + ")"
    bb_params = [
        f"DIMENSION {nb_channels}",
        "BB_MAX_BLOCK_SIZE 3",
        "EVAL_OPPORTUNISTIC yes",
        f"MAX_EVAL {max_evals}",
        f"CACHE_FILE {fname_cache}",
        f"STATS_FILE {fname_stats}",
        f"GRANULARITY {granularity}",
        # "EVAL_USE_CACHE false"
        "LH_SEARCH 10 1"
    ]

    if use_fmap_obj:
        bb_params.append("BB_OUTPUT_TYPE OBJ EB EB")
    else:
        bb_params.append("BB_OUTPUT_TYPE OBJ EB")

    print_queue.put(f"Starting Nomad instance for slice {i_slice}, in process")
    if use_surrogate:
        # bb_params.append('VNS_MADS_SEARCH_WITH_SURROGATE true')
        bb_params.append('SURROGATE_MAX_BLOCK_SIZE 10')
        bb_params.append('EVAL_QUEUE_SORT SURROGATE')
        a = PyNomad.optimize(partial(bb_block,
                                 i_slice=i_slice,
                                 sum_cstr=sum_cstr,
                                 send_to_scanner_queue=send_to_scanner_queue,
                                 f_queue=f_queue,
                                 result_queue=result_queue,
                                 print_queue=print_queue,
                                 nii_fmap=nii_fmap,
                                 nii_coil_profiles=nii_coil_profiles,
                                 use_fmap_obj=use_fmap_obj,
                                 obj_fmap_queue=obj_fmap_queue,
                                 channels_to_shim=channels_to_shim,
                                 use_f0_offset_from_gradients=use_f0_offset_from_gradients,
                                 slice_loc_queue=slice_loc_queue),
                         x0, lb, ub, bb_params, partial(surrogate, i_slice=i_slice,
                                                        surr_queue=surr_queue,
                                                        nii_fmap=nii_fmap,
                                                        nii_coil_profiles=nii_coil_profiles,
                                                        print_queue=print_queue,
                                                        sum_cstr=sum_cstr,
                                                        use_fmap_obj=use_fmap_obj,
                                                        channels_to_shim=channels_to_shim,
                                                        use_f0_offset_from_gradients=use_f0_offset_from_gradients))
    else:

        a = PyNomad.optimize(partial(bb_block,
                                 i_slice=i_slice,
                                 sum_cstr=sum_cstr,
                                 send_to_scanner_queue=send_to_scanner_queue,
                                 f_queue=f_queue,
                                 result_queue=result_queue,
                                 print_queue=print_queue,
                                 nii_fmap=nii_fmap,
                                 nii_coil_profiles=nii_coil_profiles,
                                 use_fmap_obj=use_fmap_obj,
                                 obj_fmap_queue=obj_fmap_queue,
                                 channels_to_shim=channels_to_shim,
                                 use_f0_offset_from_gradients=use_f0_offset_from_gradients,
                                 slice_loc_queue=slice_loc_queue),
                         x0, lb, ub, bb_params)
    print_queue.put(f"Exited Nomad instance {i_slice}: {a}")


def bb_block(block, i_slice: int, sum_cstr, send_to_scanner_queue, f_queue, result_queue: Queue, print_queue, nii_fmap,
             nii_coil_profiles, use_fmap_obj, obj_fmap_queue, channels_to_shim, use_f0_offset_from_gradients, slice_loc_queue):
    nb_points = block.size()
    eval_ok = [False for _ in range(nb_points)]
    skips = [False for _ in range(nb_points)]

    global is_sigint
    if is_sigint:
        print_queue.put(f"Nomad instance {i_slice} interrupted")
        return eval_ok

    b0_rmse_list = []
    if use_fmap_obj:
        global is_obj_with_fmap_ready, obj_mask, obj_time_ordering_to_nifti_slices, list_max_rmse
        if not is_obj_with_fmap_ready:
            obj_mask = obj_fmap_queue.get()
            obj_time_ordering_to_nifti_slices = obj_fmap_queue.get()
            list_max_rmse = obj_fmap_queue.get()
            is_obj_with_fmap_ready = True
        a_slice = obj_time_ordering_to_nifti_slices[i_slice][0]
        print_queue.put(f"Objective function using fmap {i_slice}")

    if use_f0_offset_from_gradients:
        global is_slice_loc_ready, slice_loc, target_iso
        if not is_slice_loc_ready:
            slice_loc = slice_loc_queue.get()
            target_iso = slice_loc_queue.get()
            is_slice_loc_ready = True


    # Send the points to the scanner
    for i in range(nb_points):
        x = block.get_x(i)
        coefs = [x.get_coord(i) for i in range(x.size())]
        
        coefs_to_send_to_scanner_mt = fill_non_opt_currents_with_0s(coefs, channels_to_shim)
        if use_f0_offset_from_gradients:
            f0_offset = calculate_f0_from_gradients(coefs_to_send_to_scanner_mt[1:], target_iso, slice_loc)
            coefs_to_send_to_scanner_mt = [round(-f0_offset + coefs_to_send_to_scanner_mt[0])] + coefs_to_send_to_scanner_mt[1:]
        # Scale from mT/m to uT/m for the gradients
        coefs_to_send_to_scanner_ut = mult_grads_by_1000(coefs_to_send_to_scanner_mt)

        # Skip if it does not respect constraints
        skips[i] = np.sum(np.abs(coefs)) > sum_cstr
        if use_fmap_obj:
            res = nii_fmap.get_fdata()[obj_mask[..., a_slice] > 0] + nii_coil_profiles.get_fdata()[obj_mask[..., a_slice] > 0] @ coefs_to_send_to_scanner_ut
            b0_rmse_coef = np.sqrt(np.mean(res ** 2))
            b0_rmse_list.append(b0_rmse_coef)
            skips[i] = skips[i] or (b0_rmse_coef > list_max_rmse[i_slice])  # Example threshold, adjust as needed
        if skips[i]:
            continue
        currents = Currents(coefs_to_send_to_scanner_mt, i_slice)
        # print_queue.put(f"Sending currents to Scanner: {currents.currents}, Slice: {currents.slice}")
        send_to_scanner_queue.put(currents)

    for k in range(nb_points):
        x = block.get_x(k)
        coefs = [x.get_coord(i) for i in range(x.size())]
        cst = np.sum(np.abs(coefs)) - sum_cstr
        if skips[k]:
            # skipped because it did not respect the constraints
            ans = str(1) + " " + str(cst)
            if use_fmap_obj:
                ans += " " + str(b0_rmse_list[k] - list_max_rmse[i_slice])
            print_queue.put(f"Skipping since constraints not respected for slice {i_slice} with coefs: {coefs}, {ans=}")
            x.setBBO(ans.encode("UTF-8"))
            eval_ok[k] = 1
            continue

        a_tuple = f_queue.get()
        if a_tuple is None:
            # We put None in the queue to tell the nomad instance to stop
            is_sigint = True
            print_queue.put(f"Killing Nomad instance {i_slice}")
            os.kill(os.getpid(), signal.SIGINT)
            return eval_ok

        f = a_tuple[0]
        repetition = a_tuple[1]
        ans = str(-f) + " " + str(cst)
        if use_fmap_obj:
            ans += " " + str(b0_rmse_list[k] - list_max_rmse[i_slice])
        x.setBBO(ans.encode("UTF-8"))
        if result_queue.full():
            print_queue.put(f"Queue is full from nomad instance:{i_slice}")
        print_queue.put(f"Adding to result_queue for slice:{i_slice}")

        # Add 0s for channels that are not optimized
        coefs_to_send_to_scanner_mt = fill_non_opt_currents_with_0s(coefs, channels_to_shim)
        if use_f0_offset_from_gradients:
            f0_offset = calculate_f0_from_gradients(coefs_to_send_to_scanner_mt[1:], target_iso, slice_loc)
            coefs_to_send_to_scanner_mt = [round(-f0_offset + coefs_to_send_to_scanner_mt[0])] + coefs_to_send_to_scanner_mt[1:]

        result_queue.put(Currents(coefs_to_send_to_scanner_mt, i_slice, f, repetition))
        eval_ok[k] = 1

    return eval_ok


def surrogate(block, i_slice, surr_queue, nii_fmap, nii_coil_profiles, print_queue, sum_cstr, use_fmap_obj, channels_to_shim, use_f0_offset_from_gradients):

    global is_surrogate_ready, surr_mask, time_ordering_to_nifti_slices, list_max_rmse, slice_loc, target_iso
    if not is_surrogate_ready:
        surr_mask = surr_queue.get()
        time_ordering_to_nifti_slices = surr_queue.get()
        list_max_rmse = surr_queue.get()
        if use_f0_offset_from_gradients:
            slice_loc = surr_queue.get()
            target_iso = surr_queue.get()
        is_surrogate_ready = True

    a_slice = time_ordering_to_nifti_slices[i_slice][0]

    nb_points = block.size()
    eval_ok = [False for _ in range(nb_points)]
    for i in range(nb_points):
        x = block.get_x(i)
        coefs = [x.get_coord(i) for i in range(x.size())]
        coefs_to_send_to_scanner_mt = fill_non_opt_currents_with_0s(coefs, channels_to_shim)
        if use_f0_offset_from_gradients:
            f0_offset = calculate_f0_from_gradients(coefs_to_send_to_scanner_mt[1:], target_iso, slice_loc)
            coefs_to_send_to_scanner_mt = [-f0_offset + coefs_to_send_to_scanner_mt[0]] + coefs_to_send_to_scanner_mt[1:]
        # Scale from mT/m to uT/m for the gradients
        coefs_to_send_to_scanner_ut = mult_grads_by_1000(coefs_to_send_to_scanner_mt)
        # Constraints
        cst = np.sum(np.abs(coefs)) - sum_cstr

        res = nii_fmap.get_fdata()[surr_mask[..., a_slice] > 0] + nii_coil_profiles.get_fdata()[surr_mask[..., a_slice] > 0] @ coefs_to_send_to_scanner_ut
        b0_rmse_coef = np.sqrt(np.mean(res ** 2))
        print_queue.put(f"Surrogate: Evaluating point {i} with coefs: {coefs}, Slice: {i_slice}, B0 RMSE: {b0_rmse_coef}")
        ans = str(b0_rmse_coef) + " " + str(cst)
        if use_fmap_obj:
            ans += " " + str(b0_rmse_coef - list_max_rmse[i_slice])
        x.setBBO(ans.encode("UTF-8"))
        eval_ok[i] = True

    print_queue.put(f"Surrogate: {eval_ok}")
    return eval_ok


def fill_non_opt_currents_with_0s(currents, channels_to_shim):
    i = 0
    coefs_to_send_to_scanner = []
    channels_to_check = ['f', 'x', 'y', 'z']
    for ch in channels_to_check:
        if ch in channels_to_shim:
            coefs_to_send_to_scanner.append(currents[i])
            i += 1
        else:
            coefs_to_send_to_scanner.append(0.0)
    return coefs_to_send_to_scanner


def mult_grads_by_1000(a_list):
    if len(a_list) != 4:
        raise ValueError("Input list must have 4 elements")
    # Scale gradients, not f0
    return [a_list[0],] + [a_list[i] * 1000 for i in range(1, 4)]


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


def load_anat_in_epi_space(nii_epi):
    nii_anat = find_and_read_anat_within_data()
    if nii_anat is None:
        raise ValueError("Anatomical image not found in data folder")
    
    if nii_anat.ndim != 3:
        raise ValueError("Anatomical image must be 3D")
    elif not np.all(nii_anat.shape == nii_epi.shape) or not np.all(nii_anat.affine == nii_epi.affine):
        nii_anat_resampled = resample_from_to(nii_anat, nii_epi, order=2, mode='grid-constant', cval=0)
    else:
        logging.info("No need to resample the anatomical image")
        nii_anat_resampled = nii_anat

    return nii_anat_resampled


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


def get_isocenter(json_data):
    """ Get the isocenter location in RAS coordinates from the json file.

    The patient position is used to infer the table position in the patient coordinate system.
    When the table is at (0,0,0), the origin is at the isocenter. We can therefore infer
    the isocenter as -table_position when the table_position is in RAS coordinates.

    Args:
        json_data (dict): Dictionary containing the BIDS sidecar information

    Returns:
        numpy.ndarray: Isocenter location in RAS coordinates
    """
    table_position = json_data.get('TablePosition')

    patient_position = json_data.get('PatientPosition')

    if table_position is None or patient_position is None:
        raise ValueError("TablePosition or PatientPosition not found in json data")

    table_position = np.array(table_position)

    # Define coordinate transformations for each patient position
    position_transforms = {
        'HFS': [0, 1, 2],      # x=x, y=y, z=z
        'HFP': [0, 1, 2],      # x=-x, y=-y, z=z
        'FFS': [0, 1, 2],      # x=-x, y=y, z=-z
        'FFP': [0, 1, 2],      # x=x, y=-y, z=-z
        'LFP': [2, 1, 0],      # x=-z, y=-y, z=-x
        'LFS': [2, 1, 0],      # x=-z, y=y, z=x
        'RFP': [2, 1, 0],      # x=z, y=-y, z=x
        'RFS': [2, 1, 0],      # x=z, y=y, z=-x
        'HFDR': [1, 0, 2],     # x=-y, y=x, z=z
        'HFDL': [1, 0, 2],     # x=y, y=-x, z=z
        'FFDR': [1, 0, 2],     # x=-y, y=-x, z=-z
        'FFDL': [1, 0, 2],     # x=y, y=x, z=-z
        'AFDR': [1, 2, 0],     # x=-y, y=z, z=-x
        'AFDL': [1, 2, 0],     # x=y, y=z, z=x
        'PFDR': [1, 2, 0],     # x=-y, y=-z, z=x
        'PFDL': [1, 2, 0],     # x=y, y=-z, z=-x
    }

    # Define sign patterns for each patient position
    position_signs = {
        'HFS': [1, 1, 1],      'HFP': [-1, -1, 1],
        'FFS': [-1, 1, -1],    'FFP': [1, -1, -1],
        'LFP': [-1, -1, -1],   'LFS': [-1, 1, 1],
        'RFP': [1, -1, 1],     'RFS': [1, 1, -1],
        'HFDR': [-1, 1, 1],    'HFDL': [1, -1, 1],
        'FFDR': [-1, -1, -1],  'FFDL': [1, 1, -1],
        'AFDR': [-1, 1, -1],   'AFDL': [1, 1, 1],
        'PFDR': [-1, -1, 1],   'PFDL': [1, -1, -1],
    }

    if patient_position not in position_transforms:
        raise ValueError(f"Patient position {patient_position} not implemented")

    # Transform table position to RAS coordinates
    indices = position_transforms[patient_position]
    signs = position_signs[patient_position]

    table_position_ras = np.zeros(3)
    for i in range(3):
        table_position_ras[i] = signs[i] * table_position[indices[i]]

    # The isocenter is located at -table_position
    return -table_position_ras


def calculate_f0_from_gradients(grads, iso, location) -> float:
    """ Calculate the f0 offset induced by the gradients at the isocenter location.

    Args:
        grads (list): List of gradient strengths in mT/m [Gx, Gy, Gz]
        iso (list): Isocenter location in RAS coordinates [x, y, z] in mm
        location (list): Location where to calculate the f0 offset in RAS coordinates [x, y, z] in mm

    Returns:
        float: f0 offset in Hz
    """
    gamma = 42.577478518e6  # Hz/T
    f0_offset = 0

    # Convert mT/m to T/m and multiply by position difference in m
    x_offset = grads[0] * 1e-3 * ((location[0] - iso[0]) * 1e-3)
    y_offset = grads[1] * 1e-3 * ((location[1] - iso[1]) * 1e-3)
    z_offset = grads[2] * 1e-3 * ((location[2] - iso[2]) * 1e-3)

    # F0 offset depends on shim coordinate system, this is on Siemens, which is LAI:
    # x gradient increases field when going to the Left (negative x in RAS)
    # y gradient increases field when going to the Anterior (positive y in RAS)
    # z gradient increases field when going to the Inferior (negative z in RAS)
    f0_offset = -x_offset + y_offset - z_offset
    f0_offset_hz = f0_offset * gamma

    return f0_offset_hz


def get_slice_mask_center_of_mass(nii_mask, slice_number) -> list:
    """ Get the slice location in RAS coordinates from the nifti image and slice index.

    Args:
        nii_mask (nib.Nifti1Image): Nifti image
        slice_number (int): Slice index along the 3rd dimension
    Returns:
        list: Slice location in RAS coordinates [x, y, z] in meters
    """
    # Get the affine matrix
    affine = nii_mask.affine
    com = center_of_mass(nii_mask.get_fdata()[..., slice_number])
    voxel_coords = np.array([com[0], com[1], slice_number, 1])
    # Convert to RAS coordinates
    ras_coords = affine @ voxel_coords
    return ras_coords[:3].tolist()


def get_mututal_information_obj(epi_rect_mask, anat_epi_space_rect_mask, coords, metric='ssim'):
    xmin, xmax, ymin, ymax = coords
    epi_rect_mask = epi_rect_mask[xmin:xmax, ymin:ymax]
    anat_epi_space_rect_mask = anat_epi_space_rect_mask[xmin:xmax, ymin:ymax]
    if metric == 'ssim':
        mi = ssim(epi_rect_mask, anat_epi_space_rect_mask, data_range=anat_epi_space_rect_mask.max() - anat_epi_space_rect_mask.min(),
                gaussian_weights=True,
                sigma=0.5)
    elif metric == 'nmi':
        mi = (nmi(epi_rect_mask, anat_epi_space_rect_mask) - 1) * 1  # Normalize to be between 0 and 1
    else:
        raise ValueError(f"Metric {metric} not supported for mutual information objective function")
    
    return mi


def extract_rect_mask_coords(data_mask, extra_pad=1000):
    """
    Create a square mask based on the center of mass of data_mask for each slice
    and extract a 2D array from nii_anat_epi_space. The square size is inferred
    from the bounding box of the oval in nii_anat_epi_space.

    Args:
        data_mask (np.ndarray): 3D mask array (X, Y, Z).
        extra_pad (int): Extra padding to add to the square size.

    Returns:
        list: Extracted 2D arrays for each slice.
    """
    nz = data_mask.shape[2]
    extracted_coords = []

    for slice_idx in range(nz):
        # Compute the center of mass for the current slice
        com = center_of_mass(data_mask[:, :, slice_idx])
        if np.isnan(com[0]) or np.isnan(com[1]):  # Skip if the slice is empty
            raise ValueError(f"Slice {slice_idx} in data_mask is empty.")

        center_x, center_y = round(com[0]), round(com[1])

        # Infer square size from the bounding box of the oval
        non_zero_coords = np.argwhere(data_mask[:, :, slice_idx] > 0)  # Find non-zero (oval) coordinates
        if non_zero_coords.size == 0:  # Skip if no oval is found
            raise ValueError(f"No non-zero region found in slice {slice_idx} of data_mask.")

        x_min_oval, y_min_oval = non_zero_coords.min(axis=0)
        x_max_oval, y_max_oval = non_zero_coords.max(axis=0)
        rect_width = x_max_oval - x_min_oval
        rect_height = y_max_oval - y_min_oval

        half_width = rect_width // 2 + extra_pad
        half_height = rect_height // 2 + extra_pad

        # Define the rectangular region
        x_min = max(center_x - half_width, 0)
        x_max = min(center_x + half_width, data_mask.shape[0])
        y_min = max(center_y - half_height, 0)
        y_max = min(center_y + half_height, data_mask.shape[1])

        # append the coordinates
        extracted_coords.append((x_min, x_max, y_min, y_max))

    return extracted_coords



def create_obj_interactive_plot(fname_currents_per_volume, path_output, opt_channels):

    if not os.path.exists(path_output):
        os.makedirs(path_output)

    curr_per_volume = np.load(fname_currents_per_volume)  # shape (n_slices, n_iterations, n_channels + 1) with the last dimension being the objective value
    n_slices = curr_per_volume.shape[0]

    # Reorder according to the last column (objective value) 
    if opt_channels == ['f', 'z']:
        workers = []
        for i_slice in range(n_slices):
            worker = Process(target=plot_shim_path_wrapper, args=(i_slice, curr_per_volume, path_output))
            worker.start()
            workers.append(worker)
        
        for worker in workers:
            worker.join()
    else:
        print("Shim path plotting is only implemented for ['f', 'z'] optimization channels.")


def plot_shim_path_wrapper(slice_idx, curr_per_volume, path_output):

    signal_intensity = curr_per_volume[slice_idx, :, -1]
    fname_fig = os.path.join(path_output, f"slice_{slice_idx}_obj_through_time.png")
    plot_obj_through_time(signal_intensity, fname_fig)
    
    currents = curr_per_volume[slice_idx, :, :-1]
    plot_shim_path(currents, signal_intensity, os.path.join(path_output, f"slice_{slice_idx}_shim_path.gif"))


def plot_shim_path(currents: np.ndarray, obj_values: np.ndarray, fname_fig:str):

    idx_rem = []
    for idx in range(4):
        if np.allclose(np.nan_to_num(currents[:, idx]), 0):
            idx_rem.append(idx)
    currents_plot = np.delete(currents, idx_rem, axis=1)

    if currents_plot.shape[1] != 2:
        raise ValueError("Currents array must have exactly two current channels.")

    fig = plt.Figure(figsize=(5, 5))
    ax = fig.add_subplot(1, 1, 1)

    norm = plt.Normalize(np.nan_to_num(obj_values).min(), np.nan_to_num(obj_values).max())
    scat = ax.scatter(currents_plot[0, 0], currents_plot[0, 1], c=obj_values[1], norm=norm, cmap='jet')
    ax.set(xlim=[np.nan_to_num(currents_plot)[:, 0].min(), np.nan_to_num(currents_plot)[:, 0].max()], ylim=[np.nan_to_num(currents_plot)[:, 1].min(), np.nan_to_num(currents_plot)[:, 1].max()], xlabel='f', ylabel='z')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cax.set_title("Sig int.", fontsize=14)
    fig.colorbar(scat, cax=cax)

    def update(frame):
        # for each frame, update the data stored on each artist.
        x = currents_plot[:frame+1, 0]
        y = currents_plot[:frame+1, 1]
        # z = obj_values[frame]
        # update the scatter plot:
        data = np.stack([x, y]).T
        scat.set_offsets(data)
        scat.set_array(obj_values[:frame+1])
        # scat.set_clim(vmin=min(obj_values), vmax=max(obj_values))
        return (scat,)

    ani = animation.FuncAnimation(fig=fig, func=update, frames=currents.shape[0], interval=30)
    ani.save(filename=fname_fig, writer='pillow')
    #plt.close('all')


def plot_obj_through_time(obj_values: np.ndarray, fname_fig):
    """ Plot objective function values through time.

    Args:
        obj_values (np.ndarray): Array of objective function values.
    """
    fig = plt.Figure(figsize=(5, 5))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(obj_values, marker='o')
    ax.set_title("Objective Function Values Through Time")
    ax.set_xlabel("Time Point")
    ax.set_ylabel("Objective Function Value")
    ax.grid()
    fig.savefig(fname_fig)


class TerminationException(Exception):
    """Exception to be raised when the optimization is terminated by the user."""
    pass
