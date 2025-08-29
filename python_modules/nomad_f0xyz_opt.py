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
from scipy.ndimage import binary_dilation, binary_erosion, generate_binary_structure, iterate_structure
import signal
import shutil
from threading import Thread
import traceback

from mrd2nii.mrd2nii_main import mrd2nii_volume, mrd2nii_stack

# Folder for debug output files
debugFolder = os.path.abspath("/tmp/share/debug/nomad_files")
dataFolder = os.path.abspath("/tmp/share/saved_data")
fname_solution = os.path.join(debugFolder, "solution.txt")

NB_CHANNELS = 4

# Signal interrupt for bb_block
is_sigint = False

# Surrogate global variables
is_surrogate_ready = False
surr_mask = None  # 4d fmap coord (X, Y, Z, ishim)
time_ordering_to_nifti_slices = None

# obj with fmap
is_obj_with_fmap_ready = None
obj_mask = None
obj_time_ordering_to_nifti_slices = None
MAX_RMSE = 100


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
        ('currents',  ctypes.c_float * NB_CHANNELS) #          NB_CHANNELS bytes
    ]                                               #       =  NB_CHANNELS + 1 bytes total


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


def find_and_read_mask_within_data():
    # For now, take the latest created mask in dataFolder
    mask_files = [f for f in os.listdir(dataFolder) if "mask" in f and f.endswith('.nii.gz')]
    if not mask_files:
        raise RuntimeError(f"No mask files found. Is it in {dataFolder}?")

    latest_mask_file = max(mask_files, key=lambda f: os.path.getctime(os.path.join(dataFolder, f)))
    logging.info(f"Using latest mask file: {latest_mask_file}")
    return nib.load(os.path.join(dataFolder, latest_mask_file))


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
            print_queue.put(f"Sending currents to Teensy: {currents.currents}, Slice: {currents.slice}")
            feedback_data = MyFeedbackData()
            feedback_data.slice_nb = currents.slice
            for i in range(NB_CHANNELS):
                feedback_data.currents[i] = currents.currents[i]
            connection.send_feedback("MyFeedback", feedback_data)
        else:
            print_queue.put(f"Number of currents sent to Teensy: {n_currents}")
            break

def stop(send_to_scanner_thread, nomad_instance_stopped, workers, result_queue, result_thread,
         print_queue, print_queue_thread, nb_slices, nb_repetitions, chronological_to_nifti_ordering, fatsat, send_to_scanner_queue, f_queues):
    
    def process_results(currents, nb_slices, nb_repetitions, chronological_to_nifti_ordering):

        logging.info("Processing results")
        opt_currents = np.zeros((nb_slices, NB_CHANNELS))
        max_sig_int = np.zeros((nb_slices,))
        logging.info(f"Number of processed currents received in result_queue: {len(currents)}")

        all_currents = np.zeros((nb_slices, nb_repetitions + 1, NB_CHANNELS + 1))

        for current in currents:
            if current.f > max_sig_int[current.slice]:
                max_sig_int[current.slice] = current.f
                opt_currents[current.slice, :] = current.currents
            all_currents[current.slice, current.repetition, :-1] = current.currents
            all_currents[current.slice, current.repetition, -1] = current.f

        currents_per_volume = np.empty(shape=(nb_slices, nb_repetitions + 1, NB_CHANNELS + 1))
        currents_per_volume[:] = np.nan

        for current in currents:
            a_cur = copy.deepcopy(current.currents)
            a_cur.append(current.f)
            nifti_slice = chronological_to_nifti_ordering[current.slice][0]
            currents_per_volume[nifti_slice, current.repetition, :] = a_cur

        logging.info("Processed currents")
        np.save(os.path.join(debugFolder, f"all_currents.npy"), all_currents)
        np.save(os.path.join(debugFolder, f"currents_per_volume.npy"), currents_per_volume)
        return opt_currents
    
    def write_text_file(data, nb_slices, fatsat):
        logging.info("Writing currents to text file")
        with open(fname_solution, 'w') as f:
            for i_slice in range(nb_slices):
                if fatsat:
                    for i_channel in range(NB_CHANNELS):
                        f.write(f"0.0, ")
                    f.write("\n")
                for i_channel in range(NB_CHANNELS):
                    f.write(f"{data[i_slice, i_channel]}, ")
                f.write("\n")
        logging.info(f"Currents written to {fname_solution}")
    
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

    # stop worker threads
    logging.info("Stopping worker thread")
    for worker in workers:
        worker.join()

    result_queue.put(None)
    result_currents = result_thread.join()
    opt_currents = process_results(result_currents, nb_slices, nb_repetitions, chronological_to_nifti_ordering)

    # Write the currents to a text file
    write_text_file(opt_currents, nb_slices, fatsat)
    logging.info("Processing print queue")
    print_queue.put(None)
    messages = print_queue_thread.join()
    for message in messages:
        logging.info(message)
    print_queue.join()
    logging.info(opt_currents)


def process(connection, config, metadata):
    if os.path.exists(debugFolder):
        shutil.rmtree(debugFolder)
    os.makedirs(debugFolder)
    logging.info("Optimizing")
    is_processing = True
    
    logging.info(f"Config: {config}\n")

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

    # Constants
    fatsat = False

    # Options
    max_evals = 500
    use_surrogate = False
    objective = "Sig int"


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

    use_fmap_obj = True if objective == "Sig int + fmap" else False
    nii_fmap = None
    nii_coil = None
    nii_orig_coil_profiles = None
    nii_coil_profiles = None

    if use_surrogate or use_fmap_obj:
        nii_fmap = nii_fmap
        nii_orig_coil_profiles = nii_coil
        nii_coil_profiles = resample_from_to(nii_orig_coil_profiles,
                                             nii_fmap,
                                             order=1,
                                             mode='grid-constant',
                                             cval=0)

    if use_surrogate:
        logging.info("Using surrogate model")
        surr_queues = [Queue() for _ in range(nb_slices)]
    else:
        surr_queues = None

    if use_fmap_obj:
        logging.info("Using fmap obj function")
        obj_fmap_queues = [Queue() for _ in range(nb_slices)]
    else:
        obj_fmap_queues = None

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

    x0 = [0 for _ in range(NB_CHANNELS)]
    # Prisma has 2300 uT/m, the input needs to be in mT/m
    # lb = [-1000, -2.3, -2.3, -2.3]
    lb = [-1000, -1, -1, -1]
    ub = [1000, 1, 1, 1]
    sum_cstr = np.inf

    # start the optimizer threads
    for i in range(nb_slices):
        logging.info(f"Starting nomad instance for slice {i}")
        nomad_instance_stopped[i] = False
        worker = Process(target=nomad_instance, args=(
            i, x0, lb, ub,
            NB_CHANNELS,
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
            use_fmap_obj,
            obj_fmap_queues[i] if obj_fmap_queues is not None else None))
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

                            # We need to convert the mask in the same geometry as the fmap to use the surrogate
                            if use_fmap_obj or use_surrogate:
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
                                nii_mask_sliced = nib.Nifti1Image(mask_per_slice, nii_fmap.affine,
                                                                    header=nii_fmap.header)
                                nib.save(nii_mask_sliced, f"{debugFolder}/nii_sliced_mask.nii.gz")
                                time_ordering = parse_slices(json_data)

                                if use_surrogate:
                                    # Send the mask when we have it so that the surrogate opt can start processing
                                    for i in range(nb_slices):
                                        surr_queues[i].put(mask_per_slice)
                                        surr_queues[i].put(time_ordering)
                                if use_fmap_obj:
                                    # Send the mask to the objective function
                                    for i in range(nb_slices):
                                        obj_fmap_queues[i].put(mask_per_slice)
                                        obj_fmap_queues[i].put(time_ordering)

                            for a_item in volume_images:
                                ret_code, nomad_instance_stopped = process_image(a_item, mrd_idx_to_order_idx, metadata, nomad_instance_stopped, is_processing, f_queues, result_queue, nii_mask)
                                if ret_code:
                                    img_processed += 1
                        continue

                    ret_code, nomad_instance_stopped = process_image(item, mrd_idx_to_order_idx, metadata, nomad_instance_stopped, is_processing, f_queues, result_queue, nii_mask)

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
        connection.send_close()
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
             f_queues)

def process_image(item, mrd_idx_to_order_idx, metadata, nomad_instance_stopped, is_processing, f_queues, result_queue, nii_mask):
    # The slices are already in order of acquisition time
    slice_mrd = item.getHead().slice

    if nomad_instance_stopped[slice_mrd]:
        logging.info("Nomad instance %d is stopped, not processing image", slice_mrd)
        return 0, nomad_instance_stopped

    if not is_processing:
        logging.info("Not processing image since 'is_processing' is false")
        return 0, nomad_instance_stopped

    currents_applied = True if item.user_int[5] == 1 else False
    logging.info(f"user_int: {item.user_int[:]}")
    logging.info(f"user_float: {item.user_float[:]}")

    if currents_applied:
        mask = nii_mask.get_fdata()[..., mrd_idx_to_order_idx[slice_mrd]]

        if np.sum(np.abs(mask)) == 0:
            logging.info(f"Mask is empty, not processing image for slice: {slice_mrd}")
            # If nothing to shim in that slice, stop the instance
            f_queues[slice_mrd].put(None)
            # "optimal" currents should be 0s
            cur = Currents([0 for _ in range(NB_CHANNELS)], slice_mrd, -np.inf, item.getHead().repetition)
            result_queue.put(cur)
            nomad_instance_stopped[slice_mrd] = True
            return 1, nomad_instance_stopped

        # Use mrd2nii
        nii_tmp = mrd2nii_stack(metadata, item, include_slice_gap=True)
        if item.getHead().repetition == 20:
            nib.save(nii_tmp, f"{debugFolder}/slice{mrd_idx_to_order_idx[slice_mrd]}_rep{item.getHead().repetition}.nii.gz")
            nib.Nifti1Image(mask, nii_tmp.affine, header=nii_tmp.header)
            nib.save(nib.Nifti1Image(mask, nii_tmp.affine, header=nii_tmp.header), f"{debugFolder}/mask_slice{mrd_idx_to_order_idx[slice_mrd]}.nii.gz")
        avg_sig_int = np.average(nii_tmp.get_fdata(), weights=mask)
        f_queues[slice_mrd].put((avg_sig_int, item.getHead().repetition))
    else:
        logging.info(f"Image received for slice: {slice_mrd}, but no currents were applied")
        return 0, nomad_instance_stopped

    return 1, nomad_instance_stopped


def nomad_instance(i_slice, x0, lb, ub, nb_channels, max_evals, sum_cstr, f_queue, result_queue, send_to_scanner_queue,
                   print_queue, use_surrogate, surr_queue, nii_fmap, nii_coil_profiles, use_fmap_obj, obj_fmap_queue):
    fname_cache = os.path.join(debugFolder, f"cache_slice{i_slice}.txt")
    fname_stats = os.path.join(debugFolder, f"stats_slice{i_slice}.txt")
    bb_params = [
        f"DIMENSION {nb_channels}",
        # "BB_OUTPUT_TYPE OBJ EB"
        "BB_MAX_BLOCK_SIZE 3",
        "EVAL_OPPORTUNISTIC yes",
        f"MAX_EVAL {max_evals}",
        f"CACHE_FILE {fname_cache}",
        f"STATS_FILE {fname_stats}",
        "GRANULARITY * 0.001",
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
        PyNomad.optimize(partial(bb_block,
                                 i_slice=i_slice,
                                 sum_cstr=sum_cstr,
                                 send_to_scanner_queue=send_to_scanner_queue,
                                 f_queue=f_queue,
                                 result_queue=result_queue,
                                 print_queue=print_queue,
                                 nii_fmap=nii_fmap,
                                 nii_coil_profiles=nii_coil_profiles,
                                 use_fmap_obj=use_fmap_obj,
                                 obj_fmap_queue=obj_fmap_queue),
                         x0, lb, ub, bb_params, partial(surrogate, i_slice=i_slice,
                                                        surr_queue=surr_queue,
                                                        nii_fmap=nii_fmap,
                                                        nii_coil_profiles=nii_coil_profiles,
                                                        print_queue=print_queue,
                                                        sum_cstr=sum_cstr,
                                                        use_fmap_obj=use_fmap_obj))
    else:

        PyNomad.optimize(partial(bb_block,
                                 i_slice=i_slice,
                                 sum_cstr=sum_cstr,
                                 send_to_scanner_queue=send_to_scanner_queue,
                                 f_queue=f_queue,
                                 result_queue=result_queue,
                                 print_queue=print_queue,
                                 nii_fmap=nii_fmap,
                                 nii_coil_profiles=nii_coil_profiles,
                                 use_fmap_obj=use_fmap_obj,
                                 obj_fmap_queue=obj_fmap_queue),
                         x0, lb, ub, bb_params)
    print_queue.put(f"Exited Nomad instance {i_slice}")


def bb_block(block, i_slice: int, sum_cstr, send_to_scanner_queue, f_queue, result_queue: Queue, print_queue, nii_fmap,
             nii_coil_profiles, use_fmap_obj, obj_fmap_queue):

    nb_points = block.size()
    eval_ok = [False for _ in range(nb_points)]
    skips = [False for _ in range(nb_points)]

    global is_sigint
    if is_sigint:
        print_queue.put(f"Nomad instance {i_slice} interrupted")
        return eval_ok

    b0_rmse_list = []
    if use_fmap_obj:
        global is_obj_with_fmap_ready, obj_mask, obj_time_ordering_to_nifti_slices
        if not is_obj_with_fmap_ready:
            obj_mask = obj_fmap_queue.get()
            obj_time_ordering_to_nifti_slices = obj_fmap_queue.get()
            is_obj_with_fmap_ready = True
        a_slice = obj_time_ordering_to_nifti_slices[i_slice][0]
        print_queue.put(f"Objective function using fmap {i_slice}")

    # Send the points to the scanner
    for i in range(nb_points):
        x = block.get_x(i)
        coefs = [x.get_coord(i) for i in range(x.size())]

        # Skip if it does not respect constraints
        skips[i] = np.sum(np.abs(coefs)) > sum_cstr
        if use_fmap_obj:
            res = nii_fmap.get_fdata()[obj_mask[..., a_slice] > 0] + nii_coil_profiles.get_fdata()[obj_mask[..., a_slice] > 0] @ coefs
            b0_rmse_coef = np.sqrt(np.mean(res ** 2))
            b0_rmse_list.append(b0_rmse_coef)
            skips[i] = skips[i] or (b0_rmse_coef > MAX_RMSE)  # Example threshold, adjust as needed
        if skips[i]:
            continue
        currents = Currents(coefs, i_slice)
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
                ans += " " + str(b0_rmse_list[k] - MAX_RMSE)
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
            ans += " " + str(b0_rmse_list[k] - MAX_RMSE)
        x.setBBO(ans.encode("UTF-8"))
        # print_queue.put(f"result_queue.maxsize: {result_queue._maxsize}")
        if result_queue.full():
            print_queue.put(f"Queue is full from nomad instance:{i_slice}")
        print_queue.put(f"Adding to result_queue for slice:{i_slice}")
        result_queue.put(Currents(coefs, i_slice, f, repetition))
        eval_ok[k] = 1

    return eval_ok


def surrogate(block, i_slice, surr_queue, nii_fmap, nii_coil_profiles, print_queue, sum_cstr, use_fmap_obj):

    global is_surrogate_ready, surr_mask, time_ordering_to_nifti_slices
    if not is_surrogate_ready:
        surr_mask = surr_queue.get()
        time_ordering_to_nifti_slices = surr_queue.get()
        is_surrogate_ready = True

    a_slice = time_ordering_to_nifti_slices[i_slice][0]

    nb_points = block.size()
    eval_ok = [False for _ in range(nb_points)]
    for i in range(nb_points):
        x = block.get_x(i)
        coefs = [x.get_coord(i) for i in range(x.size())]

        # Constraints
        cst = np.sum(np.abs(coefs)) - sum_cstr

        res = nii_fmap.get_fdata()[surr_mask[..., a_slice] > 0] + nii_coil_profiles.get_fdata()[surr_mask[..., a_slice] > 0] @ coefs
        b0_rmse_coef = np.sqrt(np.mean(res ** 2))
        print_queue.put(f"Surrogate: Evaluating point {i} with coefs: {coefs}, Slice: {i_slice}, B0 RMSE: {b0_rmse_coef}")
        ans = str(b0_rmse_coef) + " " + str(cst)
        if use_fmap_obj:
            ans += " " + str(b0_rmse_coef - MAX_RMSE)
        x.setBBO(ans.encode("UTF-8"))
        eval_ok[i] = True

    print_queue.put(f"Surrogate: {eval_ok}")
    return eval_ok



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


class TerminationException(Exception):
    """Exception to be raised when the optimization is terminated by the user."""
    pass
