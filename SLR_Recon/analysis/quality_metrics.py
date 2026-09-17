import argparse
import csv
import os
import shutil
import subprocess
import sys
import nibabel as nib
import numpy as np


def run_cmd(cmd_list, description=""):
    print(f"\n[SCT] {description}...")
    try:
        subprocess.run(cmd_list, check=True, text=True)
    except subprocess.CalledProcessError:
        print(f"\n[ERROR] Command failed: {' '.join(cmd_list)}")
        sys.exit(1)


def load_volume_data(image_path, echo_idx=0):
    img = nib.load(image_path)
    data = np.nan_to_num(np.abs(img.get_fdata()))

    if data.ndim == 4:
        print(f"--> Input is 4D {data.shape}, extracting Echo index {echo_idx}")
        data = data[:, :, :, echo_idx]
    elif data.ndim != 3:
        raise ValueError(f"Expected 3D or 4D volume, got shape {data.shape}")

    return data, img.affine, img.header


def compute_metrics_and_write_csv(image_path, cord_seg_path, gm_seg_path, csf_seg_path, csv_out_path, echo_idx=0, thresh=0.5):
    img_data, _, _ = load_volume_data(image_path, echo_idx=echo_idx)
    
    cord_mask = nib.load(cord_seg_path).get_fdata() > thresh
    
    has_gm = gm_seg_path is not None and os.path.exists(gm_seg_path)
    if has_gm:
        gm_mask = nib.load(gm_seg_path).get_fdata() > thresh
    else:
        print("--> Warning: GM segmentation file not found. GM metrics will be skipped.")
        gm_mask = np.zeros_like(cord_mask, dtype=bool)

    wm_mask = cord_mask & (~gm_mask)
    nib.save(nib.Nifti1Image(wm_mask.astype(np.uint8), np.eye(4)), os.path.join(os.path.dirname(csv_out_path), "wm_seg.nii.gz"))
    has_csf = csf_seg_path is not None and os.path.exists(csf_seg_path)
    if has_csf:
        csf_mask = nib.load(csf_seg_path).get_fdata() > thresh
    else:
        print("Warning: CSF segmentation file not found. CSF metrics will be skipped.")
        csf_mask = np.zeros_like(cord_mask, dtype=bool)

    nz = img_data.shape[2]
    slice_rows = []

    gm_snrs, wm_snrs, csf_snrs = [], [], []
    cnrs_gm_wm, cnrs_wm_csf = [], []

    for z in range(nz):
        slice_img = img_data[:, :, z]
        s_gm_mask = gm_mask[:, :, z]
        s_wm_mask = wm_mask[:, :, z]
        s_csf_mask = csf_mask[:, :, z]

        v_gm = slice_img[s_gm_mask]
        v_wm = slice_img[s_wm_mask]
        v_csf = slice_img[s_csf_mask]

        # GM Stats
        mean_gm = np.mean(v_gm) if len(v_gm) > 0 else 0.0
        std_gm = np.std(v_gm) if len(v_gm) > 0 else 0.0
        snr_gm = mean_gm / std_gm if std_gm > 0 else 0.0

        # WM Stats
        mean_wm = np.mean(v_wm) if len(v_wm) > 0 else 0.0
        std_wm = np.std(v_wm) if len(v_wm) > 0 else 0.0
        snr_wm = mean_wm / std_wm if std_wm > 0 else 0.0

        # CSF Stats
        mean_csf = np.mean(v_csf) if len(v_csf) > 0 else 0.0
        std_csf = np.std(v_csf) if len(v_csf) > 0 else 0.0
        snr_csf = mean_csf / std_csf if std_csf > 0 else 0.0

        # GM / WM CNR
        denom_gm_wm = np.sqrt(std_gm**2 + std_wm**2)
        cnr_gm_wm = abs(mean_gm - mean_wm) / denom_gm_wm if denom_gm_wm > 0 else 0.0

        # WM / CSF CNR
        denom_wm_csf = np.sqrt(std_wm**2 + std_csf**2)
        cnr_wm_csf = abs(mean_wm - mean_csf) / denom_wm_csf if denom_wm_csf > 0 else 0.0

        slice_rows.append({
            "Slice": f"Slice_{z:02d}",
            "GM_Mean": round(float(mean_gm), 3),
            "GM_Std": round(float(std_gm), 3),
            "GM_SNR": round(float(snr_gm), 3),
            "WM_Mean": round(float(mean_wm), 3),
            "WM_Std": round(float(std_wm), 3),
            "WM_SNR": round(float(snr_wm), 3),
            "CSF_Mean": round(float(mean_csf), 3),
            "CSF_Std": round(float(std_csf), 3),
            "CSF_SNR": round(float(snr_csf), 3),
            "CNR_GM_WM": round(float(cnr_gm_wm), 3),
            "CNR_WM_CSF": round(float(cnr_wm_csf), 3)
        })

        if snr_gm > 0: gm_snrs.append(snr_gm)
        if snr_wm > 0: wm_snrs.append(snr_wm)
        if snr_csf > 0: csf_snrs.append(snr_csf)
        if cnr_gm_wm > 0: cnrs_gm_wm.append(cnr_gm_wm)
        if cnr_wm_csf > 0: cnrs_wm_csf.append(cnr_wm_csf)

    v_gm_3d = img_data[gm_mask]
    v_wm_3d = img_data[wm_mask]
    v_csf_3d = img_data[csf_mask]
    
    snr_gm_3d = np.mean(v_gm_3d) / np.std(v_gm_3d) if len(v_gm_3d) > 0 and np.std(v_gm_3d) > 0 else 0.0
    snr_wm_3d = np.mean(v_wm_3d) / np.std(v_wm_3d) if len(v_wm_3d) > 0 and np.std(v_wm_3d) > 0 else 0.0
    snr_csf_3d = np.mean(v_csf_3d) / np.std(v_csf_3d) if len(v_csf_3d) > 0 and np.std(v_csf_3d) > 0 else 0.0

    cnr_gm_wm_3d = abs(np.mean(v_gm_3d) - np.mean(v_wm_3d)) / np.sqrt(np.std(v_gm_3d)**2 + np.std(v_wm_3d)**2) if len(v_gm_3d) > 0 and len(v_wm_3d) > 0 else 0.0
    cnr_wm_csf_3d = abs(np.mean(v_wm_3d) - np.mean(v_csf_3d)) / np.sqrt(np.std(v_wm_3d)**2 + np.std(v_csf_3d)**2) if len(v_wm_3d) > 0 and len(v_csf_3d) > 0 else 0.0

    headers = [
        "Slice", 
        "GM_Mean", "GM_Std", "GM_SNR", 
        "WM_Mean", "WM_Std", "WM_SNR", 
        "CSF_Mean", "CSF_Std", "CSF_SNR", 
        "CNR_GM_WM", "CNR_WM_CSF"
    ]
    
    with open(csv_out_path, mode="w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(slice_rows)
        
        # Summary rows
        writer.writerow({
            "Slice": "MEAN_ACROSS_SLICES",
            "GM_SNR": round(float(np.mean(gm_snrs)), 3) if gm_snrs else 0,
            "WM_SNR": round(float(np.mean(wm_snrs)), 3) if wm_snrs else 0,
            "CSF_SNR": round(float(np.mean(csf_snrs)), 3) if csf_snrs else 0,
            "CNR_GM_WM": round(float(np.mean(cnrs_gm_wm)), 3) if cnrs_gm_wm else 0,
            "CNR_WM_CSF": round(float(np.mean(cnrs_wm_csf)), 3) if cnrs_wm_csf else 0
        })

        writer.writerow({
            "Slice": "STD_ACROSS_SLICES",
            "GM_SNR": round(float(np.std(gm_snrs)), 3) if gm_snrs else 0,
            "WM_SNR": round(float(np.std(wm_snrs)), 3) if wm_snrs else 0,
            "CSF_SNR": round(float(np.std(csf_snrs)), 3) if csf_snrs else 0,
            "CNR_GM_WM": round(float(np.std(cnrs_gm_wm)), 3) if cnrs_gm_wm else 0,
            "CNR_WM_CSF": round(float(np.std(cnrs_wm_csf)), 3) if cnrs_wm_csf else 0
        })
        writer.writerow({
            "Slice": "OVERALL_3D_VOLUME",
            "GM_SNR": round(float(snr_gm_3d), 3),
            "WM_SNR": round(float(snr_wm_3d), 3),
            "CSF_SNR": round(float(snr_csf_3d), 3),
            "CNR_GM_WM": round(float(cnr_gm_wm_3d), 3),
            "CNR_WM_CSF": round(float(cnr_wm_csf_3d), 3)
        })

    print(f"\n--> Successfully saved CSV metrics report: {csv_out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Automated Spinal Cord MRI Pipeline: Takes a 3D/4D NIfTI volume, runs SCT segmentations, and computes slice-wise metrics."
    )
    parser.add_argument("-i", "--input", required=True, help="Path to input 3D or 4D NIfTI file (e.g. stacked_e1.nii.gz)")
    parser.add_argument("-o", "--output_dir", required=True, help="Folder where results and CSV will be saved")
    parser.add_argument("--echo_idx", type=int, default=0, help="Echo index if input volume is 4D (default: 0)")
    parser.add_argument("--seg", default=None, help="Optional pre-made whole cord segmentation (.nii.gz)")
    parser.add_argument("--gm_seg", default=None, help="Optional pre-made grey matter segmentation (.nii.gz)")
    parser.add_argument("--csf_seg", default=None, help="Optional pre-made CSF segmentation (.nii.gz)")
    parser.add_argument("--csv_name", default="metrics_summary.csv", help="Output CSV filename (default: metrics_summary.csv)")
    
    args = parser.parse_args()

    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file does not exist: {args.input}")

    os.makedirs(args.output_dir, exist_ok=True)

    input_volume_path = args.input

    seg_path = os.path.join(args.output_dir, "sc_seg.nii.gz")
    if args.seg and os.path.exists(args.seg):
        print(f"--> Using provided whole cord segmentation: {args.seg}")
        shutil.copy(args.seg, seg_path)
    else:
        run_cmd(
            ["sct_deepseg", "spinalcord", "-i", input_volume_path, "-o", seg_path],
            description="Running SCT Whole Spinal Cord Deep Segmentation"
        )

    gm_seg_path = os.path.join(args.output_dir, "gm_seg.nii.gz")
    if args.gm_seg and os.path.exists(args.gm_seg):
        print(f"--> Using provided GM segmentation: {args.gm_seg}")
        shutil.copy(args.gm_seg, gm_seg_path)
    else:
        run_cmd(
            ["sct_deepseg", "graymatter", "-i", input_volume_path, "-o", gm_seg_path],
            description="Running SCT Grey Matter Deep Segmentation"
        )

    csf_seg_path = None
    if args.csf_seg and os.path.exists(args.csf_seg):
        print(f"--> Using provided CSF segmentation: {args.csf_seg}")
        csf_seg_path = os.path.join(args.output_dir, "csf_seg.nii.gz")
        shutil.copy(args.csf_seg, csf_seg_path)

    csv_out_path = os.path.join(args.output_dir, args.csv_name)
    compute_metrics_and_write_csv(
        image_path=input_volume_path,
        cord_seg_path=seg_path,
        gm_seg_path=gm_seg_path,
        csf_seg_path=csf_seg_path,
        csv_out_path=csv_out_path,
        echo_idx=args.echo_idx
    )


if __name__ == "__main__":
    main()