import argparse
import glob
import os
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns


def extract_metadata_from_filename(filename):
    """Extracts Subject ID and Anatomical Region from filename."""
    fname_lower = filename.lower()
    
    region = "Unknown"
    for r in ["cervical", "thoracic", "lumbar"]:
        if r in fname_lower:
            region = r.capitalize()
            break

    sub_match = re.search(r"(sub-?\d+)", fname_lower)
    subject = sub_match.group(1) if sub_match else "Unknown"

    return subject, region


def load_technique_directory(folder_path, technique_name):
    """Loads all CSVs from a directory and tags them with technique, subject, and region."""
    if not os.path.exists(folder_path):
        raise FileNotFoundError(f"Directory not found: {folder_path}")

    csv_files = glob.glob(os.path.join(folder_path, "*.csv"))
    if not csv_files:
        raise FileNotFoundError(f"No CSV files found in {folder_path}")

    summary_labels = ["MEAN_ACROSS_SLICES", "STD_ACROSS_SLICES", "OVERALL_3D_VOLUME"]
    dfs = []

    for path in csv_files:
        fname = os.path.basename(path)
        subject, region = extract_metadata_from_filename(fname)

        df = pd.read_csv(path)
        df = df[~df["Slice"].isin(summary_labels)].copy()

        if df.empty:
            continue

        def extract_idx(val):
            match = re.search(r"\d+", str(val))
            return int(match.group(0)) if match else None

        df["Slice_Idx"] = df["Slice"].apply(extract_idx)
        df["Subject"] = subject
        df["Region"] = region
        df["Technique"] = technique_name
        df["Source_File"] = fname

        dfs.append(df)

    if not dfs:
        raise ValueError(f"No valid slice data found in {folder_path}")

    return pd.concat(dfs, ignore_index=True)


def get_pvalue_annotation(p_val):
    """Converts p-value float into scientific star annotation."""
    if p_val < 0.0001:
        return "****"
    elif p_val < 0.001:
        return "***"
    elif p_val < 0.01:
        return "**"
    elif p_val < 0.05:
        return "*"
    return "ns"


def plot_metric_violins(df, output_dir, metrics=None):
    """Generates styled violin plots with significance bars and paired connecting lines."""
    os.makedirs(output_dir, exist_ok=True)

    if metrics is None:
        metrics = ["GM_SNR", "WM_SNR", "CSF_SNR", "CNR_GM_WM", "CNR_WM_CSF"]

    region_order = [r for r in ["Cervical", "Thoracic", "Lumbar"] if r in df["Region"].unique()]
    techniques = ["Siemens", "SLR"]
    
    # Matching palette: Orange & Purple
    palette = {"Siemens": "#D4819E", "SLR": "#679469"}
    
    # Base layout style
    sns.set_theme(style="ticks", font="DejaVu Sans", font_scale=1.1)

    for metric in metrics:
        if metric not in df.columns:
            continue

        df[metric] = pd.to_numeric(df[metric], errors="coerce")
        df_clean = df.dropna(subset=[metric]).copy()

        fig, ax = plt.subplots(figsize=(8, 6.5), dpi=300)

        # 1. Side-by-side full violins with inner boxplots
        sns.violinplot(
            data=df_clean,
            x="Region",
            y=metric,
            hue="Technique",
            order=region_order,
            hue_order=techniques,
            palette=palette,
            split=False,
            inner="box",
            cut=0,
            linewidth=1.3,
            width=0.7,
            ax=ax
        )

        # Make violin bodies semi-transparent
        for poly in ax.collections:
            poly.set_alpha(0.7)

        # 2. Add individual slice points & connecting lines between pairs
        x_offsets = {"Siemens": -0.18, "SLR": 0.18}
        
        for r_idx, region in enumerate(region_order):
            region_df = df_clean[df_clean["Region"] == region]

            pts_t1 = region_df[region_df["Technique"] == techniques[0]]
            pts_t2 = region_df[region_df["Technique"] == techniques[1]]

            # Jittered scatter points for Technique 1 (Siemens / Orange)
            x_pos_t1 = np.random.normal(r_idx + x_offsets["Siemens"], 0.02, size=len(pts_t1))
            ax.scatter(
                x_pos_t1, pts_t1[metric], 
                color="#D4819E", edgecolor="black", linewidth=0.5, s=28, zorder=4, alpha=0.9
            )

            # Jittered scatter points for Technique 2 (SLR / Purple)
            x_pos_t2 = np.random.normal(r_idx + x_offsets["SLR"], 0.02, size=len(pts_t2))
            ax.scatter(
                x_pos_t2, pts_t2[metric], 
                color="#679469", edgecolor="black", linewidth=0.5, s=28, zorder=4, alpha=0.9
            )

            # Draw light connecting lines between matching slice pairs
            paired_merge = pd.merge(
                pts_t1, pts_t2, 
                on=["Subject", "Region", "Slice_Idx"], 
                suffixes=("_t1", "_t2")
            )
            for _, pair in paired_merge.iterrows():
                y1 = pair[f"{metric}_t1"]
                y2 = pair[f"{metric}_t2"]
                ax.plot(
                    [r_idx + x_offsets["Siemens"], r_idx + x_offsets["SLR"]],
                    [y1, y2],
                    color="gray",
                    alpha=0.22,
                    linewidth=0.8,
                    zorder=3
                )

            # 3. Compute Stats & Draw Significance Brackets
            # 3. Compute Stats, Color Brackets, and Print Results to Terminal
            vals_t1 = pts_t1[metric].values  # Siemens
            vals_t2 = pts_t2[metric].values  # SLR

            if len(vals_t1) > 0 and len(vals_t2) > 0:
                if len(vals_t1) != len(vals_t2):
                    raise ValueError("Paired data must have the exact same number of observations.")

                all_differences_zero = all(x == y for x, y in zip(vals_t1, vals_t2))
                if all_differences_zero:
                    p_val = 1.0
                else:
                    _, p_val = stats.wilcoxon(vals_t2, vals_t1)

                annot = get_pvalue_annotation(p_val)

                # Determine direction: calculate medians to see who performs better
                med_siemens = np.median(vals_t1)
                med_slr = np.median(vals_t2)
                diff = med_slr - med_siemens

                # Define status & colors based on performance and statistical significance
                if p_val < 0.05:
                    if diff > 0:
                        winner = "SLR (Custom) is BETTER"
                        annotation_color = "#2E7D32"  # Green
                    else:
                        winner = "Siemens (Standard) is BETTER"
                        annotation_color = "#D4819E"  # Red
                else:
                    winner = "No Significant Difference (Equivalent)"
                    annotation_color = "black"

                # Print explicit verdict to console
                print(f"  [{metric} | {region}] Siemens Med: {med_siemens:.2f} | SLR Med: {med_slr:.2f} | p={p_val:.4f} ({annot}) --> {winner}")

                # Bracket position
                y_max_pair = max(vals_t1.max(), vals_t2.max())
                bar_y = y_max_pair * 1.08
                bar_h = bar_y * 0.02

                x1 = r_idx + x_offsets["Siemens"]
                x2 = r_idx + x_offsets["SLR"]

                # Bracket line with dynamic color
                ax.plot([x1, x1, x2, x2], [bar_y - bar_h, bar_y, bar_y, bar_y - bar_h], color=annotation_color, lw=1.4)
                
                # Asterisk / 'ns' text with dynamic color
                ax.text(
                    (x1 + x2) * 0.5, bar_y + bar_h * 0.3, 
                    annot, ha="center", va="bottom", 
                    fontsize=12, fontweight="bold", color=annotation_color
                )

        # 4. Spines & Headroom
        sns.despine(top=True, right=True)
        ax.spines["left"].set_linewidth(1.2)
        ax.spines["bottom"].set_linewidth(1.2)

        global_max = df_clean[metric].max()
        ax.set_ylim(bottom=-0.05 * global_max, top=global_max * 1.45)
        ax.margins(x=0.1)

        # 5. Labels & Formatting
        formatted_title = metric.replace("_", " ")
        ax.set_ylabel(formatted_title, fontsize=12, fontweight="bold", labelpad=8)
        ax.set_xlabel("", fontsize=12)
        ax.tick_params(labelsize=11, length=5, width=1.2)

        # Legend
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(
            handles[:2], labels[:2], 
            frameon=False, 
            loc="upper right",
            fontsize=11
        )

        plt.tight_layout()

        out_fname = os.path.join(output_dir, f"violin_{metric.lower()}.png")
        plt.savefig(out_fname, dpi=300)
        plt.close(fig)
        print(f"--> Saved plot: {out_fname}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate journal-style violin plots comparing Siemens vs SLR metrics."
    )
    parser.add_argument(
        "--siemens_dir",
        required=True,
        help="Path to Siemens CSV folder (e.g., qualitymetricsresults/siemens)"
    )
    parser.add_argument(
        "--slr_dir",
        required=True,
        help="Path to SLR CSV folder (e.g., qualitymetricsresults/SLR)"
    )
    parser.add_argument(
        "-o", "--output_dir",
        default="./plots",
        help="Directory to save the plots"
    )

    args = parser.parse_args()

    df_siemens = load_technique_directory(args.siemens_dir, technique_name="Siemens")
    df_slr = load_technique_directory(args.slr_dir, technique_name="SLR")

    df_all = pd.concat([df_siemens, df_slr], ignore_index=True)

    metrics_to_plot = ["GM_SNR", "WM_SNR", "CSF_SNR", "CNR_GM_WM", "CNR_WM_CSF"]
    plot_metric_violins(df_all, args.output_dir, metrics=metrics_to_plot)


if __name__ == "__main__":
    main()