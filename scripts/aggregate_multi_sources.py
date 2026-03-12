import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    df2 = df.dropna(subset=["volume_ml"]).copy()
    if "strain" not in df2.columns:
        df2["strain"] = "unknown"
    group_cols = (
        ["experiment", "strain", "bead_present", "volume_ml"]
        if "experiment" in df2.columns
        else ["strain", "bead_present", "volume_ml"]
    )
    g = df2.groupby(group_cols, observed=True)
    out = g.agg(
        length_mean=("length", "mean"),
        length_stdev=("length", "std"),
        area_mean=("area", "mean"),
        area_stdev=("area", "std"),
        n_cells=("length", "count"),
    ).reset_index()
    return out


def plot_facets(summary_csv, output_png):

    df = pd.read_csv(summary_csv)
    df["volume_ml"] = pd.to_numeric(df["volume_ml"], errors="coerce")
    df["bead_present"] = df["bead_present"].astype(bool)
    strains = df["strain"].unique()
    experiments = df["experiment"].unique() if "experiment" in df.columns else [""]
    bead_labels = [False, True]
    bead_names = {False: "no bead", True: "bead"}
    x_volumes = [1, 2, 3, 4, 5]
    facets = [
        ("length_mean", "length", "length_stdev"),
        ("area_mean", "area", "area_stdev"),
        ("OD_morning_mean", "OD morning", "OD_morning_stdev"),
        ("OD_afternoon_mean", "OD afternoon", "OD_afternoon_stdev"),
    ]
    for strain in strains:
        fig, axes = plt.subplots(2, 2, figsize=(14, 8))
        df_strain = df[df["strain"] == strain]
        for idx, (mean_col, title, stdev_col) in enumerate(facets):
            ax = axes[idx // 2, idx % 2]
            for exp in experiments:
                for bead in bead_labels:
                    d = df_strain[
                        (df_strain["experiment"] == exp) & (df_strain["bead_present"] == bead)
                    ]
                    d = d[d["volume_ml"].isin(x_volumes)]
                    if mean_col in d.columns and stdev_col in d.columns and not d.empty:
                        ax.errorbar(
                            d["volume_ml"],
                            d[mean_col],
                            yerr=d[stdev_col],
                            marker="o",
                            label=f"{exp} ({bead_names[bead]})",
                            capsize=4,
                        )
            ax.set_title(title)
            ax.set_xlabel("Volume (mL)")
            ax.set_xticks(x_volumes)
            ax.grid(True)
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper right")
        fig.suptitle(f"Strain: {strain}")
        fig.tight_layout()
        out_png = output_png.replace(".png", f"_{strain}.png")
        plt.savefig(out_png)
        plt.close(fig)


def main():
    """Aggregate microscopy data from multiple sources into a single CSV file.

    Usage: python aggregate_multi_sources.py <processed_folder1> <processed_folder2> ... <output_csv>

    Args:
        processed_folder1: Path to the first processed folder containing microscopy data.
        processed_folder2: Path to the second processed folder containing microscopy data.
        ...: Additional processed folders containing microscopy data.
        output_csv: Path to the output CSV file.
    """
    if len(sys.argv) < 3:
        print(
            "Usage: python aggregate_multi_sources.py <processed_folder1> <processed_folder2> ... <output_csv>"
        )
        sys.exit(1)
    *input_folders, output_csv = sys.argv[1:]
    dfs = []
    experiment_labels = ["96-well", "24-well", "ttubes"]
    for i, folder in enumerate(input_folders):
        folder_path = Path(folder)
        experiment = experiment_labels[i] if i < len(experiment_labels) else f"exp_{i + 1}"
        for csv_file in folder_path.glob("*.csv"):
            df = pd.read_csv(csv_file)
            df["experiment"] = experiment
            if experiment == "96-well":
                if "treatment" in df.columns:
                    t = df["treatment"].astype(str).str.lower()
                    df["bead_present"] = t.apply(lambda x: True if "1 mm bead" in x else False)
                elif "bead_present" in df.columns:
                    df["bead_present"] = (
                        df["bead_present"]
                        .astype(str)
                        .str.lower()
                        .map(
                            {
                                "no": False,
                                "1": True,
                                "yes": True,
                                "true": True,
                                "false": False,
                                "0": False,
                            }
                        )
                    )
                df["volume_ml"] = 1.0
            else:
                if "bead_present" in df.columns:
                    df["bead_present"] = (
                        df["bead_present"]
                        .astype(str)
                        .str.lower()
                        .map(
                            {
                                "no": False,
                                "1": True,
                                "yes": True,
                                "true": True,
                                "false": False,
                                "0": False,
                            }
                        )
                    )
            dfs.append(df)
    if not dfs:
        print("No CSV files found in input folders.")
        sys.exit(1)
    big = pd.concat(dfs, ignore_index=True)
    # Check for 96-well experiment rows
    n_96well = (big["experiment"] == "96-well").sum() if "experiment" in big.columns else 0
    print(f"Rows with experiment='96-well': {n_96well}")
    summary = summarize(big)
    # Merge with Input_ODs.csv
    od_path = "/Users/roman/Repositories/2026-pombe-beads/Data_analysis/Input_ODs.csv"
    od_df = pd.read_csv(od_path)
    od_df.columns = [c.strip() for c in od_df.columns]
    od_df = od_df.rename(columns={"beads": "bead_present", "volume": "volume_ml"})
    od_df["bead_present"] = (
        od_df["bead_present"]
        .astype(str)
        .str.lower()
        .map({"true": True, "false": False, "1": True, "0": False})
    )
    od_df["volume_ml"] = pd.to_numeric(od_df["volume_ml"], errors="coerce")
    merge_keys = ["experiment", "strain", "bead_present", "volume_ml"]
    summary = summary.merge(od_df, on=merge_keys, how="left")
    summary.to_csv(output_csv, index=False)
    print(f"Summary merged with Input_ODs.csv and exported to {output_csv}")


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1].endswith(".csv"):
        plot_facets(sys.argv[1], sys.argv[2])
    else:
        main()
