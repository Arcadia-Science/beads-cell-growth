from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent

PROCESSED_CSV_PATHS = {
    "96-well": REPO_ROOT / "data" / "microscopy" / "combined_dic_measurements_96well.csv",
    "24-well": REPO_ROOT / "data" / "microscopy" / "combined_dic_measurements_24well.csv",
    "ttubes": REPO_ROOT / "data" / "microscopy" / "combined_dic_measurements_ttubes.csv",
}
OD_PATH = REPO_ROOT / "data" / "plate-reader" / "Baseline_ODs_stdev.csv"
OUTPUT_CSV = REPO_ROOT / "data" / "aggregated_summary.csv"

STRAIN_FIG_DIRS = {
    "SP286": REPO_ROOT / "figures" / "figure-5",
}


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


def plot_facets(summary_csv: str, strain_fig_dirs: dict[str, Path]) -> None:
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
        fig_dir = strain_fig_dirs.get(strain)
        if fig_dir is None:
            print(f"No output directory configured for strain {strain}, skipping plot.")
            continue
        fig_dir.mkdir(parents=True, exist_ok=True)

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
        out_png = fig_dir / f"aggregated_summary_{strain}.png"
        plt.savefig(out_png)
        plt.close(fig)
        print(f"Saved: {out_png}")


def main():
    """Aggregate microscopy data from multiple sources into a single CSV file."""
    dfs = []
    for experiment, csv_path in PROCESSED_CSV_PATHS.items():
        if not csv_path.exists():
            print(f"Warning: {csv_path} not found, skipping.")
            continue

        df = pd.read_csv(csv_path)
        df["experiment"] = experiment
        if experiment == "96-well":
            if "treatment" in df.columns:
                t = df["treatment"].astype(str).str.lower()
                df["bead_present"] = t.apply(lambda x: True if "4.5 mm bead" in x else False)
                df = df[t.isin(["no bead", "4.5 mm bead"])]
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
        return
    big = pd.concat(dfs, ignore_index=True)
    n_96well = (big["experiment"] == "96-well").sum() if "experiment" in big.columns else 0
    print(f"Rows with experiment='96-well': {n_96well}")
    summary = summarize(big)
    od_df = pd.read_csv(OD_PATH)
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
    summary.to_csv(OUTPUT_CSV, index=False)
    print(f"Summary merged with OD data and exported to {OUTPUT_CSV}")
    plot_facets(str(OUTPUT_CSV), STRAIN_FIG_DIRS)


if __name__ == "__main__":
    main()
