from __future__ import annotations
import argparse
import re
from math import isnan
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).resolve().parent.parent

WELL_RE = re.compile(r"^Well([A-Z])(\d{2})")

# Strain order and treatment mapping (column-based bead sizes)
# Rows A-D -> SP286
STRAIN_ORDER = ["SP286"]

# Order for plotting: show NO bead first (left), then small->large beads
TREATMENT_ORDER = [
    "NO bead",
    "1 mm bead",
    "1.5 mm bead",
    "3 mm bead",
    "4.5 mm bead",
]


def violin_grouped_with_means(df: pd.DataFrame, y: str, title: str, outpath: Path) -> None:
    outpath.parent.mkdir(parents=True, exist_ok=True)

    dfp = df.dropna(subset=[y]).copy()
    if dfp.empty:
        print(f"No data to plot for {y}")
        return

    sns.set(style="whitegrid")

    fig, ax = plt.subplots(figsize=(12, 6))

    # Compute group statistics (mean, std) for each strain x treatment
    stats_df = (
        dfp.groupby(["strain", "treatment"], observed=True)[y]
        .agg(["mean", "std", "count"])
        .reset_index()
    )

    # Plot grouped bar chart: for each strain, draw a set of dodged bars for treatments
    x_centers = {strain: i for i, strain in enumerate(STRAIN_ORDER)}
    n_hue = len(TREATMENT_ORDER)
    total_width = 0.8
    bar_width = total_width / n_hue * 0.9
    step = total_width / n_hue

    for i_t, treatment in enumerate(TREATMENT_ORDER):
        bar_x = []
        bar_h = []
        bar_err = []
        for strain in STRAIN_ORDER:
            row = stats_df[(stats_df["strain"] == strain) & (stats_df["treatment"] == treatment)]
            if not row.empty and row["count"].iloc[0] > 0:
                mu = row["mean"].iloc[0]
                sd = row["std"].iloc[0] if not isnan(row["std"].iloc[0]) else 0.0
            else:
                mu = 0.0
                sd = 0.0
            base = x_centers[strain]
            offset = -total_width / 2 + (i_t + 0.5) * step
            x = base + offset
            bar_x.append(x)
            bar_h.append(mu)
            bar_err.append(sd)

        ax.bar(bar_x, bar_h, yerr=bar_err, width=bar_width, label=treatment, alpha=0.9)

    ax.set_title(title)
    ax.set_xlabel("Strain")
    ax.set_xticks(list(x_centers.values()))
    ax.set_xticklabels(list(x_centers.keys()))

    # Use user-friendly axis labels with units when possible
    if y == "length":
        ylabel = "length (µm)"
    elif y == "area":
        ylabel = "area (µm²)"
    else:
        ylabel = y

    ax.set_ylabel(ylabel)

    # ---- Mean annotations ----
    means = dfp.groupby(["strain", "treatment"], observed=True)[y].mean().reset_index()

    # track maximum y used for significance labels so we can expand the
    # axis top limit to make space for stars / letters
    max_sig_y = 0.0

    # Compute x positions that match seaborn's dodge layout
    x_centers = {strain: i for i, strain in enumerate(STRAIN_ORDER)}
    n_hue = len(TREATMENT_ORDER)
    total_width = 0.8  # seaborn category width approx
    step = total_width / n_hue

    # Global y-range for offsetting annotations
    y_min = dfp[y].min()
    y_max = dfp[y].max()
    y_span = (y_max - y_min) if (y_max is not None and y_min is not None) else 0.0
    ann_offset = y_span * 0.03 if y_span else 0.1

    # Delegate grouping-letter computation to the shared helper.
    # Replace previous external grouping-letter helper with a per-strain
    # one-way ANOVA followed by pairwise comparisons vs the NO bead control.
    # sig_map maps (strain, treatment) -> significance label (e.g. '*', '**', '***' or 'ns').
    sig_map = {}

    def pval_to_stars(p):
        if p is None or (isinstance(p, float) and (pd.isna(p) or np.isnan(p))):
            return None
        try:
            p = float(p)
        except Exception:
            return None
        if p < 0.001:
            return "***"
        if p < 0.01:
            return "**"
        if p < 0.05:
            return "*"
        return "ns"

    # For each strain, run one-way ANOVA across treatments (if possible)
    for strain in STRAIN_ORDER:
        sdf = dfp[dfp["strain"] == strain]
        if sdf.empty:
            continue

        # identify groups present
        groups_present = [t for t in TREATMENT_ORDER if not sdf[sdf["treatment"] == t].empty]

        # pairwise comparisons vs NO bead control
        control_vals = sdf.loc[sdf["treatment"] == "NO bead", y].dropna().values
        pvals = []
        treatments_tested = []
        for t in groups_present:
            if t == "NO bead":
                # control - mark as ns placeholder
                sig_map[(strain, t)] = None
                continue
            vals = sdf.loc[sdf["treatment"] == t, y].dropna().values
            if len(vals) < 2 or len(control_vals) < 2:
                p = np.nan
            else:
                try:
                    stat, p = stats.ttest_ind(
                        vals, control_vals, equal_var=False, nan_policy="omit"
                    )
                except Exception:
                    p = np.nan
            pvals.append(p)
            treatments_tested.append(t)

        # multiple testing correction (Holm) if available, else do simple Bonferroni as fallback
        if pvals:
            pvals_arr = np.array([np.nan if p is None else p for p in pvals], dtype=float)
            mask = ~np.isnan(pvals_arr)
            adj = np.full_like(pvals_arr, np.nan)
            if mask.any():
                _, p_adj, _, _ = multipletests(pvals_arr[mask], method="holm")
                adj[mask] = p_adj

            # map adjusted p-values to stars
            for t, _p_raw, p_adj in zip(treatments_tested, pvals_arr, adj, strict=False):
                star = pval_to_stars(p_adj)
                sig_map[(strain, t)] = star

        # ensure control has an entry (None)
        sig_map.setdefault((strain, "NO bead"), None)

    for _, row in means.iterrows():
        strain = row["strain"]
        treatment = row["treatment"]
        mu = row[y]

        if pd.isna(mu) or strain not in x_centers:
            continue

        base = x_centers[strain]
        j = TREATMENT_ORDER.index(str(treatment))

        # offset each hue category within the strain bin
        offset = -total_width / 2 + (j + 0.5) * step
        x = base + offset

        # Lookup std for this group so label placement avoids the error bar
        row_stat = stats_df[(stats_df["strain"] == strain) & (stats_df["treatment"] == treatment)]
        sd_val = 0.0
        if not row_stat.empty:
            try:
                sd_val = (
                    float(row_stat["std"].iloc[0]) if not isnan(row_stat["std"].iloc[0]) else 0.0
                )
            except Exception:
                sd_val = 0.0

        # place mean value inside the bar (half height) for readability
        if mu > 0:
            y_text_inside = mu * 0.5
        else:
            y_text_inside = mu + ann_offset

        ax.text(
            x,
            y_text_inside,
            f"{mu:.2f}",
            ha="center",
            va="center",
            fontsize=8,
            color="white",
            fontweight="bold",
        )

        # significance label (placed above the error bar / bar top)
        sig_label = sig_map.get((str(strain), str(treatment)), None)
        sig_y = (
            mu
            + (sd_val if sd_val > 0 else 0.0)
            + max(ann_offset * 0.6, (sd_val if sd_val > 0 else mu * 0.03))
        )
        # Skip showing 'ns' on the NO bead (control) bars to reduce clutter
        if sig_label and not (str(treatment) == "NO bead" and sig_label == "ns"):
            ax.text(
                x,
                sig_y,
                sig_label,
                ha="center",
                va="bottom",
                fontsize=8,
                fontweight="bold",
            )
            if sig_y and (not pd.isna(sig_y)):
                try:
                    max_sig_y = max(max_sig_y, float(sig_y))
                except Exception:
                    pass

    # reserve room for legend to the right and a taller top margin for stats
    ax.legend(title="Condition", bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.subplots_adjust(right=0.78, top=0.92)

    # expand y-axis if significance labels would be clipped
    try:
        cur_min, cur_max = ax.get_ylim()
        if max_sig_y > 0 and max_sig_y > cur_max:
            ax.set_ylim(bottom=max(0, cur_min), top=max_sig_y * 1.12)
    except Exception:
        pass

    fig.tight_layout()
    fig.savefig(str(outpath), dpi=300)
    plt.close(fig)

    print(f"Saved: {outpath}")


def parse_well_from_name(name: str) -> tuple[str | None, int | None]:
    m = WELL_RE.match(name)
    if not m:
        return None, None
    return m.group(1), int(m.group(2))


def map_strain_and_treatment(letter: str, num: int) -> tuple[str, str]:
    """Map a well (row letter + column number) to (strain, treatment).

    Updated rules (column determines bead size):
    - Rows A-D -> SP286
    - Cols 1-2  -> 1 mm bead
    - Cols 3-4  -> 1.5 mm bead
    - Cols 5-6  -> 3 mm bead
    - Cols 7-8  -> 4.5 mm bead
    - Cols 9-10 -> NO bead
    """
    letter = letter.upper().strip()

    if letter in {"A", "B", "C", "D"}:
        strain = "SP286"
    else:
        return "unknown", "unknown"

    # columns 1..10 only
    if not (1 <= num <= 10):
        return "unknown", "unknown"

    # treatment by explicit column ranges
    if 1 <= num <= 2:
        treatment = "1 mm bead"
    elif 3 <= num <= 4:
        treatment = "1.5 mm bead"
    elif 5 <= num <= 6:
        treatment = "3 mm bead"
    elif 7 <= num <= 8:
        treatment = "4.5 mm bead"
    else:  # 9-10
        treatment = "NO bead"

    return strain, treatment


def find_matching_csvs(processed_dir: Path) -> list[Path]:
    files: list[Path] = []
    for p in processed_dir.glob("*.csv"):
        if not p.name.startswith("Well"):
            continue
        letter, num = parse_well_from_name(p.name)
        if letter is None or num is None:
            continue
        strain, treatment = map_strain_and_treatment(letter, num)
        if strain != "unknown":
            files.append(p)
    return sorted(files)


def read_one_csv(p: Path) -> pd.DataFrame:
    df = pd.read_csv(p)
    # Accept a few common column name variants produced by different
    # processing pipelines. Normalize to the canonical names used by the
    # plotting code: 'axis_major_length' and 'area'.
    axis_candidates = ["axis_major_length", "axis_major_length_um", "axis_major_length_px"]
    area_candidates = ["area", "area_um2", "area_px2"]

    axis_col = next((c for c in axis_candidates if c in df.columns), None)
    if axis_col is None:
        raise KeyError(
            f"{p.name} missing required column: axis_major_length (or variants). Columns: "
            f"{', '.join(df.columns)}"
        )

    out = pd.DataFrame(
        {
            "axis_major_length": pd.to_numeric(df[axis_col], errors="coerce"),
        }
    )

    # area is optional but desired; accept variants
    area_col = next((c for c in area_candidates if c in df.columns), None)
    if area_col is not None:
        out["area"] = pd.to_numeric(df[area_col], errors="coerce")
    else:
        out["area"] = pd.NA

    # Provide a canonical 'length' column (user-facing name) that aliases
    # the measured axis major length. We assume values are already in
    # micrometers when the column name contains '_um' or when the
    # canonical name is used; if your files are in pixels, supply a
    # conversion factor separately (px -> µm) and we can apply it here.
    out["length"] = out["axis_major_length"]

    return out


def build_combined_dataframe(csv_paths: list[Path]) -> pd.DataFrame:
    parts = []
    for p in csv_paths:
        letter, num = parse_well_from_name(p.name)
        if letter is None or num is None:
            continue

        strain, treatment = map_strain_and_treatment(letter, num)
        if strain == "unknown":
            continue

        df = read_one_csv(p)
        df["source_file"] = p.name
        df["well_letter"] = letter
        df["well_num"] = num
        df["strain"] = strain
        df["treatment"] = treatment
        parts.append(df)

    if not parts:
        return pd.DataFrame()

    out = pd.concat(parts, ignore_index=True)

    # Expose a canonical 'length' column (alias for axis_major_length)
    out["length"] = out["axis_major_length"]

    out["strain"] = pd.Categorical(out["strain"], categories=STRAIN_ORDER, ordered=True)
    out["treatment"] = pd.Categorical(out["treatment"], categories=TREATMENT_ORDER, ordered=True)
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="DIC-only: violin plots faceted by strain and bead treatment"
    )
    parser.add_argument(
        "--processed-dir",
        "-d",
        type=Path,
        default=Path.cwd() / "processed",
        help="Path to directory containing processed CSV files",
    )
    args = parser.parse_args(argv)

    processed = args.processed_dir

    if not processed.exists() or not processed.is_dir():
        print(f"Processed directory not found: {processed}")
        return 0

    csvs = find_matching_csvs(processed)
    if not csvs:
        print(f"No matching Well*.csv files found in {processed}")
        return 0

    df = build_combined_dataframe(csvs)
    if df.empty:
        print("No data after mapping wells to strain/treatment. Check filenames + mapping rules.")
        return 0

    csv_out = REPO_ROOT / "data" / "microscopy" / "combined_dic_measurements_96well.csv"
    csv_out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_out, index=False)
    print(f"Saved combined data: {csv_out}")

    fig_dir = REPO_ROOT / "figures" / "figure-4"
    fig_dir.mkdir(parents=True, exist_ok=True)

    violin_grouped_with_means(
        df,
        y="length",
        title="DIC length (µm) by strain and bead treatment (means labeled)",
        outpath=fig_dir / "bar_length_grouped.svg",
    )

    violin_grouped_with_means(
        df,
        y="area",
        title="DIC area (µm²) by strain and bead treatment (means labeled)",
        outpath=fig_dir / "bar_area_grouped.svg",
    )

    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
