from __future__ import annotations
import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

WELL_RE = re.compile(r"^Well([A-Z])(\d{2})")

# Strain order remains the two strains used in these experiments
STRAIN_ORDER = ["dea2^", "SP286"]

# Treatments for the supplements experiment are combinations of supplements
# mapped by well column and row. We enumerate expected labels here so
# plotting can treat them as ordered categoricals.
TREATMENT_ORDER = [
    "NONE",
    "D",
    "V",
    "A",
    "D+A",
    "V+A",
    "D+V",
    "D+V+A",
]


def violin_grouped_with_means(
    df: pd.DataFrame, y: str, title: str, outpath: Path, clusters: int = 2
) -> None:
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

    # Per-strain one-way ANOVA and pairwise comparisons versus the 'NONE'
    # (no supplement) control with Holm-corrected p-values.
    sig_map: dict[tuple[str, str], str] = {}

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

    # collect per-strain results for optional CSV export
    pval_rows: list[dict] = []

    for strain in dfp["strain"].cat.categories:
        df_strain = dfp[dfp["strain"] == strain]
        if df_strain.empty:
            continue

        # groups present for ANOVA
        groups_present = [
            t for t in TREATMENT_ORDER if not df_strain[df_strain["treatment"] == t].empty
        ]

        # overall one-way ANOVA across present groups (if >=2)
        anova_p = None
        try:
            arrays = [
                df_strain.loc[df_strain["treatment"] == t, y].dropna().values
                for t in groups_present
                if len(df_strain.loc[df_strain["treatment"] == t, y].dropna()) > 0
            ]
            if len(arrays) >= 2:
                anova_res = stats.f_oneway(*arrays)
                anova_p = (
                    float(anova_res.pvalue) if hasattr(anova_res, "pvalue") else float(anova_res[1])
                )
        except Exception:
            anova_p = None

        # pairwise comparisons vs NONE control
        control_vals = df_strain.loc[df_strain["treatment"] == "NONE", y].dropna().values
        pvals = []
        treatments_tested = []
        for t in groups_present:
            if t == "NONE":
                sig_map[(strain, t)] = None
                continue
            vals = df_strain.loc[df_strain["treatment"] == t, y].dropna().values
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
            pval_rows.append(
                {
                    "strain": strain,
                    "treatment": t,
                    "anova_p": anova_p,
                    "p_raw": p,
                    "p_adj": None,
                }
            )

        # multiple testing correction across treatments for this strain
        if pvals:
            pvals_arr = np.array([np.nan if p is None else p for p in pvals], dtype=float)
            mask = ~np.isnan(pvals_arr)
            adj = np.full_like(pvals_arr, np.nan)
            if mask.any():
                _, p_adj, _, _ = multipletests(pvals_arr[mask], method="holm")
                adj[mask] = p_adj

            for t, _p_raw, p_adj in zip(treatments_tested, pvals_arr, adj, strict=False):
                star = pval_to_stars(p_adj)
                sig_map[(strain, t)] = star
                # update pval_rows with adjusted p
                for r in pval_rows:
                    if r["strain"] == strain and r["treatment"] == t:
                        r["p_adj"] = float(p_adj) if not np.isnan(p_adj) else None
                        break

    # write p-values CSV next to the output if we have any results
    try:
        if pval_rows:
            pvals_df = pd.DataFrame(pval_rows)
            csv_out = outpath.parent / f"one_way_{y}_pvalues.csv"
            pvals_df.to_csv(csv_out, index=False)
            print(f"Wrote p-values summary: {csv_out}")
    except Exception:
        pass

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
                sd = row["std"].iloc[0] if not np.isnan(row["std"].iloc[0]) else 0.0
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
    # Replace caret '^' with Greek delta 'Δ' for display only
    display_names = [s.replace("^", "Δ") for s in list(x_centers.keys())]
    ax.set_xticklabels(display_names)

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

    # Note: we purposely do not compute ANOVA or pairwise-star annotations
    # here because we already produced grouping letters for each strain.

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
                    float(row_stat["std"].iloc[0]) if not np.isnan(row_stat["std"].iloc[0]) else 0.0
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

    Mapping rules for the supplements experiment (user-specified):
    - Strain is determined by the column block: columns 1-6 -> dea2^, columns 7-12 -> SP286
    - Within each 6-column block the treatment is determined by the row letter
      and the column sub-block (1-3 vs 4-6):
        columns 1-3 (or 7-9):
          A -> "NONE", B -> "D", C -> "V", D -> "A"
        columns 4-6 (or 10-12):
          A -> "D+A", B -> "V+A", C -> "D+V", D -> "D+V+A"
    The returned `treatment` is a short label like "D+V".
    """
    letter = letter.upper().strip()

    # validate column range for this plate (expect 1..12)
    if not (1 <= num <= 12):
        return "unknown", "unknown"

    # strain by column block: 1-6 => dea2^, 7-12 => SP286
    if 1 <= num <= 6:
        strain = "dea2^"
    else:
        strain = "SP286"

    # position within the 6-column block (1..6)
    pos = ((num - 1) % 6) + 1

    # mapping for columns 1-3 (pos 1..3)
    map_1_3 = {
        "A": "NONE",
        "B": "D",
        "C": "V",
        "D": "A",
    }

    # mapping for columns 4-6 (pos 4..6)
    map_4_6 = {
        "A": "D+A",
        "B": "V+A",
        "C": "D+V",
        "D": "D+V+A",
    }

    if pos <= 3:
        treatment = map_1_3.get(letter, "unknown")
    else:
        treatment = map_4_6.get(letter, "unknown")

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
            "{', '.join(df.columns)}"
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
    # derive structured metadata from the treatment string: numeric volume and bead presence
    # treatment strings are like '1 mL no bead' or '3 mL bead' so extract leading number
    vol_match = out["treatment"].str.extract(r"^(\d+)")
    out["volume_ml"] = pd.to_numeric(vol_match[0], errors="coerce")
    # map volumes 1..5 to group letters a..e
    out["volume_group"] = out["volume_ml"].apply(
        lambda v: chr(96 + int(v)) if not pd.isna(v) and 1 <= int(v) <= 26 else pd.NA
    )
    out["bead_present"] = out["treatment"].str.contains(r"\\bbead\\b")

    out["strain"] = pd.Categorical(out["strain"], categories=STRAIN_ORDER, ordered=True)
    out["treatment"] = pd.Categorical(out["treatment"], categories=TREATMENT_ORDER, ordered=True)
    return out


def violin_faceted(df: pd.DataFrame, y: str, title: str, outpath: Path) -> None:
    outpath.parent.mkdir(parents=True, exist_ok=True)

    dfp = df.dropna(subset=[y]).copy()
    if dfp.empty:
        print(f"No data to plot for {y}")
        return

    sns.set(style="whitegrid")

    # Facet rows = volume (1-5 mL), cols = strain
    volumes = [f"{v} mL" for v in (1, 2, 3, 4, 5)]
    strains = STRAIN_ORDER
    n_rows = len(volumes)
    n_cols = len(strains)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 2.5), squeeze=False)

    # Precompute stats
    stats_df = (
        dfp.groupby(["strain", "treatment"], observed=True)[y]
        .agg(["mean", "std", "count"])
        .reset_index()
    )

    def get_stats(strain: str, treatment: str) -> tuple[float, float, int]:
        row = stats_df[(stats_df["strain"] == strain) & (stats_df["treatment"] == treatment)]
        if not row.empty and row["count"].iloc[0] > 0:
            return (
                float(row["mean"].iloc[0]),
                float(row["std"].iloc[0]) if not np.isnan(row["std"].iloc[0]) else 0.0,
                int(row["count"].iloc[0]),
            )
        return 0.0, 0.0, 0

    # Compute pairwise bead vs no-bead p-values per strain-volume
    sig_map: dict[tuple[str, str], str] = {}
    for strain in dfp["strain"].cat.categories:
        df_strain = dfp[dfp["strain"] == strain]
        for vol in volumes:
            t_no = f"{vol} no bead"
            t_yes = f"{vol} bead"
            g_no = df_strain[df_strain["treatment"] == t_no][y].dropna().values
            g_yes = df_strain[df_strain["treatment"] == t_yes][y].dropna().values
            label_no = "ns"
            label_yes = "ns"
            if len(g_no) > 0 and len(g_yes) > 0:
                try:
                    _, pval = stats.ttest_ind(g_no, g_yes, equal_var=False, nan_policy="omit")
                    pval = float(pval)
                    if pval <= 0.001:
                        label_yes = "***"
                    elif pval <= 0.01:
                        label_yes = "**"
                    elif pval <= 0.05:
                        label_yes = "*"
                    else:
                        label_yes = "ns"
                except Exception:
                    label_yes = "ns"
            sig_map[(str(strain), t_no)] = label_no
            sig_map[(str(strain), t_yes)] = label_yes

    for i_row, vol in enumerate(volumes):
        for j_col, strain in enumerate(strains):
            ax = axes[i_row][j_col]

            t_no = f"{vol} no bead"
            t_yes = f"{vol} bead"

            mu_no, sd_no, cnt_no = get_stats(strain, t_no)
            mu_yes, sd_yes, cnt_yes = get_stats(strain, t_yes)

            if cnt_no == 0 and cnt_yes == 0:
                ax.text(0.5, 0.5, "no data", ha="center", va="center", fontsize=8)
                ax.set_xticks([])
                ax.set_ylim(bottom=0)
            else:
                x = [0, 1]
                heights = [mu_no, mu_yes]
                errs = [sd_no, sd_yes]
                ax.bar(x, heights, yerr=errs, width=0.6, color=["C0", "C1"], alpha=0.9)
                ax.set_xticks([])
                ax.set_ylim(bottom=0)

                # mean inside bars
                for xi, h in zip(x, heights, strict=False):
                    y_text_inside = (h * 0.5) if h > 0 else (h + 0.05)
                    ax.text(
                        xi,
                        y_text_inside,
                        f"{h:.2f}",
                        ha="center",
                        va="center",
                        fontsize=7,
                        color="white",
                        fontweight="bold",
                    )

                # significance: show above the bead bar if present
                sig_yes = sig_map.get((str(strain), t_yes), None)
                sig_no = sig_map.get((str(strain), t_no), None)
                if sig_yes and sig_yes != "ns":
                    sig_y = (
                        mu_yes
                        + (sd_yes if sd_yes > 0 else 0.0)
                        + max(mu_yes * 0.03, (sd_yes if sd_yes > 0 else 0.1))
                    )
                    ax.text(
                        1, sig_y, sig_yes, ha="center", va="bottom", fontsize=8, fontweight="bold"
                    )
                if sig_no and sig_no != "ns" and not t_no.endswith("no bead"):
                    sig_y = (
                        mu_no
                        + (sd_no if sd_no > 0 else 0.0)
                        + max(mu_no * 0.03, (sd_no if sd_no > 0 else 0.1))
                    )
                    ax.text(
                        0, sig_y, sig_no, ha="center", va="bottom", fontsize=8, fontweight="bold"
                    )

            # Titles and labels
            if i_row == 0:
                ax.set_title(strain.replace("^", "Δ"), fontsize=9)
            if j_col == 0:
                if y == "length":
                    ax.set_ylabel(f"{vol}\nlength (µm)")
                elif y == "area":
                    ax.set_ylabel(f"{vol}\narea (µm²)")
                else:
                    ax.set_ylabel(vol)
            else:
                ax.set_ylabel("")

    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    fig.savefig(str(outpath), dpi=300)
    plt.close(fig)
    print(f"Saved: {outpath}")


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
    parser.add_argument(
        "--out-dir",
        "-o",
        type=Path,
        default=None,
        help="Output directory (default: processed-dir)",
    )
    parser.add_argument(
        "--clusters",
        "-k",
        type=int,
        default=2,
        help="Number of clusters to form per-strain for grouping letters (default: 2)",
    )
    args = parser.parse_args(argv)

    processed = args.processed_dir
    out_dir = args.out_dir or processed

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

    combined_path = out_dir / "combined_dic_measurements.csv"
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(combined_path, index=False)
    print(f"Saved combined data: {combined_path}")

    violin_grouped_with_means(
        df,
        y="length",
        title="DIC length (µm) by strain and bead treatment (means labeled)",
        outpath=out_dir / "bar_length_grouped.svg",
        clusters=args.clusters,
    )

    violin_grouped_with_means(
        df,
        y="area",
        title="DIC area (µm²) by strain and bead treatment (means labeled)",
        outpath=out_dir / "bar_area_grouped.svg",
        clusters=args.clusters,
    )

    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
