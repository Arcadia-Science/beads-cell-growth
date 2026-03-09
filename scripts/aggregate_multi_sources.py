from __future__ import annotations
import argparse
import re
from pathlib import Path

import arcadia_pycolor as ap
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def find_combined_csv(base: Path) -> Path | None:
    p1 = base / "processed" / "combined_dic_measurements.csv"
    p2 = base / "combined_dic_measurements.csv"
    if p1.exists():
        return p1
    if p2.exists():
        return p2
    return None


def _arc_hex(name: str, fallback: str) -> str:
    """Return a hex color for a named arcadia color, or fallback."""
    val = getattr(ap, name, None)
    if val is None:
        return fallback
    try:
        return str(val)
    except Exception:
        return fallback


# Choose font: prefer 'Chateau' if installed, otherwise fall back to Arcadia default or
# a sensible matplotlib default. This sets the global rcParams so all plots use it.
def _choose_font(preferred: str = "Chateau") -> str:
    try:
        from matplotlib import font_manager as fm

        installed = [f.name for f in fm.fontManager.ttflist]
        for n in installed:
            if preferred.lower() in n.lower():
                return n
    except Exception:
        pass
    try:
        return getattr(ap.style_defaults, "DEFAULT_FONT", "DejaVu Sans")
    except Exception:
        return "DejaVu Sans"


chosen_font = _choose_font("Chateau")
plt.rcParams["font.family"] = chosen_font
try:
    plt.rcParams["font.size"] = ap.style_defaults.BASE_FONT_SIZE
except Exception:
    pass


def load_source(path: Path, label: str) -> pd.DataFrame:
    p = find_combined_csv(path)
    if p is None:
        raise FileNotFoundError(
            f"combined_dic_measurements.csv not found under {path} or {path}/processed"
        )
    df = pd.read_csv(p)
    df["__source"] = label
    return df


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    # ensure canonical columns exist: length, area, volume_ml, bead_present
    if "length" not in df.columns:
        # try axis_major_length
        for c in ("axis_major_length", "axis_major_length_um", "axis_major_length_px"):
            if c in df.columns:
                df["length"] = pd.to_numeric(df[c], errors="coerce")
                break
    if "area" not in df.columns:
        for c in ("area", "area_um2", "area_px2"):
            if c in df.columns:
                df["area"] = pd.to_numeric(df[c], errors="coerce")
                break
    if "volume_ml" not in df.columns:
        # attempt to parse from treatment column
        if "treatment" in df.columns:
            vol = df["treatment"].astype(str).str.extract(r"(\d+(?:\.\d+)?)\s*mL")
            df["volume_ml"] = pd.to_numeric(vol[0], errors="coerce")
    # Recompute bead_present from treatment to avoid relying on possibly
    # inconsistent precomputed values in source CSVs.
    if "treatment" in df.columns:
        t = df["treatment"].astype(str).str.lower()
        # detect explicit 'no bead' labels first
        no_bead = t.str.contains(r"\bno\s*bead\b", regex=True, na=False)
        has_bead = t.str.contains(r"\bbead\b", regex=True, na=False)
        # present only when 'bead' is present and not explicitly 'no bead'
        df["bead_present"] = has_bead & (~no_bead)
    else:
        df["bead_present"] = False
    return df


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    # group by strain, source, bead_present, volume_ml and compute mean and sd
    df2 = df.dropna(subset=["volume_ml"]).copy()
    # ensure strain exists
    if "strain" not in df2.columns:
        df2["strain"] = "unknown"
    g = df2.groupby(["strain", "__source", "bead_present", "volume_ml"], observed=True)
    out = g.agg(
        length_mean=("length", "mean"),
        length_sd=("length", "std"),
        length_n=("length", "count"),
        area_mean=("area", "mean"),
        area_sd=("area", "std"),
        area_n=("area", "count"),
    ).reset_index()
    return out


def plot_strain_facets(df_summary: pd.DataFrame, strain: str, outpath: Path) -> None:
    df = df_summary[df_summary["strain"] == strain].copy()
    if df.empty:
        print(f"No summary rows for strain {strain}, skipping plot")
        return

    df["bead_label"] = df["bead_present"].apply(lambda b: "no bead" if not b else "bead")
    sources = ["96-well", "24-well", "ttubes"]
    df["__source"] = pd.Categorical(df["__source"], categories=sources, ordered=True)
    df = df.sort_values(["__source", "bead_present", "volume_ml"])

    phenotypes = ["length", "area"]
    # Always produce a 2x2 panel grid: length, area, OD morning, OD afternoon.
    # If OD columns are absent, the OD panels will show a 'No OD data' note.
    # adjusted figure size: reduced ~30% from the previous 20x12 to give
    # panels more room than the original but avoid being overly large
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(14, 8.4))
    axes = axes.reshape(2, 2)
    # set parchment background for figure and axes (Arcadia palette)
    parchment = _arc_hex("parchment", "#F7F3EA")
    try:
        fig.patch.set_facecolor(parchment)
    except Exception:
        pass
    for ax in axes.flat:
        try:
            ax.set_facecolor(parchment)
        except Exception:
            pass
    # show only left and bottom spines (black axis lines); hide top/right
    for ax in axes.flat:
        try:
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)
            for spine in ("left", "bottom"):
                ax.spines[spine].set_visible(True)
                ax.spines[spine].set_color("black")
                ax.spines[spine].set_linewidth(1.0)
        except Exception:
            pass
    has_od = all(
        c in df.columns
        for c in ("OD_morning_mean", "OD_morning_stdev", "OD_afternoon_mean", "OD_afternoon_stdev")
    )

    bead_labels = ["no bead", "bead"]
    # explicit color mapping per (source, bead_label) so the two treatments per source are
    # visually similar but distinct. Example palette:
    # 96-well -> brown (no bead), orange (bead)
    # 24-well -> pink (no bead), red (bead)
    # ttubes  -> dark green (no bead), lime green (bead)
    color_map = {
        ("96-well", "no bead"): _arc_hex("canary", "#E5C272"),
        ("96-well", "bead"): _arc_hex("mustard", "#CF9202"),
        ("24-well", "no bead"): _arc_hex("vital", "#FF7A7A"),
        ("24-well", "bead"): _arc_hex("aegean", "#FF2A00"),
        ("ttubes", "no bead"): _arc_hex("tangerine", "#00C100"),
        ("ttubes", "bead"): _arc_hex("amber", "#006400"),
    }
    marker_map = {"no bead": "o", "bead": "s"}
    base_offsets = {"96-well": -0.12, "24-well": 0.0, "ttubes": 0.12}

    for j_col, phen in enumerate(phenotypes):
        ax = axes[0, j_col]
        for src in sources:
            for bead_label in bead_labels:
                dsrc = df[(df["__source"] == src) & (df["bead_label"] == bead_label)]
                if dsrc.empty:
                    continue
                x = dsrc["volume_ml"].astype(float).values
                y = dsrc[f"{phen}_mean"].values
                yerr = dsrc.get(f"{phen}_sd", pd.Series([0] * len(dsrc))).values
                order = np.argsort(np.asarray(x))
                x = np.asarray(x)[order]
                y = np.asarray(y)[order]
                yerr = np.asarray(yerr)[order]
                x_plot = x + base_offsets.get(src, 0.0)
                color = color_map.get((src, bead_label), None)
                label = f"{src} ({bead_label})"
                ax.errorbar(
                    x_plot,
                    y,
                    yerr=yerr,
                    fmt=marker_map.get(bead_label, "o"),
                    color=color,
                    capsize=0,
                    label=label,
                )
                ax.plot(x_plot, y, linestyle="-", color=color, alpha=0.8)

        ax.set_title(f"{strain} — {phen}")
        ax.set_xlabel("Volume (mL)")
        ax.set_xticks([1, 2, 3, 4, 5])
        ax.set_xlim(0.8, 5.2)
        ax.set_ylabel(f"{phen} (mean)")
        # disable background grid lines for a cleaner parchment-style figure
        ax.grid(False)

        # OD morning (lower-left)
        ax_od_m = axes[1, 0]
        if has_od:
            for src in sources:
                for bead_label in bead_labels:
                    dsrc = df[(df["__source"] == src) & (df["bead_label"] == bead_label)]
                    if dsrc.empty:
                        continue
                    x = dsrc["volume_ml"].astype(float).values
                    y = dsrc["OD_morning_mean"].values
                    yerr = dsrc["OD_morning_stdev"].values
                    order = np.argsort(np.asarray(x))
                    x = np.asarray(x)[order]
                    y = np.asarray(y)[order]
                    yerr = np.asarray(yerr)[order]
                    x_plot = x + base_offsets.get(src, 0.0)
                    color = color_map.get((src, bead_label), None)
                    ax_od_m.errorbar(
                        x_plot,
                        y,
                        yerr=yerr,
                        fmt=marker_map.get(bead_label, "o"),
                        color=color,
                        capsize=0,
                    )
                    ax_od_m.plot(x_plot, y, linestyle="-", color=color, alpha=0.8)
            ax_od_m.set_title(f"{strain} — OD morning")
            ax_od_m.set_xlabel("Volume (mL)")
            ax_od_m.set_xticks([1, 2, 3, 4, 5])
            ax_od_m.set_xlim(0.8, 5.2)
            ax_od_m.set_ylabel("OD (mean)")
            ax_od_m.grid(False)
        else:
            ax_od_m.text(
                0.5,
                0.5,
                "No OD morning data",
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax_od_m.transAxes,
            )
            ax_od_m.set_axis_off()

        # OD afternoon (lower-right)
        ax_od_a = axes[1, 1]
        if has_od:
            for src in sources:
                for bead_label in bead_labels:
                    dsrc = df[(df["__source"] == src) & (df["bead_label"] == bead_label)]
                    if dsrc.empty:
                        continue
                    x = dsrc["volume_ml"].astype(float).values
                    y = dsrc["OD_afternoon_mean"].values
                    yerr = dsrc["OD_afternoon_stdev"].values
                    order = np.argsort(np.asarray(x))
                    x = np.asarray(x)[order]
                    y = np.asarray(y)[order]
                    yerr = np.asarray(yerr)[order]
                    x_plot = x + base_offsets.get(src, 0.0)
                    color = color_map.get((src, bead_label), None)
                    ax_od_a.errorbar(
                        x_plot,
                        y,
                        yerr=yerr,
                        fmt=marker_map.get(bead_label, "o"),
                        color=color,
                        capsize=0,
                    )
                    ax_od_a.plot(x_plot, y, linestyle="-", color=color, alpha=0.8)
            ax_od_a.set_title(f"{strain} — OD afternoon")
            ax_od_a.set_xlabel("Volume (mL)")
            ax_od_a.set_xticks([1, 2, 3, 4, 5])
            ax_od_a.set_xlim(0.8, 5.2)
            ax_od_a.set_ylabel("OD (mean)")
            ax_od_a.grid(False)
        else:
            ax_od_a.text(
                0.5,
                0.5,
                "No OD afternoon data",
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax_od_a.transAxes,
            )
            ax_od_a.set_axis_off()

    import matplotlib.lines as mlines

    # create legend entries with consistent colors by source and markers by bead label
    handles = [
        mlines.Line2D(
            [],
            [],
            color=color_map.get((s, b)),
            marker=marker_map[b],
            linestyle="None",
            label=f"{s} ({b})",
        )
        for s in sources
        for b in bead_labels
    ]
    # place legend outside on the right and leave space in layout
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.suptitle(f"{strain} — phenotypes by volume (overlayed bead/no-bead)")
    outpath.parent.mkdir(parents=True, exist_ok=True)
    # tighten layout and reserve space on the right for the legend so titles/labels do not overlap
    try:
        fig.tight_layout(rect=(0, 0, 0.78, 0.95))
    except Exception:
        # fallback to a reasonable subplots_adjust if tight_layout fails
        fig.subplots_adjust(right=0.78, hspace=0.35, wspace=0.25, top=0.92)
    # save PNG and SVG (use bbox_inches='tight' so the legend outside the axes is included)
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    try:
        svg_path = outpath.with_suffix(".svg")
        fig.savefig(svg_path, bbox_inches="tight")
    except Exception:
        pass
    plt.close(fig)


# The metric-specific 'by_volume' plotting function was removed per user request.
# If you later want to restore metric-specific plots, re-add a `plot_metric` function
# or call `plot_metric` from a separate script.


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Aggregate multiple processed folders and plot mean length/area by volume"
    )
    parser.add_argument(
        "--paths",
        "-p",
        nargs=3,
        required=True,
        help="Three base folders: 96-well, 24-well, ttubes (in that order)",
    )
    parser.add_argument(
        "--map-96-volume",
        type=float,
        default=1.0,
        help=(
            "Volume (mL) to assign to 96-well 4.5 mm bead samples when mapping to the 1..5 mL "
            "x-axis",
        ),
    )
    parser.add_argument("--out-dir", "-o", default=Path.cwd(), type=Path)
    parser.add_argument(
        "--plots",
        action="store_true",
        help="Generate per-strain facetted plots (length, area, OD) into out-dir",
    )
    args = parser.parse_args(argv)

    labels = ["96-well", "24-well", "ttubes"]
    parts = []
    for path, label in zip(args.paths, labels, strict=False):
        p = Path(path)
        df = load_source(p, label)
        df = normalize_columns(df)
        df["__source_label"] = label
        parts.append(df)

    big = pd.concat(parts, ignore_index=True)

    # For 96-well, assign the mapped --map-96-volume to 4.5mm bead samples and explicit 'no bead'
    # rows but do NOT drop other 96-well rows (they may contain DIC measurements for other bead
    # sizes).
    mask_96 = big["__source"] == "96-well"
    if mask_96.any() and "treatment" in big.columns:
        treatment = big["treatment"].astype(str)
        # match '4.5 mm bead' (allow 4.5 or 45 with optional dot/space variants)
        cond_bead45 = treatment.str.contains(r"4\.?5\s*mm|4\.5", na=False, case=False)
        cond_nobead = treatment.str.contains(r"\bno\s*bead\b", na=False, case=False)
        # assign mapped volume to the matched 96-well rows (both bead and no-bead), keep all others
        big.loc[mask_96 & cond_bead45, "volume_ml"] = float(args.map_96_volume)
        big.loc[mask_96 & cond_nobead, "volume_ml"] = float(args.map_96_volume)

    summary = summarize(big)

    outdir = args.out_dir
    # Attempt to find an existing expanded summary that includes OD morning/afternoon
    # columns (for example: 'aggregate_summary_by_volume_plus my data.csv') in the
    # output directory and merge OD columns into our summary if present.
    try:
        # look for any summary file starting with the base name in the outdir
        candidates = list(outdir.glob("aggregate_summary_by_volume*.csv"))
        for c in candidates:
            try:
                ext = pd.read_csv(c)
            except Exception:
                continue
            od_cols = {
                "OD_morning_mean",
                "OD_morning_stdev",
                "OD_afternoon_mean",
                "OD_afternoon_stdev",
            }
            if od_cols.issubset(set(ext.columns)):
                # merge ext's OD columns into summary on the grouping keys
                merge_keys = ["strain", "__source", "bead_present", "volume_ml"]
                # keep existing summary columns and bring in OD columns from ext
                ext_small = ext[merge_keys + list(od_cols)].copy()
                # ensure consistent types for merge keys
                for k in merge_keys:
                    if k in ext_small.columns and k in summary.columns:
                        try:
                            ext_small[k] = ext_small[k].astype(summary[k].dtype)
                        except Exception:
                            pass
                summary = summary.merge(ext_small, on=merge_keys, how="left")
                break
    except Exception:
        # non-fatal if merging fails; plotting will simply show empty OD panels
        pass
    # write the summary CSV
    summary.to_csv(outdir / "aggregate_summary_by_volume.csv", index=False)

    print(f"Wrote summary CSV to {outdir}")

    # optionally produce per-strain facetted plots (length, area, OD morning/afternoon)
    if args.plots:
        strains = summary["strain"].dropna().unique().tolist()
        preferred = ["SP286", "dea2^"]
        for s in preferred:
            if s in strains:
                safe_name = re.sub(r"[^A-Za-z0-9_-]", "_", s)
                outpath = outdir / f"aggregate_{safe_name}_by_bead_and_metric.png"
                plot_strain_facets(summary, s, outpath)
        # also produce for any remaining strains (non-preferred)
        for s in strains:
            if s in preferred or s is None:
                continue
            safe_name = re.sub(r"[^A-Za-z0-9_-]", "_", s)
            outpath = outdir / f"aggregate_{safe_name}_by_bead_and_metric.png"
            plot_strain_facets(summary, s, outpath)


if __name__ == "__main__":
    raise SystemExit(main())
