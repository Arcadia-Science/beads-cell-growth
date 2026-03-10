"""Create OD heatmaps from Excel files in Data/.

For each of the two Excel files ("Baseline_ODs.xlsx" and "Supplement_ODs.xlsx") this
script builds a grid of small heatmaps facetted by Experiment (rows) and
genotype (columns). Each small heatmap has two rows: morning (row 0) and
afternoon (row 1). Columns are ordered in the same sequence as they appear in
the Excel sheet for that Experiment+genotype (samples are grouped by the other
sample-identifying columns).

Color map: 'magma' (matplotlib). Outputs saved to Data_analysis/heatmaps.

Usage: python make_od_heatmaps.py

"""

from pathlib import Path

import arcadia_pycolor as ap
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

rc = getattr(ap.style_defaults, "ARCADIA_MATPLOTLIB_RC_PARAMS", None)
if rc:
    plt.rcParams.update(rc)
plt.rcParams["savefig.transparent"] = False
plt.rcParams["image.cmap"] = "magma"
plt.rcParams["font.size"] = ap.style_defaults.BASE_FONT_SIZE

DATA_DIR = Path("data")
OUT_DIR = DATA_DIR / "heatmaps"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FILES = ["Supplement_ODs.xlsx", "Baseline_ODs.xlsx"]


def _read_first_sheet(path: Path) -> pd.DataFrame:
    return pd.read_excel(path, sheet_name=0)


def _ordered_unique(seq):
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def build_matrices(df: pd.DataFrame):
    # normalize column names for robust access
    colmap = {c: c.strip() for c in df.columns}
    df = df.rename(columns=colmap)
    # required columns: Experiment, genotype, time, OD
    cols = [c for c in df.columns]
    key_cols = [c for c in cols if c.lower() not in ("experiment", "genotype", "time", "od")]

    experiments = _ordered_unique(df["Experiment"].astype(str).tolist())
    genotypes = _ordered_unique(df["genotype"].astype(str).tolist())

    matrices = {}  # (exp,geno) -> (2 x n_samples array, col_labels)
    all_vals = []
    for exp in experiments:
        for geno in genotypes:
            sub = df[(df["Experiment"].astype(str) == exp) & (df["genotype"].astype(str) == geno)]
            if sub.empty:
                continue
            # build sample keys preserving order of appearance
            sample_keys = []
            for _, row in sub.iterrows():
                key = " | ".join(str(row[c]) for c in key_cols)
                if key not in sample_keys:
                    sample_keys.append(key)
            n = len(sample_keys)
            mat = np.full((2, n), np.nan)
            for _, row in sub.iterrows():
                key = " | ".join(str(row[c]) for c in key_cols)
                try:
                    col = sample_keys.index(key)
                except ValueError:
                    continue
                t = str(row["time"]).lower() if not pd.isna(row["time"]) else ""
                if "morning" in t:
                    r = 0
                elif "afternoon" in t:
                    r = 1
                else:
                    # unknown time; skip
                    continue
                try:
                    val = float(row["OD"])
                except Exception:
                    val = np.nan
                mat[r, col] = val
                if not np.isnan(val):
                    all_vals.append(val)
            matrices[(exp, geno)] = (mat, sample_keys)
    return experiments, genotypes, matrices, all_vals


def plot_file(path: Path):
    df = _read_first_sheet(path)
    experiments, genotypes, matrices, all_vals = build_matrices(df)
    if len(experiments) == 0 or len(genotypes) == 0:
        print("No experiments/genotypes found in", path)
        return

    nrow = len(experiments)
    ncol = len(genotypes)
    # choose figure size so each heatmap 'cell' is a perfect square.
    # Use a cell size in inches and set figsize = (ncol*cell, nrow*cell).
    cell_size_in = 0.5  # inches per cell (tweakable)
    fig_w = max(4, ncol * cell_size_in)
    fig_h = max(2, nrow * cell_size_in)
    # scale the entire figure by SCALE (1.0 = original). Set to 3.6 to
    # make the figure 360% larger in both x and y (previously 300% + 20%).
    SCALE = 3.6
    fig_w *= SCALE
    fig_h *= SCALE
    # add extra horizontal space (in inches) to accommodate per-facet
    # colorbars so they don't overlap neighboring facets
    per_cbar_in = cell_size_in * 0.08
    extra_cbar_total = ncol * per_cbar_in
    fig_w += extra_cbar_total
    fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(fig_w, fig_h))
    if nrow == 1 and ncol == 1:
        axes = np.array([[axes]])
    elif nrow == 1:
        axes = axes.reshape(1, -1)
    elif ncol == 1:
        axes = axes.reshape(-1, 1)

    # compute vmin/vmax from all values; if empty fallback
    if len(all_vals) > 0:
        vmin = np.nanmin(all_vals)
        vmax = np.nanmax(all_vals)
    else:
        vmin, vmax = 0.0, 1.0

    # NOTE: we keep the global vmin/vmax computed above and by default use
    # the same scale across all strains/genres so different facets are
    # comparable. (If you want per-experiment scales instead, we can
    # re-enable that behavior.)

    # prefer the Arcadia 'magma' gradient converted to a Matplotlib Colormap
    try:
        cmap_obj = ap.gradients.magma.to_mpl_cmap()
    except Exception:
        cmap_obj = plt.get_cmap("magma")

    # reverse the colormap so the gradient is flipped (user request)
    try:
        cmap_obj = cmap_obj.reversed()
    except Exception:
        # fallback: if the colormap has a name, attempt to get the '_r' variant
        try:
            name = getattr(cmap_obj, "name", None)
            if name:
                cmap_obj = mpl.colormaps.get(name + "_r")
            else:
                cmap_obj = mpl.colormaps.get("magma_r")
        except Exception:
            cmap_obj = mpl.colormaps.get("magma_r")

    parchment = str(getattr(ap, "parchment", "#F7F3EA"))
    fig.patch.set_facecolor(parchment)
    # ensure legend font matches axes
    try:
        base_fs = plt.rcParams.get("font.size", 10)
        plt.rcParams["legend.fontsize"] = base_fs
    except Exception:
        pass
    # adjust subplot spacing to avoid title/label overlap; give more room on
    # the right for per-facet colorbars
    try:
        fig.subplots_adjust(top=0.92, bottom=0.08, left=0.07, right=0.95, hspace=0.6, wspace=0.6)
    except Exception:
        pass

    for i, exp in enumerate(experiments):
        for j, geno in enumerate(genotypes):
            ax = axes[i, j]
            key = (exp, geno)
            if key not in matrices:
                ax.set_visible(False)
                continue
            mat, col_labels = matrices[key]
            # plot heatmap using imshow
            # set axes facecolor to parchment for consistent background
            try:
                ax.set_facecolor(parchment)
            except Exception:
                pass
            # use the same global vmin/vmax across all facets so strains are comparable
            ax.imshow(
                mat, aspect="equal", cmap=cmap_obj, vmin=vmin, vmax=vmax, interpolation="nearest"
            )
            # y ticks
            ax.set_yticks([0, 1])
            ax.set_yticklabels(["morning", "afternoon"])
            # set axis tick label sizes to base font so they match legends
            base_fs = plt.rcParams.get("font.size", 10)
            ax.tick_params(axis="both", labelsize=base_fs)
            # x ticks = each sample
            ncols = mat.shape[1]
            ax.set_xticks(np.arange(ncols))
            # set xtick labels as short versions of col_labels to avoid crowding
            short_labels = [
                label if len(label) <= 12 else (label[:9] + "...") for label in col_labels
            ]
            ax.set_xticklabels(short_labels, rotation=90, fontsize=8)
            # label row/column facets
            if j == 0:
                # leftmost column: show experiment name as y-label on left
                ax.set_ylabel(exp, fontsize=base_fs)
            if i == 0:
                # top row: put genotype as title
                ax.set_title(geno, fontsize=base_fs)
            # spine styling: only left/bottom visible
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)
            for spine in ("left", "bottom"):
                ax.spines[spine].set_visible(True)
                ax.spines[spine].set_color("black")
                ax.spines[spine].set_linewidth(1.0)
            # no per-facet colorbar; we'll draw a single global colorbar below
            pass
    # draw one colorbar per group: 'ttubes', '96', '24'
    groups = ["ttubes", "96", "24"]

    # helper to detect group membership from strings
    def detect_group(exp, geno, col_labels):
        s = " ".join([str(exp), str(geno)] + [str(x) for x in col_labels]).lower()
        if "ttube" in s or "ttubes" in s:
            return "ttubes"
        if "96" in s:
            return "96"
        if "24" in s:
            return "24"
        # fallback: try genotype
        g = str(geno).lower()
        if "ttube" in g:
            return "ttubes"
        if "96" in g:
            return "96"
        if "24" in g:
            return "24"
        return "96"

    # collect axes positions for each group
    group_axes = {g: [] for g in groups}
    for i, exp in enumerate(experiments):
        for j, geno in enumerate(genotypes):
            key = (exp, geno)
            ax = axes[i, j]
            if key not in matrices:
                continue
            _, col_labels = matrices[key]
            grp = detect_group(exp, geno, col_labels)
            group_axes.setdefault(grp, []).append(ax)

    # compute per-group vmin/vmax from global data but limited to group members
    group_ranges = {}
    for g in groups:
        vals = []
        for ax in group_axes.get(g, []):
            # find which axes index this is
            # search in axes grid
            found = False
            for ii in range(axes.shape[0]):
                for jj in range(axes.shape[1]):
                    if axes[ii, jj] is ax:
                        key = (experiments[ii], genotypes[jj])
                        if key in matrices:
                            mat, _ = matrices[key]
                            f = mat[np.isfinite(mat)]
                            if f.size > 0:
                                vals.append(f)
                        found = True
                        break
                if found:
                    break
        if vals:
            allv = np.concatenate(vals)
            group_ranges[g] = (float(np.nanmin(allv)), float(np.nanmax(allv)))
        else:
            group_ranges[g] = (vmin, vmax)

    # draw colorbar for each group positioned to the right of the rightmost
    # axis in that group and spanning the vertical span of its axes
    for g in groups:
        axes_list = group_axes.get(g, [])
        if not axes_list:
            continue
        y0 = min(a.get_position().y0 for a in axes_list)
        y1 = max(a.get_position().y1 for a in axes_list)
        x1 = max(a.get_position().x1 for a in axes_list)
        pad = 0.01
        cbar_w = 0.02
        try:
            cax = fig.add_axes((x1 + pad, y0, cbar_w, y1 - y0))
            gvmin, gvmax = group_ranges[g]
            gradient = np.linspace(gvmin, gvmax, 256).reshape(-1, 1)
            cax.imshow(
                gradient,
                aspect="auto",
                cmap=cmap_obj,
                origin="lower",
                extent=(0, 1, gvmin, gvmax),
                interpolation="nearest",
            )
            cax.set_xticks([])
            ticks = np.linspace(gvmin, gvmax, num=5)
            cax.set_yticks(ticks)
            cax.yaxis.tick_right()
            cax.yaxis.set_label_position("right")
            try:
                cax.set_facecolor(parchment)
            except Exception:
                pass
            base_fs = plt.rcParams.get("font.size", 10)
            cax.tick_params(axis="y", labelsize=base_fs)
            # add a small label identifying the group
            cax.set_title(g, fontsize=base_fs, pad=2)
        except Exception:
            pass

    out_png = OUT_DIR / (path.stem.replace(" ", "_") + "_heatmap.png")
    out_svg = OUT_DIR / (path.stem.replace(" ", "_") + "_heatmap.svg")
    fig.tight_layout(rect=(0, 0, 0.9, 0.95))
    # save PNG and SVG with explicit non-transparent background matching the
    # figure patch (parchment)
    try:
        fig.savefig(out_png, dpi=300, transparent=False, facecolor=fig.get_facecolor())
    except Exception:
        fig.savefig(out_png, dpi=300)
    try:
        fig.savefig(out_svg, transparent=False, facecolor=fig.get_facecolor())
    except Exception:
        pass
    plt.close(fig)
    print("Wrote", out_png)


if __name__ == "__main__":
    for fname in FILES:
        p = DATA_DIR / fname
        if p.exists():
            plot_file(p)
        else:
            print("Missing file", p)
