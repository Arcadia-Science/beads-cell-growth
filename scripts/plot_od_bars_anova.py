from pathlib import Path

import arcadia_pycolor.colors as colors
import arcadia_pycolor.mpl as apc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scikit_posthocs as sp
from scipy.stats import f_oneway

font_dirpath = Path.home() / "Library/Fonts"
apc.setup(font_dirpath=str(font_dirpath))


input_path = "/Users/roman/Repositories/2026-pombe-beads/Data/All_ODs.xlsx"
df = pd.read_excel(input_path)


# Map beads: 0.0 (FALSE) stays 0.0, 1.0 (TRUE) becomes 4.5, others unchanged
def map_beads(val):
    # Only map boolean True or string 'TRUE' to 4.5, not numeric 1.0
    if isinstance(val, bool) and val is True:
        return 4.5
    if isinstance(val, str) and val.strip().upper() == "TRUE":
        return 4.5
    if isinstance(val, bool) and val is False:
        return 0.0
    if isinstance(val, str) and val.strip().upper() == "FALSE":
        return 0.0
    return val


df["beads"] = df["beads"].apply(map_beads)
output_dir = Path(input_path).parent

group_cols = ["Experiment", "strain", "time", "volume", "beads"]
summ = df.groupby(group_cols)["OD"].agg(["mean", "std", "count"]).reset_index()


def pval_to_stars(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return ""


def plot_bar_with_anova(df, summary, experiment, time, output_dir):
    plt.close("all")
    parchment = str(colors.parchment)
    fig, ax = plt.subplots(figsize=(8, 6), facecolor=parchment)
    ax.set_facecolor(parchment)
    fig.patch.set_facecolor(parchment)
    sub = summary[(summary["Experiment"] == experiment) & (summary["time"] == time)]
    raw = df[(df["Experiment"] == experiment) & (df["time"] == time)]
    if sub.empty:
        print(f"No data for {experiment} {time}")
        return
    if experiment == "96-well":
        # Always include these bead values on x-axis
        bead_values = [0.0, 1.0, 1.5, 3.0, 4.5]
        bead_map = {0.0: "0", 1.0: "1", 1.5: "1.5", 3.0: "3", 4.5: "4.5"}
        means = []
        sems = []
        for b in bead_values:
            vals = sub[sub["beads"] == b]
            means.append(vals["mean"].values[0] if not vals.empty else float("nan"))
            # Calculate SEM: std / sqrt(n)
            if not vals.empty and vals["count"].values[0] > 0:
                sems.append(vals["std"].values[0] / np.sqrt(vals["count"].values[0]))
            else:
                sems.append(float("nan"))
        x = range(len(bead_values))
        colors_list = [
            str(colors.vital) if b == 0 or b == 0.0 else str(colors.aegean) for b in bead_values
        ]
        ax.bar(x, means, yerr=sems, color=colors_list, edgecolor=str(colors.charcoal), capsize=4)
        ax.set_xticks(x)
        ax.set_xticklabels([bead_map.get(b, str(b)) for b in bead_values])
        ax.set_xlabel("Bead Amount")
        # ANOVA: overall test
        groups = [raw[raw["beads"] == b]["OD"].values for b in bead_values]
        if all(len(g) > 1 for g in groups) and len(groups) > 1:
            _, pval = f_oneway(*groups)
            # Post-hoc: Dunnett's test (vs. control beads==0)
            dunnett_data = raw[["OD", "beads"]].dropna()
            if dunnett_data["beads"].nunique() > 1:
                dunnett_data["beads_str"] = dunnett_data["beads"].astype(str)
                dunnett = sp.posthoc_dunnett(
                    dunnett_data, val_col="OD", group_col="beads_str", control="0.0"
                )
                for i, b in enumerate(bead_values):
                    if b == 0 or b == 0.0:
                        continue
                    pval = dunnett.loc[str(b), "0.0"] if str(b) in dunnett.index else None
                    if pval is not None:
                        y_max = max(means[0], means[i])
                        y = y_max + 0.15 + 0.1 * i
                        # Only show the star, no lines at all
                        ax.text(
                            (0 + i) / 2,
                            y + 0.01,
                            pval_to_stars(pval),
                            ha="center",
                            va="bottom",
                            color="black",
                            fontsize=16,
                        )
    else:
        volumes = sorted(sub["volume"].unique())
        beads_vals = sorted(sub["beads"].unique())
        width = 0.35
        for i, beads in enumerate(beads_vals):
            vals = sub[sub["beads"] == beads]
            means = []
            sems = []
            for v in volumes:
                vvals = vals[vals["volume"] == v]
                means.append(vvals["mean"].values[0] if not vvals.empty else float("nan"))
                if not vvals.empty and vvals["count"].values[0] > 0:
                    sems.append(vvals["std"].values[0] / np.sqrt(vvals["count"].values[0]))
                else:
                    sems.append(float("nan"))
            offset = (i - 0.5) * width
            color = str(colors.vital) if beads == 0 or beads == 0.0 else str(colors.aegean)
            ax.bar(
                [v + offset for v in range(len(volumes))],
                means,
                width,
                yerr=sems,
                color=color,
                edgecolor=str(colors.charcoal),
                capsize=4,
            )
        ax.set_xticks(range(len(volumes)))
        ax.set_xticklabels([str(v) for v in volumes])
        ax.set_xlabel("Volume")
        # ANOVA and Dunnett's test for each volume
        for idx, v in enumerate(volumes):
            raw_v = raw[raw["volume"] == v]
            beads_in_v = sorted(raw_v["beads"].unique())
            groups = [raw_v[raw_v["beads"] == b]["OD"].values for b in beads_in_v]
            if all(len(g) > 1 for g in groups) and len(groups) > 1:
                _, pval = f_oneway(*groups)
                dunnett_data = raw_v[["OD", "beads"]].dropna()
                if dunnett_data["beads"].nunique() > 1:
                    dunnett_data["beads_str"] = dunnett_data["beads"].astype(str)
                    dunnett = sp.posthoc_dunnett(
                        dunnett_data, val_col="OD", group_col="beads_str", control="0.0"
                    )
                    for i, b in enumerate(beads_in_v):
                        if b == 0 or b == 0.0:
                            continue
                        pval = dunnett.loc[str(b), "0.0"] if str(b) in dunnett.index else None
                        if pval is not None:
                            mean1 = sub[(sub["volume"] == v) & (sub["beads"] == 0)]["mean"].values[
                                0
                            ]
                            mean2 = sub[(sub["volume"] == v) & (sub["beads"] == b)]["mean"].values[
                                0
                            ]
                            y_max = max(mean1, mean2)
                            y = y_max + 0.15 + 0.1 * i
                            pos1 = (
                                idx + (beads_in_v.index(0) - 0.5) * width
                                if 0 in beads_in_v
                                else idx
                            )
                            pos2 = idx + (i - 0.5) * width
                            # Only show the star, no lines at all
                            ax.text(
                                (pos1 + pos2) / 2,
                                y + 0.01,
                                pval_to_stars(pval),
                                ha="center",
                                va="bottom",
                                color="black",
                                fontsize=16,
                            )
    ax.set_ylabel("OD")
    ax.set_title(f"{experiment} - {time}")
    # ax.legend()  # Legend removed as requested
    fig.tight_layout()
    out_path_png = output_dir / f"od_bar_{experiment}_{time}.png"
    out_path_svg = output_dir / f"od_bar_{experiment}_{time}.svg"
    plt.savefig(out_path_png, facecolor=parchment, edgecolor=parchment, bbox_inches="tight")
    plt.savefig(out_path_svg, facecolor=parchment, edgecolor=parchment, bbox_inches="tight")
    print(f"Saved: {out_path_png}")
    print(f"Saved: {out_path_svg}")
    plt.close(fig)


for experiment in summ["Experiment"].unique():
    for time in summ["time"].unique():
        if str(time).lower() != "afternoon":
            continue
        plot_bar_with_anova(df, summ, experiment, time, output_dir)
