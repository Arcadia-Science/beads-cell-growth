from pathlib import Path

import arcadia_pycolor as apc
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

SUPP_ORDER = ["none", "D", "V", "A", "D+A", "V+A", "D+V", "D+V+A"]


def pval_to_stars(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return ""


def correct_pvalues(pvals):
    """Apply Holm-Bonferroni correction, tolerating NaNs."""
    valid_mask = [not np.isnan(p) for p in pvals]
    if not any(valid_mask):
        return list(pvals)
    valid_pvals = [p for p, v in zip(pvals, valid_mask, strict=True) if v]
    _, corrected, _, _ = multipletests(valid_pvals, method="holm")
    corrected_iter = iter(corrected)
    return [next(corrected_iter) if v else np.nan for v in valid_mask]


def plot_od_by_bead(df, title, output_path):
    """Bar plot of mean OD by bead size with Holm-Bonferroni-corrected t-tests vs. control."""
    bead_values = sorted(df["beads"].unique())
    means = df.groupby("beads")["OD"].mean()
    stds = df.groupby("beads")["OD"].std()
    counts = df.groupby("beads")["OD"].count()
    cis = 1.96 * stds / np.sqrt(counts)

    x = np.arange(len(means))
    fig, ax = plt.subplots(figsize=(8, 6), facecolor=apc.parchment)
    ax.bar(x, means, yerr=cis, capsize=6)
    ax.set_xticks(x)
    ax.set_xticklabels([str(b) for b in means.index])
    ax.set_xlabel("Bead size (mm)")
    ax.set_ylabel("OD")
    ax.set_title(title)

    control_od = df[df["beads"] == 0.0]["OD"].values
    treatment_beads = [b for b in bead_values if b != 0.0]
    raw_pvals = []
    for b in treatment_beads:
        treat_od = df[df["beads"] == b]["OD"].values
        if len(control_od) > 1 and len(treat_od) > 1:
            _, p = ttest_ind(control_od, treat_od, equal_var=False)
        else:
            p = np.nan
        raw_pvals.append(p)

    adj_pvals = correct_pvalues(raw_pvals)
    y_offset = means.max() * 0.05
    for b, adj_p in zip(treatment_beads, adj_pvals, strict=True):
        stars = pval_to_stars(adj_p)
        if stars:
            bar_idx = list(means.index).index(b)
            ax.text(
                bar_idx,
                means[b] + cis[b] + y_offset,
                stars,
                ha="center",
                va="bottom",
                fontsize=14,
            )

    apc.mpl.style_plot(ax)
    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(Path(output_path), bbox_inches="tight", facecolor=apc.parchment)
    return fig


def plot_od_by_volume(df, title, output_path):
    """Grouped bar plot of mean OD by volume and bead size with significance annotations."""
    volumes = sorted(df["volume"].unique())
    beads = sorted(df["beads"].unique())

    means = df.groupby(["volume", "beads"])["OD"].mean().unstack()
    stds = df.groupby(["volume", "beads"])["OD"].std().unstack()
    counts = df.groupby(["volume", "beads"])["OD"].count().unstack()
    cis = 1.96 * stds / np.sqrt(counts)

    x = np.arange(len(volumes))
    bar_w = 0.8 / len(beads)

    fig, ax = plt.subplots(figsize=(10, 6), facecolor=apc.parchment)
    for i, bead in enumerate(beads):
        ax.bar(
            x + (i - (len(beads) - 1) / 2) * bar_w,
            means[bead],
            width=bar_w,
            label=f"{bead}",
            yerr=cis[bead],
            capsize=6,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(volumes)
    ax.set_xlabel("Volume (mL)")
    ax.set_ylabel("OD")
    ax.set_title(title)
    ax.legend(title="Beads (mg/mL)", bbox_to_anchor=(1.05, 1), loc="upper left")

    treatment_beads = [b for b in beads if b != 0.0]
    comparisons, raw_pvals = [], []
    for v in volumes:
        control_od = df[(df["volume"] == v) & (df["beads"] == 0.0)]["OD"].values
        for b in treatment_beads:
            treat_od = df[(df["volume"] == v) & (df["beads"] == b)]["OD"].values
            if len(control_od) > 1 and len(treat_od) > 1:
                _, p = ttest_ind(control_od, treat_od, equal_var=False)
            else:
                p = np.nan
            comparisons.append((v, b))
            raw_pvals.append(p)

    adj_pvals = correct_pvalues(raw_pvals)
    y_offset = means.max().max() * 0.03
    for (v, b), adj_p in zip(comparisons, adj_pvals, strict=True):
        stars = pval_to_stars(adj_p)
        if stars:
            vi = volumes.index(v)
            ctrl_pos = x[vi] + (beads.index(0.0) - (len(beads) - 1) / 2) * bar_w
            treat_pos = x[vi] + (beads.index(b) - (len(beads) - 1) / 2) * bar_w
            y_max = max(means.loc[v, 0.0] + cis.loc[v, 0.0], means.loc[v, b] + cis.loc[v, b])
            ax.text(
                (ctrl_pos + treat_pos) / 2,
                y_max + y_offset,
                stars,
                ha="center",
                va="bottom",
                fontsize=14,
            )

    apc.mpl.style_plot(ax)
    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(Path(output_path), bbox_inches="tight", facecolor=apc.parchment)
    return fig


def plot_od_by_supplement(df, title, output_path, supp_order=None):
    """Bar plot of mean OD by supplement with Holm-Bonferroni-corrected t-tests vs. control."""
    if supp_order is None:
        supp_order = SUPP_ORDER

    means = df.groupby("supplement")["OD"].mean().reindex(supp_order)
    stds = df.groupby("supplement")["OD"].std().reindex(supp_order)
    counts = df.groupby("supplement")["OD"].count().reindex(supp_order)
    cis = 1.96 * stds / np.sqrt(counts)

    x = np.arange(len(supp_order))
    fig, ax = plt.subplots(figsize=(10, 6), facecolor=apc.parchment)
    ax.bar(x, means, yerr=cis, capsize=6)
    ax.set_xticks(x)
    ax.set_xticklabels(supp_order, rotation=45, ha="right")
    ax.set_xlabel("Supplement")
    ax.set_ylabel("OD")
    ax.set_title(title)

    control_od = df[df["supplement"] == "none"]["OD"].values
    treatment_supps = [s for s in supp_order if s != "none"]
    raw_pvals = []
    for s in treatment_supps:
        treat_od = df[df["supplement"] == s]["OD"].values
        if len(control_od) > 1 and len(treat_od) > 1:
            _, p = ttest_ind(control_od, treat_od, equal_var=False)
        else:
            p = np.nan
        raw_pvals.append(p)

    adj_pvals = correct_pvalues(raw_pvals)
    y_offset = means.max() * 0.05
    for s, adj_p in zip(treatment_supps, adj_pvals, strict=True):
        stars = pval_to_stars(adj_p)
        if stars:
            bar_idx = supp_order.index(s)
            ax.text(
                bar_idx,
                means[s] + cis[s] + y_offset,
                stars,
                ha="center",
                va="bottom",
                fontsize=14,
            )

    apc.mpl.style_plot(ax)
    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(Path(output_path), bbox_inches="tight", facecolor=apc.parchment)
    return fig


def plot_od_cross_experiment(df_agg, output_path=None):
    """Line plot of mean OD vs. volume across experiments, split by bead presence."""
    fig, ax = plt.subplots(figsize=(10, 6), facecolor=apc.parchment)

    for (exp, bead), grp in df_agg.groupby(["experiment", "bead_present"]):
        stats = grp.groupby("volume")["OD"].agg(["mean", "std", "count"])
        stats["ci"] = 1.96 * stats["std"] / np.sqrt(stats["count"])
        stats = stats.sort_index()
        label = f"{exp}, {'bead' if bead else 'no bead'}"
        if len(stats) == 1:
            ax.errorbar(
                stats.index,
                stats["mean"],
                yerr=stats["ci"],
                marker="o",
                capsize=4,
                label=label,
            )
        else:
            (line,) = ax.plot(stats.index, stats["mean"], marker="o", label=label)
            ax.fill_between(
                stats.index,
                stats["mean"] - stats["ci"],
                stats["mean"] + stats["ci"],
                alpha=0.2,
                color=line.get_color(),
            )

    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    ax.set_xlabel("Volume (mL)")
    ax.set_ylabel("OD")
    ax.set_title("Mean OD vs. volume by experiment and bead presence")

    apc.mpl.style_plot(ax)
    plt.tight_layout()
    if output_path is not None:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(Path(output_path), bbox_inches="tight", facecolor=apc.parchment)
    return fig
