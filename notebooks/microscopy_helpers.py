from pathlib import Path

import arcadia_pycolor as apc
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

METRIC_LABELS = {
    "area": "area (\u03bcm\u00b2)",
    "length": "length (\u03bcm)",
}

BEAD_LABELS = {False: "no bead", True: "bead"}


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


def plot_metric_by_bead(df, metric, output_path):
    """Bar plot of mean `metric` by bead size with Holm-Bonferroni-corrected t-tests vs. control."""
    bead_values = df["bead_size_mm"].unique()

    # Aggregate to per-image means to avoid pseudoreplication:
    # cells within the same image are not independent observations.
    image_means = df.groupby(["bead_size_mm", "source_file"])[metric].mean().reset_index()
    means = image_means.groupby("bead_size_mm")[metric].mean()
    stds = image_means.groupby("bead_size_mm")[metric].std()
    counts = image_means.groupby("bead_size_mm")[metric].count()
    cis = 1.96 * stds / np.sqrt(counts)

    x = np.arange(len(means))
    fig, ax = plt.subplots(figsize=(10, 6), facecolor=apc.parchment)
    ax.bar(x, means, yerr=cis, capsize=6)
    ax.set_xticks(x)
    ax.set_xticklabels([str(b) for b in means.index])
    ax.set_xlabel("Bead size (mm)")
    ax.set_ylabel(METRIC_LABELS.get(metric, metric))
    ax.set_title(f"Mean cell {metric} by bead treatment (96-well)")

    control_vals = image_means[image_means["bead_size_mm"] == 0.0][metric].values
    treatment_beads = sorted(b for b in bead_values if b != 0.0)
    raw_pvals = []
    for b in treatment_beads:
        treat_vals = image_means[image_means["bead_size_mm"] == b][metric].values
        if len(control_vals) > 1 and len(treat_vals) > 1:
            _, p = ttest_ind(control_vals, treat_vals, equal_var=False)
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
                bar_idx, means[b] + cis[b] + y_offset, stars, ha="center", va="bottom", fontsize=14
            )

    ax.set_ylim(bottom=0, top=(means + cis).max() + y_offset * 6)
    apc.mpl.style_plot(ax)
    plt.tight_layout()
    fig.savefig(Path(output_path), bbox_inches="tight", facecolor=apc.parchment)
    return fig


def plot_metric_by_volume(df, metric, title, output_path):
    """Grouped bar plot of mean `metric` by volume with bead vs. no-bead,
    annotated with Holm-Bonferroni-corrected Welch's t-tests."""
    volumes = sorted(df["volume_ml"].unique())
    conditions = [False, True]

    # Aggregate to per-image means to avoid pseudoreplication:
    # cells within the same image are not independent observations.
    image_means = (
        df.groupby(["volume_ml", "bead_present", "source_file"])[metric].mean().reset_index()
    )
    means = image_means.groupby(["volume_ml", "bead_present"])[metric].mean().unstack()
    stds = image_means.groupby(["volume_ml", "bead_present"])[metric].std().unstack()
    counts = image_means.groupby(["volume_ml", "bead_present"])[metric].count().unstack()
    cis = 1.96 * stds / np.sqrt(counts)

    x = np.arange(len(volumes))
    bar_w = 0.8 / len(conditions)

    fig, ax = plt.subplots(figsize=(10, 6), facecolor=apc.parchment)
    for i, cond in enumerate(conditions):
        ax.bar(
            x + (i - (len(conditions) - 1) / 2) * bar_w,
            means[cond],
            width=bar_w,
            label=BEAD_LABELS[cond],
            yerr=cis[cond],
            capsize=6,
        )

    ax.set_xticks(x)
    ax.set_xticklabels([str(v) for v in volumes])
    ax.set_xlabel("Volume (mL)")
    ax.set_ylabel(METRIC_LABELS.get(metric, metric))
    ax.set_title(f"Mean cell {metric} by volume ({title})")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

    raw_pvals = []
    for v in volumes:
        ctrl = image_means[(image_means["volume_ml"] == v) & (~image_means["bead_present"])][
            metric
        ].values
        treat = image_means[(image_means["volume_ml"] == v) & (image_means["bead_present"])][
            metric
        ].values
        if len(ctrl) > 1 and len(treat) > 1:
            _, p = ttest_ind(ctrl, treat, equal_var=False)
        else:
            p = np.nan
        raw_pvals.append(p)

    adj_pvals = correct_pvalues(raw_pvals)
    y_offset = means.max().max() * 0.03
    for v, adj_p in zip(volumes, adj_pvals, strict=True):
        stars = pval_to_stars(adj_p)
        if stars:
            vi = volumes.index(v)
            ctrl_pos = x[vi] + (conditions.index(False) - (len(conditions) - 1) / 2) * bar_w
            treat_pos = x[vi] + (conditions.index(True) - (len(conditions) - 1) / 2) * bar_w
            y_max = max(
                means.loc[v, False] + cis.loc[v, False],
                means.loc[v, True] + cis.loc[v, True],
            )
            ax.text(
                (ctrl_pos + treat_pos) / 2,
                y_max + y_offset,
                stars,
                ha="center",
                va="bottom",
                fontsize=14,
            )

    y_max = means.max().max() + cis.max().max() + y_offset * 6
    ax.set_ylim(0, y_max)

    apc.mpl.style_plot(ax)
    plt.tight_layout()

    fig.savefig(Path(output_path), bbox_inches="tight", facecolor=apc.parchment)

    return fig


def plot_metric_cross_experiment(df_agg, metric, output_path):
    """Line plot of mean `metric` vs. volume across experiments, split by bead presence."""
    fig, ax = plt.subplots(figsize=(10, 6), facecolor=apc.parchment)

    for (exp, bead), grp in df_agg.groupby(["experiment", "bead_present"]):
        # Aggregate to per-image means to avoid pseudoreplication:
        # cells within the same image are not independent observations.
        image_means = grp.groupby(["volume_ml", "source_file"])[metric].mean().reset_index()
        stats = image_means.groupby("volume_ml")[metric].agg(["mean", "std", "count"])
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
    ax.set_ylabel(METRIC_LABELS.get(metric, metric))
    ax.set_title(f"Mean cell {metric} vs. volume by experiment and bead presence")

    apc.mpl.style_plot(ax)
    plt.tight_layout()
    fig.savefig(Path(output_path), bbox_inches="tight", facecolor=apc.parchment)
    return fig
