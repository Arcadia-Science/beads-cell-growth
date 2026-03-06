from pathlib import Path
import re
import pandas as pd
import aggregate_multi_sources as ag


def main():
    # default choices (user accepted)
    paths = [
        Path("/Users/roman/Repositories/2026-pombe-beads/Data/2016-01-16/20260116_094944_372"),
        Path("/Users/roman/Repositories/2026-pombe-beads/Data/2016-01-22/20260122_111821_521"),
        Path("/Users/roman/Repositories/2026-pombe-beads/Data/2016-01-22/20260122_113404_129"),
    ]
    labels = ["96-well", "24-well", "ttubes"]
    parts = []
    for p, label in zip(paths, labels):
        print(f"Loading {p} as {label}")
        df = ag.load_source(Path(p), label)
        df = ag.normalize_columns(df)
        parts.append(df)

    big = pd.concat(parts, ignore_index=True)

    # apply same 96-well mapping as the main script (map_96_volume=1.0)
    mask_96 = big["__source"] == "96-well"
    if mask_96.any() and "treatment" in big.columns:
        treatment = big["treatment"].astype(str)
        cond_bead45 = treatment.str.contains(r"4\.?5\s*mm|4\.5", na=False, case=False)
        cond_nobead = treatment.str.contains(r"\bno\s*bead\b", na=False, case=False)
        big.loc[mask_96 & cond_bead45, 'volume_ml'] = 1.0
        big.loc[mask_96 & cond_nobead, 'volume_ml'] = 1.0

    summary = ag.summarize(big)

    # write outputs into Data_analysis
    outdir = Path.cwd() / "Data_analysis"
    outdir.mkdir(parents=True, exist_ok=True)

    # If there is an existing expanded summary in Data_analysis that contains
    # OD columns, merge those OD columns into our summary so the plots include
    # morning/afternoon OD panels.
    try:
        candidates = list(outdir.glob("aggregate_summary_by_volume*.csv"))
        for c in candidates:
            try:
                ext = pd.read_csv(c)
            except Exception:
                continue
            od_cols = {"OD_morning_mean", "OD_morning_stdev", "OD_afternoon_mean", "OD_afternoon_stdev"}
            if od_cols.issubset(set(ext.columns)):
                merge_keys = ["strain", "__source", "bead_present", "volume_ml"]
                ext_small = ext[merge_keys + list(od_cols)].copy()
                for k in merge_keys:
                    if k in ext_small.columns and k in summary.columns:
                        try:
                            ext_small[k] = ext_small[k].astype(summary[k].dtype)
                        except Exception:
                            pass
                summary = summary.merge(ext_small, on=merge_keys, how='left')
                break
    except Exception:
        pass

    csv_path = outdir / "aggregate_summary_by_volume.csv"
    summary.to_csv(csv_path, index=False)
    print(f"Wrote summary CSV to {csv_path}")

    # generate per-strain figures using the plotting functions in the module
    strains = summary["strain"].dropna().unique().tolist()
    preferred = ["SP286", "dea2^"]
    for s in preferred:
        if s in strains:
            safe_name = re.sub(r"[^A-Za-z0-9_-]", "_", s)
            outpath = outdir / f"aggregate_{safe_name}_by_bead_and_metric.png"
            print(f"Plotting strain {s} -> {outpath}")
            ag.plot_strain_facets(summary, s, outpath)

    # metric-specific 'by_volume' plots removed per user request


if __name__ == '__main__':
    main()
