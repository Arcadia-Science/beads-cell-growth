"""Compile per-well microscopy CSVs into a single combined table.

This script is **optional**. The compiled CSVs it produces are already committed
to ``data/microscopy/``, so the analysis notebooks can be run immediately after
cloning. It is included for full reproducibility of the data-compilation step.

Usage
-----
Download processed microscopy data from Zenodo (https://zenodo.org/records/18927821)
into ``data/zenodo/``, then run ``make compile-microscopy`` or invoke per experiment::

    python scripts/compile_microscopy.py -e 96-well \
        -d data/zenodo/20260116_094944_372/processed
    python scripts/compile_microscopy.py -e ttubes \
        -d data/zenodo/20260122_111821_521/processed
    python scripts/compile_microscopy.py -e 24-well \
        -d data/zenodo/20260122_113404_129/processed
    python scripts/compile_microscopy.py -e supplements \
        -d data/zenodo/20260123_113447_096/processed

Each invocation writes a combined CSV to ``data/microscopy/``.
"""

from __future__ import annotations
import argparse
import re
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
WELL_RE = re.compile(r"^Well([A-Z])(\d{2})")

# ---------------------------------------------------------------------------
# Plate-layout mapping functions
# ---------------------------------------------------------------------------


def _map_96well(letter: str, num: int) -> tuple[str, str]:
    """Rows A-D -> SP286; column pairs encode bead diameter."""
    letter = letter.upper().strip()
    if letter not in {"A", "B", "C", "D"}:
        return "unknown", "unknown"
    if not (1 <= num <= 10):
        return "unknown", "unknown"

    if num <= 2:
        treatment = "1 mm bead"
    elif num <= 4:
        treatment = "1.5 mm bead"
    elif num <= 6:
        treatment = "3 mm bead"
    elif num <= 8:
        treatment = "4.5 mm bead"
    else:
        treatment = "NO bead"

    return "SP286", treatment


def _map_volume_bead(letter: str, num: int) -> tuple[str, str]:
    """Rows D-F -> SP286; cols 1-5 no bead, 6-10 bead, volume by position."""
    letter = letter.upper().strip()
    if letter not in {"D", "E", "F"}:
        return "unknown", "unknown"
    if not (1 <= num <= 10):
        return "unknown", "unknown"

    vol = ((num - 1) % 5) + 1
    bead = "bead" if num >= 6 else "no bead"
    return "SP286", f"{vol} mL {bead}"


def _map_supplements(letter: str, num: int) -> tuple[str, str]:
    """Cols 7-12 -> SP286; row x column sub-block determines supplement combo."""
    letter = letter.upper().strip()
    if not (7 <= num <= 12):
        return "unknown", "unknown"

    pos = ((num - 1) % 6) + 1
    left = {"A": "NONE", "B": "D", "C": "V", "D": "A"}
    right = {"A": "D+A", "B": "V+A", "C": "D+V", "D": "D+V+A"}
    treatment = (left if pos <= 3 else right).get(letter, "unknown")
    return "SP286", treatment


# ---------------------------------------------------------------------------
# Post-processing helpers (add derived columns after concat)
# ---------------------------------------------------------------------------


def _add_volume_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Derive volume_ml, volume_group, and bead_present from treatment strings."""
    vol_match = df["treatment"].str.extract(r"^(\d+)")
    df["volume_ml"] = pd.to_numeric(vol_match[0], errors="coerce")
    df["volume_group"] = df["volume_ml"].apply(
        lambda v: chr(96 + int(v)) if not pd.isna(v) and 1 <= int(v) <= 26 else pd.NA
    )
    no_bead = df["treatment"].str.contains(r"\bno\s*bead\b", case=False, regex=True)
    has_bead = df["treatment"].str.contains(r"\bbead\b", case=False, regex=True)
    df["bead_present"] = has_bead & ~no_bead
    return df


# ---------------------------------------------------------------------------
# Experiment configurations
# ---------------------------------------------------------------------------

_VOLUME_BEAD_TREATMENTS: list[str] = []
for _v in (1, 2, 3, 4, 5):
    _VOLUME_BEAD_TREATMENTS.append(f"{_v} mL no bead")
    _VOLUME_BEAD_TREATMENTS.append(f"{_v} mL bead")


@dataclass
class ExperimentConfig:
    name: str
    map_well: Callable[[str, int], tuple[str, str]]
    treatment_order: list[str]
    output_csv: Path
    strain_order: list[str] = field(default_factory=lambda: ["SP286"])
    post_process: Callable[[pd.DataFrame], pd.DataFrame] | None = None


EXPERIMENTS: dict[str, ExperimentConfig] = {
    "96-well": ExperimentConfig(
        name="96-well",
        map_well=_map_96well,
        treatment_order=[
            "NO bead",
            "1 mm bead",
            "1.5 mm bead",
            "3 mm bead",
            "4.5 mm bead",
        ],
        output_csv=REPO_ROOT / "data" / "microscopy" / "combined_dic_measurements_96well.csv",
    ),
    "24-well": ExperimentConfig(
        name="24-well",
        map_well=_map_volume_bead,
        treatment_order=list(_VOLUME_BEAD_TREATMENTS),
        output_csv=REPO_ROOT / "data" / "microscopy" / "combined_dic_measurements_24well.csv",
        post_process=_add_volume_columns,
    ),
    "ttubes": ExperimentConfig(
        name="ttubes",
        map_well=_map_volume_bead,
        treatment_order=list(_VOLUME_BEAD_TREATMENTS),
        output_csv=REPO_ROOT / "data" / "microscopy" / "combined_dic_measurements_ttubes.csv",
        post_process=_add_volume_columns,
    ),
    "supplements": ExperimentConfig(
        name="supplements",
        map_well=_map_supplements,
        treatment_order=["NONE", "D", "V", "A", "D+A", "V+A", "D+V", "D+V+A"],
        output_csv=REPO_ROOT / "data" / "microscopy" / "combined_dic_measurements_supplements.csv",
    ),
}

# ---------------------------------------------------------------------------
# Shared pipeline
# ---------------------------------------------------------------------------


def parse_well_from_name(name: str) -> tuple[str | None, int | None]:
    m = WELL_RE.match(name)
    if not m:
        return None, None
    return m.group(1), int(m.group(2))


def read_one_csv(p: Path) -> pd.DataFrame:
    df = pd.read_csv(p)

    axis_candidates = ["axis_major_length", "axis_major_length_um", "axis_major_length_px"]
    area_candidates = ["area", "area_um2", "area_px2"]

    axis_col = next((c for c in axis_candidates if c in df.columns), None)
    if axis_col is None:
        raise KeyError(
            f"{p.name} missing required column: axis_major_length (or variants). "
            f"Columns: {', '.join(df.columns)}"
        )

    out = pd.DataFrame({"axis_major_length": pd.to_numeric(df[axis_col], errors="coerce")})

    area_col = next((c for c in area_candidates if c in df.columns), None)
    if area_col is not None:
        out["area"] = pd.to_numeric(df[area_col], errors="coerce")
    else:
        out["area"] = pd.NA

    out["length"] = out["axis_major_length"]
    return out


def find_matching_csvs(
    processed_dir: Path,
    map_well: Callable[[str, int], tuple[str, str]],
) -> list[Path]:
    files: list[Path] = []
    for p in processed_dir.glob("*.csv"):
        if not p.name.startswith("Well"):
            continue
        letter, num = parse_well_from_name(p.name)
        if letter is None or num is None:
            continue
        strain, _ = map_well(letter, num)
        if strain != "unknown":
            files.append(p)
    return sorted(files)


def build_combined_dataframe(
    csv_paths: list[Path],
    config: ExperimentConfig,
) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    for p in csv_paths:
        letter, num = parse_well_from_name(p.name)
        if letter is None or num is None:
            continue

        strain, treatment = config.map_well(letter, num)
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
    out["length"] = out["axis_major_length"]

    if config.post_process is not None:
        out = config.post_process(out)

    out["strain"] = pd.Categorical(out["strain"], categories=config.strain_order, ordered=True)
    out["treatment"] = pd.Categorical(
        out["treatment"], categories=config.treatment_order, ordered=True
    )
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Compile per-well microscopy CSVs into a combined table."
    )
    parser.add_argument(
        "--experiment",
        "-e",
        required=True,
        choices=sorted(EXPERIMENTS),
        help="Experiment type (determines plate-layout mapping and output path)",
    )
    parser.add_argument(
        "--processed-dir",
        "-d",
        type=Path,
        default=Path.cwd() / "processed",
        help="Directory containing per-well CSV files (default: ./processed)",
    )
    args = parser.parse_args(argv)

    config = EXPERIMENTS[args.experiment]
    processed: Path = args.processed_dir

    if not processed.exists() or not processed.is_dir():
        print(f"Error: processed directory not found: {processed}")
        return 1

    csvs = find_matching_csvs(processed, config.map_well)
    if not csvs:
        print(f"Error: no matching Well*.csv files found in {processed}")
        return 1

    df = build_combined_dataframe(csvs, config)
    if df.empty:
        print(
            "Error: no data after mapping wells to strain/treatment."
            " Check filenames + mapping rules."
        )
        return 1

    config.output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(config.output_csv, index=False)
    print(f"Saved combined data ({len(df)} rows): {config.output_csv}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
