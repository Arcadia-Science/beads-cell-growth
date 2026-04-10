# 2026-bead-growth

## Purpose

This repository contains the plotting and statistical analysis scripts used to generate figures for the pub ["Adding a single glass bead to liquid cultures improves cell culture quality."](https://doi.org/10.57844/arcadia-mu7g-v6rz) It also includes raw absorbance (ABS) and processed optical density (OD) data from an iD5 plate reader. Microscopy data (processed with [arcadia-microscopy-tools](https://github.com/Arcadia-Science/arcadia-microscopy-tools) and [Cellpose](https://github.com/MouseLand/cellpose)) is available on [Zenodo](https://zenodo.org/records/18927821).

## Installation and Setup

This repository uses [uv](https://docs.astral.sh/uv/getting-started/installation/) to manage software environments and installations. Run the following commands to create the development environment:

```bash
git clone https://github.com/Arcadia-Science/2026-bead-growth.git
cd 2026-bead-growth

uv venv --python 3.12
source .venv/bin/activate

uv sync --all-extras
```


## Overview

### Folder Structure

| Path | Contents |
|---|---|
| `scripts/` | Python scripts for plotting and statistical analysis |
| `data/microscopy/` | Example ND2 microscopy images of *S. pombe* cells and per-experiment cell morphology CSVs (`combined_dic_measurements_*.csv`) |
| `data/plate-reader/` | Per-experiment XLSX files containing raw ABS and processed OD measurements |
| `data/plate-reader/od-measurements.csv` | Aggregated OD data across all baseline and supplement experiments |
| `data/zenodo/` | (optional, git-ignored) Download location for processed microscopy data from Zenodo |


### Mapping Figures to Analysis Scripts

All data needed to reproduce figures is included in the repository. Each analysis notebook reads pre-compiled CSVs from `data/microscopy/` and `data/plate-reader/` directly.

| Figure | Experiment | Analysis script |
|---|---|---|
| Figure 1 | Experimental setup | - |
| Figure 2 | Test tubes | `notebooks/microscopy-analysis-ttube-24well.ipynb` |
| Figure 3 | 24-well plates | `notebooks/microscopy-analysis-ttube-24well.ipynb` |
| Figure 4 | 96-well plates | `notebooks/microscopy-analysis-96-well.ipynb` |
| Figure 5 | Aggregate WT | `scripts/aggregate_multi_sources.py` |
| Figure 6 | Supplements | `notebooks/microscopy-analysis-96-well.ipynb` |

OD analysis is performed in `notebooks/od-analysis.ipynb` using `data/plate-reader/od-measurements.csv`.

The cell morphology CSVs in `data/microscopy/` were compiled by `scripts/compile_microscopy_data.py`, which aggregates per-well measurement files produced by `scripts/segment_cells.py` (Cellpose segmentation). Running `compile_microscopy_data.py` is **not required** to reproduce figures — it is included for full reproducibility of the data-compilation step. See its docstring or `make compile-microscopy` for usage. All analysis is restricted to the wild type strain (SP286). Each experiment has its own dedicated folder on [Zenodo](https://zenodo.org/records/18927821):

- Test tubes: `20260122_111821_521`
- 24-well plates: `20260122_113404_129`
- 96-well plates: `20260116_094944_372`
- Supplements: `20260123_113447_096`


### Experimental OD Data

For more information, see the [plate reader data README](data/plate-reader/README.md).

### Experimental Microscopy Data

For more information, see the [microscopy data README](data/microscopy/README.md).



### Reproducing Figures

All data needed to reproduce figures is included in the repository. After cloning and installing (see above), run the analysis notebooks and scripts listed in the table above.

Optionally, to regenerate the compiled microscopy CSVs from the per-well segmentation outputs:

1. Download processed microscopy data from [Zenodo](https://zenodo.org/records/18927821) into `data/zenodo/`.
2. Run `make compile-microscopy` (this looks for Zenodo data in `data/zenodo/` by default).

### Compute Specifications

All scripts, except for `segment_cells.py`, were executed on a MacBook Pro with an M3 chip, running macOS Tahoe 26.3.1, and equipped with 36 GB of RAM.

The `segment_cells.py` script was executed via [Modal](https://modal.com) using NVIDIA T4 GPUs on a minimal Debian environment. The setup leveraged up to 10 concurrent containers and mounted an S3 bucket to access the high-throughput microscopy data.

## Contributors

- Román Ramos Báez
- Brae Bigge
- Ben Braverman
- Ryan Lane

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide--credit-for-contributions.md).
