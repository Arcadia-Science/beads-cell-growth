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


### Mapping Experiments and Figures to Scripts and Input Data Sources

Each microscopy script reads CSV files containing single-cell morphology measurements (found in processed microscopy data files). These CSVs are produced by processing microscopy images with `scripts/segment_cells.py`, which uses Cellpose for cell segmentation. There is one CSV per field of view, and multiple fields of view per well. The scripts aggregate all CSVs for each well, map well positions to corresponding strains and treatments, compute descriptive statistics, perform pairwise comparisons (using Welch's t-test with Holm correction), and generate the figures included in the publication.

| Figure | Experiment | Data compilation | Analysis notebook |
|---|---|---|---|
| Figure 1 | Experimental setup | - | - |
| Figure 2 | Test tubes | `compile_microscopy.py -e ttubes` | `microscopy-analysis-ttube-24well.ipynb` |
| Figure 3 | 24-well plates | `compile_microscopy.py -e 24-well` | `microscopy-analysis-ttube-24well.ipynb` |
| Figure 4 | 96-well plates | `compile_microscopy.py -e 96-well` | `microscopy-analysis-96-well.ipynb` |
| Figure 5 | Aggregate WT | - | `aggregate_multi_sources.py` |
| Figure 6 | Supplements | `compile_microscopy.py -e supplements` | `microscopy-analysis-96-well.ipynb` |

OD analysis is performed in `notebooks/od-analysis.ipynb` using `data/plate-reader/od-measurements.csv`.

Cell length and area measurements for Figures 2–6 are compiled by `scripts/compile_microscopy.py`, which aggregates per-well CSV files containing processed microscopy data produced by `segment_cells.py`. All analysis is restricted to the wild type strain (SP286). Each experiment has its own dedicated folder, corresponding to processed microscopy data:

- Test tubes: `20260122_111821_521`
- 24-well plates: `20260122_113404_129`
- 96-well plates: `20260116_094944_372`
- Supplements: `20260123_113447_096`


### Experimental OD Data

For more information, see the [plate reader data README](data/plate-reader/README.md).

### Experimental Microscopy Data

For more information, see the [microscopy data README](data/microscopy/README.md).



### Reproducing Figures

The combined microscopy CSVs (`data/microscopy/combined_dic_measurements_*.csv`) and OD data (`data/plate-reader/od-measurements.csv`) are included in the repository, so the analysis notebooks can be run immediately after cloning and installing.

To regenerate the combined microscopy CSVs from scratch:

1. Download processed microscopy data from [Zenodo](https://zenodo.org/records/18927821).
2. Run `compile_microscopy.py` once per experiment (see `make compile-microscopy` or the docstring in `scripts/compile_microscopy.py` for exact commands).
3. Run the corresponding analysis notebook or script to generate each figure (see table above).

### Compute Specifications

All scripts, except for `segment_cells.py`, were executed on a MacBook Pro with an M3 chip, running macOS Tahoe 26.3.1, and equipped with 36 GB of RAM.

The `segment_cells.py` script was executed via [Modal](https://modal.com) using NVIDIA T4 GPUs on a minimal Debian environment. The setup leveraged up to 10 concurrent containers and mounted an S3 bucket to access the high-throughput microscopy data.

## Contributors

- Román Ramos Báez
- Brae Bigge
- Ben Braverman
- Ryan Lane

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide--credit-for-contributions.md).
