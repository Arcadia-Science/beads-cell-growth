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
| `data/plate-reader/Baseline_ODs.csv` | Manually labelled and aggregated OD data for baseline experiments (test tubes, 24-well, 96-well) |
| `data/plate-reader/Supplement_ODs.csv` | Manually labelled and aggregated OD data for supplement experiments |
| `data/plate-reader/Baseline_ODs_stdev.csv` | Reformatted version of `Baseline_ODs.csv` that also includes standard deviations |


### Mapping Experiments and Figures to Scripts and Input Data Sources

Each microscopy script reads CSV files containing single-cell morphology measurements (found in processed microscopy data files). These CSVs are produced by processing microscopy images with `scripts/segment_cells.py`, which uses Cellpose for cell segmentation. There is one CSV per field of view, and multiple fields of view per well. The scripts aggregate all CSVs for each well, map well positions to corresponding strains and treatments, compute descriptive statistics, perform pairwise comparisons (using Welch's t-test with Holm correction), and generate the figures included in the publication.

| Figure | Experiment | Script |
|---|---|---|
| Figure 1 | Experimental setup | - |
| Figure 2 | Test tubes | `ttubes_beads.py` |
| Figure 3 | 24-well plates | `24_beads.py` |
| Figure 4 | 96-well plates | `96_beads.py` |
| Figure 5 | Aggregate WT | `aggregate_multi_sources.py` |
| Figure 6 | Supplements | `morning_supplements.py` |

OD heatmaps for Figures 2–4 are generated from `Baseline_ODs.csv` and `Supplement_ODs.csv` by `make_od_heatmaps.py`.

OD line plots for Figure 5 are generated from `Baseline_ODs_stdev.csv` by `aggregate_multi_sources.py`. `Baseline_ODs_stdev.csv` contains the same measurements as `Baseline_ODs.csv`, but reformatted so that each row represents a unique experiment/strain/bead/volume combination, and it additionally includes standard deviation columns for each measurement.

Cell length and area measurements for Figures 2–6 are compiled by aggregating CSV files containing processed microscopy data produced by `segment_cells.py`. All analysis is restricted to the wild type strain (SP286). Each experiment has its own dedicated folder, corresponding to processed microscopy data:

- Test tubes: `20260122_111821_521`
- 24-well plates: `20260122_113404_129`
- 96-well plates: `20260116_094944_372`
- Supplements: `20260123_113447_096`


### Experimental OD Data

For more information, see the [plate reader data README](data/plate-reader/README.md).

### Experimental Microscopy Data

For more information, see the [microscopy data README](data/microscopy/README.md).



### Reproducing Figures

1. Download raw and processed microscopy data from [Zenodo](https://zenodo.org/records/18927821).
2. Process raw images using [arcadia-microscopy-tools](https://github.com/Arcadia-Science/arcadia-microscopy-tools).
3. Run the corresponding script to generate each figure (see table above).

### Compute Specifications

All scripts, except for `segment_cells.py`, were executed on a MacBook Pro with an M3 chip, running macOS Tahoe 26.3.1, and equipped with 36 GB of RAM.

The `segment_cells.py` script was executed via [Modal](https://modal.com) using NVIDIA T4 GPUs on a minimal Debian environment. The setup leveraged up to 10 concurrent containers and mounted an S3 bucket to access the high-throughput microscopy data.

## Contributors

- Román Ramos Báez
- Brae Bigge
- Ben Braverman
- Ryan Lane

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide--credit-for-contributions.md).
