# 2026-bead-growth

## Purpose

This repository contains the plotting and statistical analysis scripts used to generate figures for the pub ["Adding a single glass bead to liquid cultures improves cell culture quality."](https://doi.org/10.57844/arcadia-mu7g-v6rz) It also includes raw absorbance (ABS) and processed optical density (OD) data from an iD5 plate reader. Microscopy data (processed with [arcadia-microscopy-tools](https://github.com/Arcadia-Science/arcadia-microscopy-tools) and Cellpose) is available on [Zenodo](https://zenodo.org/records/18927821).

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
| `data/experimental/` | Per-experiment XLSX files containing raw ABS and processed OD measurements |
| `data/Baseline_ODs.xlsx` | Manually labelled and aggregated OD data for baseline experiments (test tubes, 24-well, 96-well) |
| `data/Supplement_ODs.xlsx` | Manually labelled and aggregated OD data for supplement experiments |


### Mapping Experiments and Figures to Scripts

Each microscopy script reads CSV files containing single-cell morphology measurements. These CSVs are produced by processing microscopy images with Cellpose; there is one CSV per field of view, and multiple fields of view per well. The scripts aggregate all CSVs for each well, map well positions to corresponding strains and treatments, compute descriptive statistics, perform pairwise comparisons (using Welch's t-test with Holm correction), and generate the figures included in the publication.

| Figure | Experiment | Script |
|---|---|---|
| Figure 1 | Experimental setup | - |
| Figure 2 | Test tubes | `ttubes_beads.py` |
| Figure 3 | 24-well plates | `24_beads.py` |
| Figure 4 | 96-well plates | `96_beads.py` |
| Figure 5 | Aggregate WT | `aggregate_multi_sources.py` |
| Figure 6 | Aggregate mutant | `aggregate_multi_sources.py` |
| Figure 7 | Supplements | `morning_supplements.py` |

OD heatmaps for Figures 2-6 are generated from `Baseline_ODs.xlsx` and `Supplement_ODs.xlsx` by `make_od_heatmaps.py`.


### Experimental OD Data

Each file in `data/experimental/` contains the raw wellscan ABS data, the converted OD values, and the averaged readings across 5 wells per condition. These measurements were manually labelled and split by experiment into `Baseline_ODs.xlsx` and `Supplement_ODs.xlsx`.

| Experiment | File | Timepoint |
|---|---|---|
| 96-well plates | `260115_rom_96well_beads_day1.xlsx` | Day 1 |
| 96-well plates | `260116_rom_96well_beads_day2.xlsx` | Day 2 |
| Test tubes | `260122_rom_ttubes_pombe_beads_pub.xlsx` | Morning |
| Test tubes | `260122_rom_ttubes_2_pombe_beads_pub.xlsx` | Afternoon |
| 24-well plates | `260122_rom_24plates_pombe_beads_pub.xlsx` | Morning |
| 24-well plates | `260122_rom_24plates_2_pombe_beads_pub.xlsx` | Afternoon |
| Supplements | `260123_rom_supplements_pombe_beads_pub.xlsx` | Morning |
| Supplements | `260123_rom_supplements_2_pombe_beads_pub.xlsx` | Afternoon |

### TODO:

1. Is the processed microscopy data (or at least the output CSV files) being uploaded to Zenodo? (then people wouldn't have to process the microscopy data themselves).
2. Update section below accordingly.
3. Add modal script and example notebook.
4. Map JOBS run to script
5. Add compute specs (Macbook + Modal GPU info)


### Reproducing Figures

1. Download microscopy data from [Zenodo](https://zenodo.org/records/18927821).
2. Process raw images using [arcadia-microscopy-tools](https://github.com/Arcadia-Science/arcadia-microscopy-tools).
3. Run the corresponding script to generate each figure (see table above).

### Compute Specifications



## Contributors

- Román Ramos Báez
- Brae Bigge
- Ben Braverman
- Ryan Lane

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide--credit-for-contributions.md).
