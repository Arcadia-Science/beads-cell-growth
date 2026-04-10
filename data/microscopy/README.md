# Microscopy Data

Two ND2 files, each a single DIC field of view of *S. pombe* cells acquired on a Nikon Ti2-E with a 40x (0.95 NA) air objective. The channel is labeled as FITC in the filenames due to a metadata error during acquisition; the actual modality is DIC.

These images are a subset from the larger 96-well experiment (`20260116_094944_372`), the data for which is available via [Zenodo](https://zenodo.org/records/18927821). The sample images here comprise two bead conditions (no bead, single 3 mm bead) for the wild type strain (WT/SP286). Each file here is one randomly chosen field of view per condition, selected for use in the example segmentation notebook (`notebooks/segmentation-example.ipynb`).

| File | Strain | Bead condition |
|------|--------|----------------|
| `WellB10_PointB10_0004_ChannelFITC BP_Seq0094.nd2` | WT | No bead |
| `WellB06_PointB06_0002_ChannelFITC BP_Seq0128.nd2` | WT | 3 mm bead |

## Cell Morphology CSVs

Pre-compiled per-cell morphology measurements (axis major length, area, cell length) extracted from DIC microscopy images via the segmentation pipeline (`scripts/segment_cells.py`). These CSVs are committed to the repository so that analysis notebooks can be run without any additional data processing.

Each CSV can optionally be regenerated from per-well segmentation outputs using `scripts/compile_microscopy.py`:

| File | Producing command |
|------|------------------|
| `combined_dic_measurements_96well.csv` | `scripts/compile_microscopy.py -e 96-well` |
| `combined_dic_measurements_24well.csv` | `scripts/compile_microscopy.py -e 24-well` |
| `combined_dic_measurements_ttubes.csv` | `scripts/compile_microscopy.py -e ttubes` |
| `combined_dic_measurements_supplements.csv` | `scripts/compile_microscopy.py -e supplements` |
