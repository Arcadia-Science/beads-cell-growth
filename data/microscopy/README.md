# Microscopy Data

Four ND2 files, each a single DIC field of view of *S. pombe* cells acquired on a Nikon Ti2-E with a 40x (0.95 NA) air objective. The channel is labeled as FITC in the filenames due to a metadata error during acquisition; the actual modality is DIC.

These images are a subset from the larger experiment, the data for which is available via [Zenodo](https://zenodo.org/records/18927821). The sample images here comprise a 2×2 design: two strains (WT, DEA2) crossed with two bead conditions (no bead, single 3 mm bead). Each file here is one randomly chosen field of view per condition, selected for use in the example segmentation notebook (`notebooks/segmentation-example.ipynb`).

| File | Strain | Bead condition |
|------|--------|----------------|
| `WellB10_PointB10_0004_ChannelFITC BP_Seq0094.nd2` | WT | No bead |
| `WellB06_PointB06_0002_ChannelFITC BP_Seq0128.nd2` | WT | 3 mm bead |
| `WellF10_PointF10_0007_ChannelFITC BP_Seq0277.nd2` | DEA2 | No bead |
| `WellF06_PointF06_0003_ChannelFITC BP_Seq0309.nd2` | DEA2 | 3 mm bead |
