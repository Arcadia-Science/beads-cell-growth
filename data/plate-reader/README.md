# OD Data

Absorbance (ABS) measurements were collected on a SpectraMax iD5 plate reader using a 595 nm wavelength in "kinetic" mode, recording absorbance values every 10 minutes throughout the time course.

Each file in `data/plate-reader/` contains the raw wellscan ABS data, the converted OD values, and the averaged readings across 5 wells per condition. These measurements were manually labelled and split by experiment into `Baseline_ODs.csv`, `Supplement_ODs.csv`, and `Baseline_ODs_stdev.csv`.

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

## Derived CSVs

These files are manually labelled and aggregated from the raw XLSX data above.

| File | Description | Columns |
|------|-------------|---------|
| `Baseline_ODs.csv` | Averaged OD values per condition for baseline experiments (test tubes, 24-well, 96-well) | `Experiment, genotype, volume, time, beads, OD` |
| `Baseline_ODs_stdev.csv` | Mean and standard deviation of OD per condition; each row is a unique experiment/strain/bead/volume combination | `experiment, strain, beads, volume, OD_morning_mean, OD_morning_stdev, OD_afternoon_mean, OD_afternoon_stdev` |
| `Supplement_ODs.csv` | OD values for the supplements experiment. The columns `A`, `V`, and `D` are binary flags indicating the presence (1) or absence (0) of the following supplements: A = Antifoam, V = Vitamin mixture, D = Dextrose. | `Experiment, genotype, time, OD, A, V, D` |
