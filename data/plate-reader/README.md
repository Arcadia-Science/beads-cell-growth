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
