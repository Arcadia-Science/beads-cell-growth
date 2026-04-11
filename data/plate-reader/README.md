# OD Data

Absorbance (ABS) measurements were collected on a SpectraMax iD5 plate reader using a 595 nm wavelength in "kinetic" mode, recording absorbance values every 10 minutes throughout the time course.

Raw wellscan data from the XLSX files below were converted to OD values, manually labelled, and compiled into a single curated CSV (`od-measurements.csv`).

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

## Compiled CSV

| File | Description | Columns |
|------|-------------|---------|
| `od-measurements.csv` | One row per replicate measurement across all experiments (test tubes, 24-well, 96-well, supplements) | `experiment, strain, volume, time, beads, rep, OD, supplement` |
