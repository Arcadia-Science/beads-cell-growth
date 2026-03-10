# 2026-bead-growth

## Purpose
The following repository has the graphing .py files that allow us to take processed data (using arcadia-microscopy-tools) and graph it as shown in the publication. .xlsx files include both the raw data ABS data and converted OD averages.

## Installation and Setup
This repository uses uv to manage software environments and installations. You can find operating system-specific instructions for installing uv [here](https://docs.astral.sh/uv/getting-started/installation/). After installing uv, run the following commands to create the pipeline run environment.

```{bash}
# Clone the repository
git clone https://github.com/Arcadia-Science/2026-bead-growth.git
cd 2026-bead-growth

# Create and activate a virtual environment
uv venv --python 3.12
source .venv/bin/activate

# Install dependencies
uv sync --all-extras
```

## Data
ODs
Input data: All raw OD data was processed, and then manually labelled and split by experiment into 'All ODs.xlsx' and 'Supplement ODs.xlsx'. 
Processing: These .xlsx files were then graphed in Figures 2-6 using 'make_od_heatmaps.py.'
Output: Figure 2-5 heatmaps and line graphs in Figure 6.

Microscopy
Input data: All raw data, avaiable in Zenodo https://zenodo.org/records/18927821?preview=1&token=eyJhbGciOiJIUzUxMiIsImlhdCI6MT[…]R4olNzyETjgfAdhh6FqKadaQRiFqIAffK-UMYqMykTpXRMu9xEf7siIjj83dHIA, was processed using cellpose. Then, the generated .csv files were used as inputs in Figures 2-6 using the following scripts:
Figure 2. 24-WELL PLATES => 24_beads.py
Figure 3. 96-WELL PLATES => 96_beads.py
Figure 5. => aggregate_multi_sources.py
Figure 6. SUPPLEMENTS => morning_supplements.py
Figure 1. TEST TUBES => ttubes_beads.py
Output: Figure 2-5 bar graphs and line graphs in Figure 6.

## Overview
Microscopy
All microscopy scripts are intuitively titled to match data from the 96-WELL PLATES, 24-WELL PLATES, SUPPLEMENTS, and TEST TUBES experiments. These scripts are for plotting and staistical analysis.

ODs
96-WELL PLATES
~260116_rom_96well_beads_day1.xlsx => This is the wellscan ABS data collected on Day1 for the 96-well plate OD measurements. Inside it is the data converted to ODs and the data averaged across the 5 well readings.
~260116_rom_96well_beads_day2.xlsx => This is the wellscan ABS data collected on Day2 for the 96-well plate OD measurements. Inside it is the data converted to ODs and the data averaged across the 5 well readings. test tubes.

TEST TUBES
~260122_rom_ttubes_pombe_beads_pub.xlsx => This is the wellscan ABS data collected in the morning for the test tube OD measurements. Inside it is the data converted to ODs and the data averaged across the 5 well readings.
~260122_rom_ttubes_2_pombe_beads_pub.xlsx => This is the wellscan ABS data collected in the afternoon for the test tube OD measurements. Inside it is the data converted to ODs and the data averaged across the 5 well readings.

24-WELL PLATES
~260122_rom_24plates_pombe_beads_pub.xlsx => This is the wellscan ABS data collected in the morning for the
24-well plate OD measurements. Inside it is the data converted to ODs and the data averaged across the 5 well readings.
~260122_rom_24plates_2_pombe_beads_pub.xlsx => This is the wellscan ABS data collected in the afternoon for the 24-well plate OD measurements. Inside it is the data converted to ODs and the data averaged across the 5 well readings.

SUPPLEMENTS
~260123_rom_supplements_pombe_beads_pub.xlsx => This is the wellscan ABS data collected in the morning for the 24-well plate supplements OD measurements. Inside it is the data converted to ODs and the data averaged across the 5 well readings.
~260123_rom_supplements_2_pombe_beads_pub.xlsx => This is the wellscan ABS data collected in the afternoon for the 24-well plate supplements OD measurements. Inside it is the data converted to ODs and the data averaged across the 5 well readings.

### Description of the folder structure
'Scripts' This folder contains all scripts used for plotting data.
'ODs' This folder contains all raw (ABS) and processed (OD) data from the id5 plate reader.

### Methods
1. Download data from ZENODO LINKS HERE.
2. Process data using arcadia-microscopy-files.
3. Generate figures using .py files listed here to map, statistically analyze, and plot processed data as shown in a given figure.

### Compute Specifications
No compute resources were used to process this data.

## Contributing
Román Ramos Báez
Brae Bigge
Ben Braverman
Ryan Lane
