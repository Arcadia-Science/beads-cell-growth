# 2026-bead-growth

## Purpose

The following repository has both the raw microscopy data .nds files and the graphing .py files that allow us to take processed data (using arcadia-microscopy-tools) and graph it as shown in the publication. .xlsx files include both the raw data ABS data and converted OD averages.

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

Input data:

TODO: Add details about the description of input / output data and links to Zenodo depositions, if applicable.

## Overview

### Description of the folder structure

### Methods

1. Download data from Zenodo: https://zenodo.org/records/18927821?preview=1&token=eyJhbGciOiJIUzUxMiIsImlhdCI6MT[…]R4olNzyETjgfAdhh6FqKadaQRiFqIAffK-UMYqMykTpXRMu9xEf7siIjj83dHIA.
2. Process data using arcadia-microscopy-files.
3. Generate figures using .py files listed here to map, statistically analyze, and plot processed data as shown in a given figure.

### Compute Specifications

No compute resources were used to process this data.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide--credit-for-contributions.md).
