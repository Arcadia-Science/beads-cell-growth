# 2026-bead-growth

## Purpose

TODO: Briefly describe the core analyses performed in the repository and the motivation behind them.

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

TODO: Add details about the description of input / output data and links to Zenodo depositions, if applicable.

## Overview

### Description of the folder structure

### Methods

TODO: Include a brief, step-wise overview of analyses performed.

> Example:
>
> 1.  Download scripts using `download.ipynb`.
> 2.  Preprocess using `./preprocessing.sh -a data/`
> 3.  Run Snakemake pipeline `snakemake --snakefile Snakefile`
> 4.  Generate figures using `pub/make_figures.ipynb`.

### Compute Specifications

TODO: Describe what compute resources were used to run the analysis. For example, you could list the operating system, number of cores, RAM, and storage space.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide--credit-for-contributions.md).
