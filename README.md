# Data-driven Hypothesis

Data-driven hypothesis is a resource for identifying novel functions for human genes developed by the [Hirschey Lab](http://www.hirscheylab.org). A typical use case starts by querying a gene, identifying genes that share similar patterns or behaviors across various measures, in order to discover new novel genes in established processes or new processes for well-studied genes.

## Querying Data

Visit [www.datadrivenhypothesis.org](http://www.datadrivenhypothesis.org) to query genes, generate data, and download reports.

## Methods

Visit [www.datadrivenhypothesis.org](http://www.datadrivenhypothesis.org) to view detailed methods for this project. Alternatively, after generating data files, knit `code/methods.Rmd` for the complete methods summary.

## Generating/Recreating Data

The .R scripts in the `/code` directory are used to generate files in the `/data` directory.
If you have [R](https://www.r-project.org), the required libraries installed you can generate the data by running:
```
mkdir data
Rscript code/generate_data.R
```

If you have an NCBI API key this can speed up portions of the data generation process. Before running `Rscript code/generate_data.R` export your NCBI API key to an environment variable named `ENTREZ_KEY` like so:
```
export ENTREZ_KEY="your-key-here"
```

#### Singularity and Slurm

On systems that support the [singularity](https://sylabs.io/singularity/) container runtime, we recommend using our container image to run the data-generating R scripts with the required dependencies.

Generating data under singularity is also demonstrated in the [submit-slurm-jobs.sh](submit-slurm-jobs.sh). It requires a `config.sh` script. This config script has three optional environment variables: `ENTREZ_KEY`, `DDH_EMAIL`, `DDH_UPLOAD_RESULTS`.
An example config.sh script that will email you about job progress and upload the results looks like so:
```
export ENTREZ_KEY="your-key-here"
export DDH_EMAIL="your-email@email.com"
export DDH_UPLOAD_RESULTS="TRUE"
```
This `DDH_UPLOAD_RESULTS` requires you to create a `.ddsclient` configuration file in your home directory on the cluster. See [instructions to setup a DukeDS config file](https://github.com/Duke-GCB/DukeDSClient/wiki/Agent-User-Keys-(setup)).
to set the `ENTREZ_KEY` environment variable. You can also speicfy a `DDH_EMAIL` email address to be emailed with job progress. 

This script can be run like so from a slurm cluster:
```
./submit-slurm-jobs.sh
```
This script will submit multiple jobs that will run one after another with differing requirements.
When this script finishes it will upload the results to DukeDS service.


