# Data-driven Hypothesis

Data-driven hypothesis is a resource for identifying novel functions for human genes developed by the [Hirschey Lab](www.hirscheylab.org). A typical use case starts by querying a gene, identifying genes that share similar patterns or behaviors across various measures, in order to discover new novel genes in established processes or new processes for well-studied genes.

## Querying Data

Visit [www.datadrivenhypothesis.org](www.datadrivenhypothesis.org) to query genes, generate data, and download reports.

## Methods

Visit [www.datadrivenhypothesis.org](www.datadrivenhypothesis.org) to view detailed methods for this project. Alternatively, after generating data files, knit `code/methods.Rmd` for the complete methods summary.

## Generating/Recreating Data

The .R scripts in the `/code` directory are used to generate files in the `/data` directory.

### Using Make

This project includes a GNU [Makefile](https://www.gnu.org/software/make/manual/html_node/Introduction.html) that can execute the `.R` scripts and generate data files as needed. If you have [R](https://www.r-project.org) and the required libraries installed, you can generate the data by running `make` on the command-line. Note that the gene\_summary expects your NCBI API key to be set in the environment variable `ENTREZ_KEY`:

```
ENTREZ_KEY="your-key-here" make
```

By default, this will make the `all` target in the [Makefile](Makefile), which produces all of the data files. To generate an individual data file or set of them (e.g. `depmap_tables`), you can simply `make` that target:

```
make depmap_data
make data/sd_threshold.rds
```

To remove all data files, use `make clean`

### Singularity

This project also supports using a container runtime like [singularity](https://sylabs.io/singularity/) to run `Rscript`. To run under Singularity, set the `RSCRIPT_CMD` environment variable as noted in [build-slurm.sh](build-slurm.sh). This scripts expects site-specific environment variables to be exported from a `config.sh` file. This file is not included in the repo, as it

The [Makefile](Makefile) also includes a `container_image` target that will download the docker image named in the `DOCKER_IMG` variable and prepare it for singularity.

### Data file generation

Alternatively, you can run the scripts manually:

To generate the data files, run:
1. code/create_gene_summary.R
2. code/generate_depmap_data.R
3. code/generate_depmap_stats.R
4. code/generate_depmap_tables.R
5. code/generate_depmap_pathways.R --type positive --num-subset-files 1 --idx 1
6. code/merge_depmap_pathways.R --type positive --num-subset-files 1
7. code/generate_depmap_pathways.R --type negative --num-subset-files 1 --idx 1
8. code/merge_depmap_pathways.R --type positive --num-subset-files 1

The files generated in steps 1-3 are required for steps 4 and 5. Step 4 takes about 60' to run locally. Steps 5 and 7 take multiple days if run as specified above. By increasing the `--num-subset-files` value you can run multiple of these processes in parallel. For example to run 4 processes of step 5:
```
code/generate_depmap_pathways.R --type positive --num-subset-files 4 --idx 1 &
code/generate_depmap_pathways.R --type positive --num-subset-files 4 --idx 2 &
code/generate_depmap_pathways.R --type positive --num-subset-files 4 --idx 3 & 
code/generate_depmap_pathways.R --type positive --num-subset-files 4 --idx 4 &
wait
```

