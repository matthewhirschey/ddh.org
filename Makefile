.PHONY = dirs

# Default command to use for Rscript shall simply be Rscript. When running under a container, this may be
# singularity exec singularity/ddh.sif RScript
RSCRIPT_CMD ?= Rscript

# Default container image registry/repository to use when making the local singularity container
# image file
DOCKER_IMG ?= docker://dukegcb/ddh:latest

# Version of depmap we are using. This value is used in some filenames in this Makefile.
# The value should match the value of the release field in code/current_release.R.
DMVER ?= 20Q1

# Defines the number intermediate pathway data files that will be created. This allow generating the pathway data in parallel.
# Changing this number requires an update to the data/master_positive.Rds and data/master_negative.Rds rules below.
NUM_SUBSET_FILES ?= 10

# Disable parallel build to handle code that creates multiple targets
.NOTPARALLEL:

# The first target is the default, it makes "all" the data. Does not include container_image
all: dirs gene_summary depmap_data depmap_stats depmap_tables depmap_pathways

# The clean target removes all the files in data/
clean:
	rm data/*

clean_container:
	rm -r singularity/*

# Map simple target names to the files on which they depend
container_image: singularity/ddh.sif
gene_summary: data/gene_summary.Rds
depmap_data: data/$(DMVER)_achilles_cor.Rds data/$(DMVER)_achilles.Rds data/$(DMVER)_expression.Rds data/$(DMVER)_expression.Rds data/$(DMVER)_expression_join.Rds
depmap_stats: data/sd_threshold.Rds data/achilles_lower.Rds data/achilles_upper.Rds data/mean_virtual_achilles.Rds data/sd_virtual_achilles.Rds
depmap_tables: data/master_top_table.Rds data/master_bottom_table.Rds
depmap_pathways: data/master_positive.Rds data/master_negative.Rds

dirs:
	mkdir -p data
	mkdir -p singularity/images

singularity/ddh.sif:
	@echo "Pulling container image"
	singularity pull singularity/images/ddh.sif $(DOCKER_IMG)

data/gene_summary.Rds: code/create_gene_summary.R
	@echo "Creating gene summary"
	$(RSCRIPT_CMD) code/create_gene_summary.R --entrezkey ${ENTREZ_KEY}

data/$(DMVER)_achilles_cor.Rds data/$(DMVER)_achilles.Rds data/$(DMVER)_expression.Rds data/$(DMVER)_expression_join.Rds: code/generate_depmap_data.R
	@echo "Creating depmap data"
	$(RSCRIPT_CMD) code/generate_depmap_data.R

data/sd_threshold.Rds data/achilles_lower.Rds data/achilles_upper.Rds data/mean_virtual_achilles.Rds data/sd_virtual_achilles.Rds: data/$(DMVER)_achilles_cor.Rds code/generate_depmap_stats.R
	@echo "Creating depmap stats"
	$(RSCRIPT_CMD) code/generate_depmap_stats.R

data/master_top_table.Rds data/master_bottom_table.Rds: code/generate_depmap_tables.R data/gene_summary.Rds data/$(DMVER)_achilles_cor.Rds data/achilles_lower.Rds data/achilles_upper.Rds
	@echo "Creating depmap tables"
	$(RSCRIPT_CMD) code/generate_depmap_tables.R

data/master_positive.Rds: code/generate_pathways.sh code/generate_depmap_pathways.R code/merge_depmap_pathways.R
	@echo "Creating positive pathways data"
	RSCRIPT_CMD="$(RSCRIPT_CMD)" PATHWAY_TYPE=positive NUM_SUBSET_FILES=$(NUM_SUBSET_FILES) ./code/generate_pathways.sh

data/master_negative.Rds: code/generate_pathways.sh code/generate_depmap_pathways.R code/merge_depmap_pathways.R
	@echo "Creating negative pathways data"
	RSCRIPT_CMD="$(RSCRIPT_CMD)" PATHWAY_TYPE=negative NUM_SUBSET_FILES=$(NUM_SUBSET_FILES) ./code/generate_pathways.sh

