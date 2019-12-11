.PHONY = dirs

# Default command to use for Rscript shall simply be Rscript. When running under a container, this may be
# singularity exec singularity/depmap.sif RScript
RSCRIPT_CMD ?= Rscript

# The first target is the default, it makes "all" the data. Does not include container_image
all: gene_summary depmap_data depmap_stats depmap_tables

# The clean target removes all the files in data/
clean:
	rm data/*

clean_container:
	rm -r singularity/*

# Map simple target names to the files on which they depend
container_image: singularity/depmap.sif
gene_summary: data/gene_summary.RData
depmap_data: data/19Q3_achilles_cor.Rdata data/19Q3_achilles.Rdata data/19Q3_expression.Rdata data/19Q3_expression_id.Rdata
depmap_stats: data/sd_threshold.rds data/achilles_lower.rds data/achilles_upper.rds data/mean_virtual_achilles.rds data/sd_virtual_achilles.rds
depmap_tables: data/master_top_table.Rdata data/master_bottom_table.RData

dirs:
	mkdir -p data
	mkdir -p singularity

singularity/depmap.sif:
	@echo "Pulling container image"
	singularity pull singularity/images/depmap.sif ${DOCKER_IMG}

data/gene_summary.RData: code/create_gene_summary.R
	@echo "Creating gene summary"
	$(RSCRIPT_CMD) code/create_gene_summary.R --entrezkey ${ENTREZ_KEY}

data/19Q3_achilles_cor.Rdata data/19Q3_achilles.Rdata data/19Q3_expression.Rdata data/19Q3_expression_id.Rdata: code/generate_depmap_data.R
	@echo "Creating depmap data"
	$(RSCRIPT_CMD) code/generate_depmap_data.R

data/sd_threshold.rds data/achilles_lower.rds data/achilles_upper.rds data/mean_virtual_achilles.rds data/sd_virtual_achilles.rds: data/19Q3_achilles_cor.Rdata code/generate_depmap_stats.R
	@echo "Creating depmap stats"
	$(RSCRIPT_CMD) code/generate_depmap_stats.R

data/master_top_table.Rdata data/master_bottom_table.RData: code/generate_depmap_tables.R data/gene_summary.RData data/19Q3_achilles_cor.Rdata data/achilles_lower.rds data/achilles_upper.rds
	@echo "Creating depmap tables"
	$(RSCRIPT_CMD) code/generate_depmap_tables.R

data/master_positive.RData data/master_negative.RData: code/generate_depmap_pathways.R data/gene_summary.RData data/gene_summary.RData data/19Q3_achilles_cor.Rdata data/achilles_lower.rds data/achilles_upper.rds
	@echo "Creating depmap pathways"
	$(RSCRIPT_CMD) code/generate_depmap_pathways.R
