.PHONY = dirs

# Default command to use for Rscript shall simply be Rscript. When running under a container, this may be
# singularity exec singularity/depmap.sif RScript
RSCRIPT_CMD ?= Rscript

# The first target is the default, it makes "all" the data. Does not include container_image
all: gene_summary depmap_data depmap_stats depmap_tables depmap_pathways

# The clean target removes all the files in data/
clean:
	rm data/*

clean_container:
	rm -r singularity/*

# Map simple target names to the files on which they depend
container_image: singularity/depmap.sif
gene_summary: data/gene_summary.RData
depmap_data: data/19Q4_achilles_cor.RData data/19Q4_achilles.RData data/19Q4_expression.RData data/19Q4_expression.RData data/19Q4_expression_join.RData
depmap_stats: data/sd_threshold.Rds data/achilles_lower.Rds data/achilles_upper.Rds data/mean_virtual_achilles.Rds data/sd_virtual_achilles.Rds
depmap_tables: data/master_top_table.RData data/master_bottom_table.RData
depmap_pathways: data/master_positive.RData data/master_negative.RData

dirs:
	mkdir -p data
	mkdir -p singularity

singularity/depmap.sif:
	@echo "Pulling container image"
	singularity pull singularity/images/depmap.sif ${DOCKER_IMG}

data/gene_summary.RData: code/create_gene_summary.R
	@echo "Creating gene summary"
	$(RSCRIPT_CMD) code/create_gene_summary.R --entrezkey ${ENTREZ_KEY}

data/19Q4_achilles_cor.RData data/19Q4_achilles.RData data/19Q4_expression.RData data/19Q4_expression_join.RData: code/generate_depmap_data.R
	@echo "Creating depmap data"
	$(RSCRIPT_CMD) code/generate_depmap_data.R

data/sd_threshold.Rds data/achilles_lower.Rds data/achilles_upper.Rds data/mean_virtual_achilles.Rds data/sd_virtual_achilles.Rds: data/19Q4_achilles_cor.RData code/generate_depmap_stats.R
	@echo "Creating depmap stats"
	$(RSCRIPT_CMD) code/generate_depmap_stats.R

data/master_top_table.RData data/master_bottom_table.RData: code/generate_depmap_tables.R data/gene_summary.RData data/19Q4_achilles_cor.RData data/achilles_lower.Rds data/achilles_upper.Rds
	@echo "Creating depmap tables"
	$(RSCRIPT_CMD) code/generate_depmap_tables.R

data/master_positive.RData data/master_negative.RData: code/generate_depmap_pathways.R data/gene_summary.RData data/gene_summary.RData data/19Q4_achilles_cor.RData data/achilles_lower.Rds data/achilles_upper.Rds
	@echo "Creating depmap pathways"
	$(RSCRIPT_CMD) code/generate_depmap_pathways.R
