.PHONY = all dirs

SINGULARITY_IMG = singularity/depmap.sif
SINGULARITY_EXEC = singularity exec ${SINGULARITY_IMG}

container: singularity/depmap.sif
gene_summary: data/gene_summary.feather
depmap_data: data/19Q3_achilles_cor.Rdata
depmap_stats: data/sd_threshold.rds

dirs:
	mkdir -p data
	mkdir -p singularity

singularity/depmap.sif:
	@echo "Pulling container image"
	singularity pull singularity/depmap.sif docker://dleehr/depmap:latest

data/gene_summary.feather: code/create_gene_summary.R
	@echo "Creating gene summary"
	${SINGULARITY_EXEC} Rscript code/create_gene_summary.R

data/19Q3_achilles_cor.Rdata: code/generate_depmap_data.R
	@echo "Creating depmap data"
	${SINGULARITY_EXEC} Rscript code/generate_depmap_data.R

data/sd_threshold.rds: code/generate_depmap_stats.R
	@echo "Creating depmap stats"
	${SINGULARITY_EXEC} Rscript code/generate_depmap_stats.R


all: dirs container depmap_data depmap_stats gene_summary
