## THIS FILE SETS PARAMETERS AND VARIABLES FOR DATA GENERATION

#update each release
release <- "21Q1"
achilles_url <- "https://ndownloader.figshare.com/files/26261293" #achilles_gene_effect.csv
ccle_url <- "https://ndownloader.figshare.com/files/26261476" #CCLE_expression.csv
cclemeta_url <- "https://ndownloader.figshare.com/files/26261569" #sample_info.csv
prism_url <- "https://ndownloader.figshare.com/files/17741420" #primary-screen-replicate-collapsed-logfold-change.csv (same as 19Q4)
prismmeta_url <- "https://ndownloader.figshare.com/files/20237715" #primary-screen-replicate-collapsed-treatment-info.csv (same as 19Q4)
achilles_log_url <- "https://ndownloader.figshare.com/files/26261431" #achilles_logfold_change.csv (for drug data)
achilles_guide_map_url <- "https://ndownloader.figshare.com/files/26261287" #achilles_guide_map.csv (for drug data)
achilles_rep_map_url <- "https://ndownloader.figshare.com/files/26261434" #achilles_replicate_mpa.csv (for drug data)

fraction_cutoff <- 0.05 #~5% FDR
na_cutoff_file = here::here("data", paste0(release, "_na_cutoff.Rds"))
if (file.exists(na_cutoff_file)) {
  na_cutoff <- readRDS(file = na_cutoff_file)   
}

#stable URLs
metabolites_url <- "https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip"
pubtator_url <- "ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/gene2pubtatorcentral.gz"
gene2pubmed_url <- "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz"
subcell_url <- "https://www.proteinatlas.org/download/subcellular_location.tsv.zip"
go_bp_url <- "https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2018" #consider updating from here: http://amigo.geneontology.org/amigo/software_list
go_def_url <- "http://purl.obolibrary.org/obo/go.obo"
proteins_url <- "https://zenodo.org/record/4007646/files/proteins.Rds?download=1"

# retry settings for enrichr library
enrichr_retries <- 3
enrichr_retry_sleep_seconds <- 30

# retry settings for rentrez library
entrez_retries <- 3
entrez_retry_sleep_seconds <- 30

# keep the following value in sync with the --array flag parameter in
# slurm/data-gene-step4-pos.subsets.sh and slurm/data-gene-step4-neg.subsets.sh
pathways_num_subset_files <- 200
