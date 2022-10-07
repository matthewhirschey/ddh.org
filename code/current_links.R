## THIS FILE SETS PARAMETERS AND VARIABLES FOR DATA GENERATION

achilles_url <- "https://ndownloader.figshare.com/files/34989889" #achilles_gene_effect.csv
ccle_url <- "https://ndownloader.figshare.com/files/34989919" #CCLE_expression.csv
cclemeta_url <- "https://ndownloader.figshare.com/files/35020903" #sample_info.csv
prism_url <- "https://ndownloader.figshare.com/files/17741420" #primary-screen-replicate-collapsed-logfold-change.csv (same since 19Q4)
prismmeta_url <- "https://ndownloader.figshare.com/files/20237715" #primary-screen-replicate-collapsed-treatment-info.csv (same since 19Q4)
achilles_log_url <- "https://ndownloader.figshare.com/files/34989955" #achilles_logfold_change.csv (for drug data)
achilles_guide_map_url <- "https://ndownloader.figshare.com/files/34989880" #achilles_guide_map.csv (for drug data)
achilles_rep_map_url <- "https://ndownloader.figshare.com/files/34989913" #achilles_replicate_map.csv (for drug data)
common_essentials_url <- "https://ndownloader.figshare.com/files/34990024" #common_essentials.csv (for gene dependency data)

fraction_cutoff <- 0.05 #~5% FDR
na_cutoff_file = here::here("data", "na_cutoff.Rds")
if (file.exists(na_cutoff_file)) {
  na_cutoff <- readRDS(file = na_cutoff_file)   
}

#stable URLs
gene_summary_url <- "https://zenodo.org/record/4604643/files/21Q1_gene_summary.Rds"
gene_location_url <- "https://zenodo.org/record/4584874/files/gene_location.Rds"
cids_url <- "https://zenodo.org/record/4719639/files/21Q1_cids.Rds"
compound_descriptions_url <- "https://zenodo.org/record/4661826/files/21Q1_compound_descriptions.Rds" #supporting file for prism_meta
compound_mesh_url <- "https://zenodo.org/record/4661826/files/21Q1_compound_mesh.Rds" #supporting file for prism_meta
metabolites_url <- "https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip"
metabolite_proteins_url <- "https://hmdb.ca/system/downloads/current/hmdb_proteins.zip"
metabolites_cell_url <- "https://data.broadinstitute.org/ccle/CCLE_metabolomics_20190502.csv"
hmdb_meta_url <- "https://zenodo.org/record/5758926/files/hmdb_meta.Rds"
hmdb_metabolites_url <- "https://zenodo.org/record/5758926/files/hmdb_metabolites.Rds"
hmdb_names_url <- "https://zenodo.org/record/5758926/files/hmdb_names.Rds"
hmdb_proteins_url <- "https://zenodo.org/record/5758926/files/hmdb_proteins.Rds"
bioconcepts2pubtator_url <- "ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/bioconcepts2pubtatorcentral.gz"
pubtator_url <- "ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/gene2pubtatorcentral.gz"
pubmed_url <- "https://ftp.ncbi.nlm.nih.gov/pub/pmc/PMC-ids.csv.gz"
gene2pubmed_url <- "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz"
cellosaurus_url <- "https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt"
subcell_url <- "https://www.proteinatlas.org/download/subcellular_location.tsv.zip"
tissue_url <- "https://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip"
go_bp_url <- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2018" #consider updating from here: http://amigo.geneontology.org/amigo/software_list
go_def_url <- "http://purl.obolibrary.org/obo/go.obo"
proteins_url <- "https://zenodo.org/record/4007646/files/proteins.Rds?download=1"
proteins_ccle_info_url <- "https://gygi.hms.harvard.edu/data/ccle/Table_S1_Sample_Information.xlsx" #https://gygi.hms.harvard.edu/publications/ccle.html
proteins_ccle_url <- "https://gygi.hms.harvard.edu/data/ccle/Table_S2_Protein_Quant_Normalized.xlsx" #https://gygi.hms.harvard.edu/publications/ccle.html

# retry settings for enrichr library
enrichr_retries <- 3
enrichr_retry_sleep_seconds <- 30

# retry settings for rentrez library
entrez_retries <- 3
entrez_retry_sleep_seconds <- 30

# keep the following value in sync with the --array flag parameter in
# slurm/data-gene-step4-pos.subsets.sh and slurm/data-gene-step4-neg.subsets.sh
pathways_num_subset_files <- 200
