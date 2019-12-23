## THIS CODE GENERATES positive_subset_*_of_*.Rds or negative_subset__*_of_*.Rds files.
## This script is meant to be run multiple times in parallel to quickly calculate the data.
## Commands to create the first file for both postiive and negative if we create 10 files for each:
##  Rscript code/generate_depmap_pathways.R --type positive --idx 1 --num-subset-files 10
##  Rscript code/generate_depmap_pathways.R --type negative --idx 1 --num-subset-files 10
## The files created by this script can be merged by running merge_depmap_pathways.R

#load libraries
library(tidyverse)
library(here)
library(corrr)
library(enrichR)
library(optparse)

#read current release information to set parameters for processing
source(here::here("code", "current_release.R"))
#read dpu_* methods for getting subset filenames and parsing command line parameters
source(here::here("code", "depmap_pathways_util.R"))

sleep_seconds_between_requests <- 0.5
focused_lib <- c("Achilles_fitness_decrease", "Achilles_fitness_increase", "Aging_Perturbations_from_GEO_down", "Aging_Perturbations_from_GEO_up", "Allen_Brain_Atlas_down", "Allen_Brain_Atlas_up", "ARCHS4_Cell-lines", "ARCHS4_IDG_Coexp", "ARCHS4_Kinases_Coexp", "ARCHS4_TFs_Coexp", "ARCHS4_Tissues", "BioCarta_2016", "BioPlex_2017", "Cancer_Cell_Line_Encyclopedia", "ChEA_2016", "Chromosome_Location_hg19", "CORUM", "Data_Acquisition_Method_Most_Popular_Genes", "Disease_Perturbations_from_GEO_down", "Disease_Perturbations_from_GEO_up", "Disease_Signatures_from_GEO_up_2014", "Drug_Perturbations_from_GEO_down", "Drug_Perturbations_from_GEO_up", "DrugMatrix", "DSigDB", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "ENCODE_Histone_Modifications_2015", "ENCODE_TF_ChIP-seq_2015", "Enrichr_Libraries_Most_Popular_Genes", "Enrichr_Submissions_TF-Gene_Coocurrence", "Epigenomics_Roadmap_HM_ChIP-seq", "ESCAPE", "GeneSigDB", "GO_Biological_Process_2018", "GO_Cellular_Component_2018", "GO_Molecular_Function_2018", "GTEx_Tissue_Sample_Gene_Expression_Profiles_down", "GTEx_Tissue_Sample_Gene_Expression_Profiles_up", "GWAS_Catalog_2019", "HMDB_Metabolites", "HomoloGene", "Human_Gene_Atlas", "Human_Phenotype_Ontology", "HumanCyc_2015", "HumanCyc_2016", "huMAP", "InterPro_Domains_2019", "Jensen_COMPARTMENTS", "Jensen_DISEASES", "Jensen_TISSUES", "KEA_2015", "KEGG_2019_Human", "KEGG_2019_Mouse", "Kinase_Perturbations_from_GEO_down", "Kinase_Perturbations_from_GEO_up", "Ligand_Perturbations_from_GEO_down", "Ligand_Perturbations_from_GEO_up", "LINCS_L1000_Chem_Pert_down", "LINCS_L1000_Chem_Pert_up", "LINCS_L1000_Kinase_Perturbations_down", "LINCS_L1000_Kinase_Perturbations_up", "LINCS_L1000_Ligand_Perturbations_down", "LINCS_L1000_Ligand_Perturbations_up", "MCF7_Perturbations_from_GEO_down", "MCF7_Perturbations_from_GEO_up", "MGI_Mammalian_Phenotype_Level_4_2019", "Microbe_Perturbations_from_GEO_down", "Microbe_Perturbations_from_GEO_up", "miRTarBase_2017", "Mouse_Gene_Atlas", "MSigDB_Computational", "MSigDB_Oncogenic_Signatures", "NCI-60_Cancer_Cell_Lines", "NURSA_Human_Endogenous_Complexome", "Old_CMAP_down", "Old_CMAP_up", "OMIM_Disease", "OMIM_Expanded", "Panther_2016", "Pfam_Domains_2019", "Pfam_InterPro_Domains", "Phosphatase_Substrates_from_DEPOD", "PPI_Hub_Proteins", "Rare_Diseases_AutoRIF_ARCHS4_Predictions", "Rare_Diseases_AutoRIF_Gene_Lists", "Rare_Diseases_GeneRIF_ARCHS4_Predictions", "Rare_Diseases_GeneRIF_Gene_Lists", "Reactome_2016", "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO", "SILAC_Phosphoproteomics", "Single_Gene_Perturbations_from_GEO_down", "Single_Gene_Perturbations_from_GEO_up", "SubCell_BarCode", "SysMyo_Muscle_Gene_Sets", "TargetScan_microRNA_2017", "TF_Perturbations_Followed_by_Expression", "TF-LOF_Expression_from_GEO", "Tissue_Protein_Expression_from_Human_Proteome_Map", "Tissue_Protein_Expression_from_ProteomicsDB", "Transcription_Factor_PPIs", "TRANSFAC_and_JASPAR_PWMs", "TRRUST_Transcription_Factors_2019", "UK_Biobank_GWAS", "Virus_Perturbations_from_GEO_down", "Virus_Perturbations_from_GEO_up", "VirusMINT", "WikiPathways_2019_Human", "WikiPathways_2019_Mouse")

#define pathway enrichment analysis loop function
enrichr_loop <- function(gene_list, databases){
  if(is_empty(gene_list)){
    return(as_tibble())
  } else {
    flat_complete <- as_tibble()
    enriched <- enrichr(gene_list, databases)
    flat_complete <- bind_rows(enriched, .id = "enrichr")
    flat_complete <- flat_complete %>% 
      arrange(Adjusted.P.value) 
    flat_complete$enrichr <- str_replace_all(flat_complete$enrichr, "\\_", " ")
    flat_complete$Term <- str_replace_all(flat_complete$Term, "\\_", " ")
    return(flat_complete)
  }
}

get_gene_names_to_process <- function(subset_file_idx, num_subset_files, achilles_cor) {
  r <- "rowname" #need to drop "rowname"
  full <- (names(achilles_cor))[!(names(achilles_cor)) %in% r] #f[!f %in% r]
  grouped_gene_names <- split(full, ntile(full, num_subset_files))
  grouped_gene_names[[subset_file_idx]]  
}

read_input_data <- function() {
  load(file = here::here("data", "gene_summary.RData"))
  load(file = here::here("data", paste0(release, "_achilles_cor.RData")))
  achilles_lower <- readRDS(file = here::here("data", "achilles_lower.Rds"))
  achilles_upper <- readRDS(file = here::here("data", "achilles_upper.Rds"))
  
  #convert cor
  class(achilles_cor) <- c("cor_df", "tbl_df", "tbl", "data.frame") 
  
  list(
    gene_summary=gene_summary,
    achilles_cor=achilles_cor,
    achilles_lower=achilles_lower,
    achilles_upper=achilles_upper
  )
}

generate_positive_data <- function(gene_group, achilles_cor, achilles_upper, gene_summary) {
  #setup containers
  subset_positive <- tibble(
    fav_gene = character(), 
    data = list()
  )

  for (fav_gene in gene_group) {
    message("Top pathway enrichment analyses for ", fav_gene)
    flat_top_complete <- achilles_cor %>%
      focus(fav_gene) %>%
      arrange(desc(.[[2]])) %>% #use column index
      filter(.[[2]] > achilles_upper) %>% #formerly top_n(20), but changed to mean +/- 3sd
      rename(approved_symbol = rowname) %>%
      left_join(gene_summary, by = "approved_symbol") %>%
      select(approved_symbol, approved_name, fav_gene) %>%
      rename(gene = approved_symbol, name = approved_name, r2 = fav_gene) %>%
      pull("gene") %>%
      c(fav_gene, .) %>%
      enrichr_loop(., focused_lib) %>%
      arrange(Adjusted.P.value) %>%
      slice(1:100)

    positive <- flat_top_complete  %>%
      mutate(fav_gene = fav_gene) %>%
      group_by(fav_gene) %>%
      nest()

    subset_positive <- subset_positive %>%
      bind_rows(positive)

    Sys.sleep(sleep_seconds_between_requests)
  }

  subset_positive
}

generate_negative_data <- function(gene_group, achilles_cor, achilles_lower, gene_summary) {
  subset_negative <- tibble(
    fav_gene = character(), 
    data = list()
  )
  
  for (fav_gene in gene_group) {
    message("Bottom pathway enrichment analyses for ", fav_gene)
    flat_bottom_complete <- achilles_cor %>%
      focus(fav_gene) %>%
      arrange(.[[2]]) %>% #use column index
      filter(.[[2]] < achilles_lower) %>% #formerly top_n(20), but changed to mean +/- 3sd
      rename(approved_symbol = rowname) %>%
      left_join(gene_summary, by = "approved_symbol") %>%
      select(approved_symbol, approved_name, fav_gene) %>%
      rename(gene = approved_symbol, name = approved_name, r2 = fav_gene) %>%
      pull("gene") %>%
      enrichr_loop(., focused_lib) %>%
      arrange(Adjusted.P.value) %>%
      slice(1:100)

    negative <- flat_bottom_complete %>%
      mutate(fav_gene = fav_gene) %>%
      group_by(fav_gene) %>%
      nest()

    subset_negative <- subset_negative %>%
      bind_rows(negative)

    Sys.sleep(sleep_seconds_between_requests)
  }

  subset_negative
}

save_subset_pathways_file <- function(pathways_type, subset_file_idx, num_subset_files, subset_data) {
  subset_filepath <- dpu_get_subset_filepath(pathways_type, subset_file_idx, num_subset_files)
  message("Saving subset ", pathways_type, " data to ", subset_filepath)
  saveRDS(subset_data, file=subset_filepath)
}

# writes subset file into /data directory based on the pathway type, index and number subset files
create_pathway_subset_file <- function(pathways_type, subset_file_idx, num_subset_files) {
  message("Generating data for pathways ", pathways_type, " subset ", subset_file_idx, ".")
  input_data <- read_input_data()
  gene_group <- get_gene_names_to_process(subset_file_idx, num_subset_files, input_data$achilles_cor)
  message("Processing ", length(gene_group), " genes.")
  if (pathways_type == dpu_pathways_positive_type) {
    subset_data <- generate_positive_data(gene_group, input_data$achilles_cor, input_data$achilles_upper, input_data$gene_summary)
  } else {
    subset_data <- generate_negative_data(gene_group, input_data$achilles_cor, input_data$achilles_lower, input_data$gene_summary)
  }
  save_subset_pathways_file(pathways_type, subset_file_idx, num_subset_files, subset_data)
}

main <- function () {
  opt <- dpu_parse_command_line(include_idx=TRUE)
  start_time <- Sys.time()
  create_pathway_subset_file(opt$pathways_type, opt$subset_file_idx, opt$num_subset_files)
  end_time <- Sys.time()
  message("Elapsed time ", (end_time - start_time))
}

if(!interactive()) {
  main()
}


