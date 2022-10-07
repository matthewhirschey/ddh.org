#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

generate_gene_location <- function(genes = NULL, regenerate = TRUE) {
  data("chromosome")
  if(regenerate == FALSE) {
    #get Rds from zenodo for protein information
    tmp_gene_location <- tempfile()
    download.file(gene_location_url, tmp_gene_location)
    gene_location <- readRDS(tmp_gene_location)
  } else {
fetch_biomart_data <- function(genes = NULL, regenerate = TRUE){
  # fetch_biomart_data fun based on John Blischak 2013-09-16
  # https://github.com/jdblischak/2013-09-19-chicago/blob/gh-pages/lessons/uchicago-r/data/create_datasets.R
  
  # Fetches Ensembl Biomart data for translated human genes.
  #
  # Args
  #  genes: character vector of ensembl gene IDs.
  #         If NULL, all genes are fetched.
  #
  # Results
  #  A data.frame where each row is the Ensembl transcipt with the
  #  longest coding sequence for a given Ensembl gene ID.
  #
  ensembl <- useMart("ensembl", 
                     dataset = "hsapiens_gene_ensembl", 
                     host = "useast.ensembl.org")
  #   listFilters(ensembl) #filters are the input, in this case ensembl_gene_id
  #   listAttributes(ensembl, "feature_page")
  #   attributePages(ensembl)
  stopifnot(is.null(genes) || is.vector(genes, mode = "character"))
  atts <- c("ensembl_transcript_id", "ensembl_gene_id", "chromosome_name", 
            "transcript_start", "transcript_end", "strand", "band", 
            "external_gene_name", "percentage_gene_gc_content", "gene_biotype",
            "source", "name_1006", "description")
  if (is.null(genes)) {
    gene_info <- getBM(attributes = atts, mart = ensembl)
  } else {
    gene_info <- getBM(attributes = atts, filters = "ensembl_gene_id",
                       values = genes, mart = ensembl)
  }
  gene_coding_len <- getBM(attributes = c("ensembl_transcript_id", "cds_length"), #dat on different 'pages' so need it sep
                           filters = "ensembl_transcript_id",
                           values = gene_info$ensembl_transcript_id,
                           mart = ensembl)
  gene_coding_len <- na.omit(gene_coding_len)
  gene_final <- 
    gene_info %>% 
    dplyr::left_join(gene_coding_len, by = "ensembl_transcript_id") %>% 
    dplyr::arrange(desc(ensembl_gene_id, cds_length)) %>% 
    dplyr::distinct(ensembl_gene_id, .keep_all = TRUE) %>%
    dplyr::rename(approved_symbol = external_gene_name) #change name here to it matches gene_summary and works with fun_text.R
  
  return(gene_final)
}  

# Testing:
# g <- c("ENSG00000198888", "ENSG00000198763", "ENSG00000198804",
#        "ENSG00000198712", "ENSG00000228253", "ENSG00000198899", "ENSG00000012048")
# 
# test_data <- fetch_biomart_data(g)

# Create full dataset
gene_location <- fetch_biomart_data()

# filter
gene_location <- 
  gene_location %>% 
  filter(chromosome_name %in% chromosome$id)
}

#save files
saveRDS(gene_location, file = here::here("data", paste0(release, "_gene_location.Rds")))
saveRDS(chromosome, file = here::here("data", paste0(release, "_chromosome.Rds")))
}
