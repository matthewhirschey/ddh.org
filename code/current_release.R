release <- "20Q1"
achilles_url <- "https://ndownloader.figshare.com/files/21521910" #achilles_gene_effect.csv
ccle_url <- "https://ndownloader.figshare.com/files/21521940" #CCLE_expression.csv
cclemeta_url <- "https://ndownloader.figshare.com/files/21781221" #sample_info.csv
fraction_cutoff <- 0.05 #~5% FDR
na_cutoff_file = here::here("data", paste0(release, "_na_cutoff.Rds"))
if (file.exists(na_cutoff_file)) {
  na_cutoff <- readRDS(file = na_cutoff_file)   
}

pubtator_url <- "ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/gene2pubtatorcentral.gz"
gene2pubmed_url <- "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz"

