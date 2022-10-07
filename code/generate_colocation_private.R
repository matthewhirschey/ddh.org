source(here::here("code", "current_release.R"))

#LOAD data 
subcell <- readRDS(file = here::here("data", paste0(release, "_subcell.Rds")))
expression_meta <- readRDS(file=here::here("data", paste0(release, "_expression_meta.Rds")))
expression_long <- readRDS(file=here::here("data", paste0(release, "_expression_long.Rds")))
TissueData <- readRDS(file = here::here("data", paste0(release, "_TissueData.Rds")))

#SUBCELL
#consider logistic regression because of so many variables to consider?

#EXPRESSION Correlation
# %>%
#   filter(!is.na(gene_expression),
#          !is.na(protein_expression)) #filter NAs? 26M filtered to 3M

# #cell lines in proteins_ccle  ##(what will eventually become) rows are filtered on semi-join
# # proteins_cell_lines <- 
# #   proteins_ccle_long %>% 
# #   dplyr::distinct(X1) %>% 
# #   pull(.)
# # 
# # genes_cell_lines <-
# #   expression_long %>% 
# #   dplyr::distinct(X1) %>% 
# #   pull(.)
# # 
# # expression_long_sub <- 
# #   expression_long %>% 
# #   filter(X1 %in% proteins_cell_lines)
# # 
# # proteins_ccle_long_sub <-
# #   proteins_ccle_long %>% 
# #   filter(X1 %in% genes_cell_lines)
# 
# #genes in proteins_ccle
# proteins_genes <-
#   proteins_ccle_long %>%
#   dplyr::distinct(gene_symbol) %>%
#   pull(.)
# 
# expression_genes <-
#   expression_long %>%
#   dplyr::distinct(gene_symbol) %>%
#   pull(.)
# 
# expression_long_sub <-
#   expression_long %>%
#   filter(gene_symbol %in% proteins_genes)
# 
# proteins_ccle_long_sub <-
#   proteins_ccle_long %>%
#   filter(gene_symbol %in% expression_genes)
# 
# #transpose to make same dim for correlation
# expression_t <- #continue with transposition to make df for corr
#   expression_long_sub %>% #make sure to use subset of data
#   tidyr::pivot_wider(names_from = gene_symbol, values_from = gene_expression) %>%  #names_glue = "{X1}_gene"
#   dplyr::arrange(X1) #switch names_from = X1, arrange(gene_symbol)
# 
# proteins_ccle_t <-
#   proteins_ccle_long_sub %>%
#   tidyr::pivot_wider(names_from = gene_symbol, values_from = protein_expression) %>%  #names_glue = "{X1}_protein",
#   dplyr::arrange(X1) #switch names_from = X1, arrange(gene_symbol)
# 
# expression_t <- #filter for same num rows
#   expression_t %>% 
#   dplyr::semi_join(proteins_ccle_t, by = "X1") #switch from "gene_symbol"
# 
# proteins_ccle_t <- #filter for same num rows
#   proteins_ccle_t %>% 
#   dplyr::semi_join(expression_t, by = "X1") #ibid
# 
# expression_t <- tibble::column_to_rownames(expression_t, var = "X1") #switch from "gene_symbol"
# proteins_ccle_t <- tibble::column_to_rownames(proteins_ccle_t, var = "X1") #ibid
# 
# dim(expression_t)
# dim(proteins_ccle_t)
# 
# ccle_gp_cor <-
#    corrr::correlate(x = proteins_ccle_t, y = expression_t)
# # gah...this gave me cell line v. cell line b/c of how I pivoted, so had to 'switch' some variables above
# 
# ccle_gp_cor %>% 
#   focus("SDHB") %>% 
#   arrange(desc(.[2]))

# Tissue correlation


#save files
# saveRDS(proteins_ccle, file = here::here("data", paste0(release, "_proteins_ccle.Rds")))
# saveRDS(ccle_gp_plot, file = here::here("data", paste0(release, "_ccle_gp_plot.Rds")))
# saveRDS(proteins_ccle_bygene, file = here::here("data", paste0(release, "_proteins_ccle_bygene.Rds")))
# saveRDS(proteins_ccle_bycell, file = here::here("data", paste0(release, "_proteins_ccle_bycell.Rds")))
