#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

generate_threshold <- function() {
achilles <- 
  readr::read_csv(achilles_url, col_names = TRUE) %>% 
  `colnames<-`(str_remove_all(names(.), "\\s\\(\\d+\\)")) %>% 
  dplyr::rename(X1 = 1)

expression <- 
  readr::read_csv(ccle_url, col_names = TRUE) %>% 
  `colnames<-`(str_remove_all(names(.), "\\s\\(\\d+\\)")) %>%
  dplyr::rename(X1 = 1)

achilles_long <- 
  achilles %>% 
  tidyr::pivot_longer(-X1, names_to = "gene", values_to = "dep_score")

#filter achilles to remove no expression dep scores(special sauce)
expression_long <- 
  expression %>% 
  dplyr::filter(expression$X1 %in% achilles$X1 == TRUE) %>% #matches cells
  tidyr::gather("gene", "gene_expression", -X1) %>% 
  dplyr::arrange(desc(gene_expression))

no_expression <- 
  expression_long %>% 
  dplyr::filter(gene_expression == 0) %>% 
  tidyr::unite(X1, gene, col = "match", sep = "-", remove = TRUE) %>% 
  dplyr::pull(match)

achilles_no0 <- 
  achilles_long %>% 
  tidyr::unite(X1, gene, col = "match", sep = "-", remove = FALSE) %>% 
  dplyr::filter(match %in% no_expression == FALSE) %>% 
  dplyr::select(-match) %>%
  tidyr::spread(gene, dep_score)

achilles_no0_plot <- 
  achilles_no0 %>% 
  dplyr::summarise_all(list(~sum(is.na(.)))) %>% 
  tidyr::gather(gene, NAs) %>% 
  dplyr::arrange(desc(NAs)) %>% 
  dplyr::mutate(pos = sum(achilles$X1 %in% expression$X1)-NAs)

#skip toomanyNAs setps

#clean Achilles correlation matrix
achilles_cor <- 
  achilles_no0 %>%
  dplyr::select(-X1) %>% 
  corrr::correlate() #(diagonal = 0) set to 0 so easy to summarize, but should be NA; so added na.rm = TRUE to fun() in EDA

#make some long files
achilles_cor_long <- 
  achilles_cor %>% 
  corrr::stretch()

#Permutation tests
virtual_achilles <- 
  achilles_cor_long %>% 
  dplyr::filter(!is.na(r)) %>%   
  moderndive::rep_sample_n(size = 20000, reps = 1000) %>%
  dplyr::group_by(replicate) %>% 
  dplyr::summarize(mean = mean(r), max = max(r), min = min(r), sd = sd(r))

mean_virtual_achilles <- mean(virtual_achilles$mean)
sd_virtual_achilles <- mean(virtual_achilles$sd)

sd_threshold <- 3

achilles_upper <- mean_virtual_achilles + sd_threshold*sd_virtual_achilles
achilles_lower <- mean_virtual_achilles - sd_threshold*sd_virtual_achilles

na_cutoff <- achilles_cor_long %>% 
  dplyr::filter(r > achilles_upper) %>%  #| achilles_correlation_raw < achilles_lower) %>% 
  dplyr::group_by(x) %>% 
  dplyr::summarize(count = n()) %>% 
  dplyr::left_join(achilles_no0_plot, by = c("x" = "gene")) %>% 
  dplyr::arrange(pos) %>% 
  dplyr::top_frac(-fraction_cutoff, pos) %>% 
  dplyr::arrange(pos) %>% 
  dplyr::pull(NAs) %>% 
  min(.)

print(na_cutoff)

#save
saveRDS(na_cutoff, file = here::here("data", paste0(release, "_na_cutoff.Rds")))
}
