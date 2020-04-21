library(tidyverse)
library(cowplot)
library(viridis)

gene_symbol <- "SDHA"
gene_symbol <- c("SDHA", "SDHB", "SDHC", "SDHD")
make_cellbins(achilles, expression_join, gene_symbol)
make_celldeps(achilles, expression_join, gene_symbol, mean_virtual_achilles)

make_cellbins <- function(table, expression_table, gene_symbol) {
  p <- table %>% #plot setup
    select(X1, any_of(gene_symbol)) %>%
    left_join(expression_table, by = "X1") %>%
    select(-X1) %>%
    pivot_longer(cols = !any_of(c("cell_line", "lineage")), names_to = "gene_symbol", values_to = "dep_score") %>% 
    ggplot() +
    geom_vline(xintercept = 1, color = "lightgray") +
    geom_vline(xintercept = -1, color = "lightgray") +
    geom_vline(xintercept = 0) +
    geom_histogram(aes(x = dep_score, fill = gene_symbol), binwidth = 0.25, color = "gray") + #fill = "#02224C"
    labs(x = "Dependency Score (binned)") +
    scale_fill_viridis(discrete = TRUE, option = "D", direction = 1, guide = "legend") +
    labs(fill = "Query Gene") +
    theme_cowplot()
  
  if(length(gene_symbol) == 1){
    p  <- p +
      guides(fill = "none")
  } else {
    p
  }
  return(p)
}

make_celldeps <- function(data_table, expression_table, gene_symbol, mean) {
  p <- data_table %>% #plot setup
    select(X1, any_of(gene_symbol)) %>%
    left_join(expression_table, by = "X1") %>%
    select(-X1) %>%
    pivot_longer(cols = !any_of(c("cell_line", "lineage")), names_to = "gene_symbol", values_to = "dep_score") %>% 
    ggplot(aes(x = fct_reorder(cell_line, dep_score, .desc = FALSE), 
               y = dep_score, 
               text = paste0("Cell Line: ", cell_line), 
               color = gene_symbol)) +
    geom_point(alpha = 0.4) + #color = "#02224C"
    labs(x = "Cell Lines", y = "Dependency Score", color = "Query Gene") +
    geom_hline(yintercept = mean) +
    geom_hline(yintercept = 1, color = "lightgray") +
    geom_hline(yintercept = -1, color = "lightgray") +
    geom_hline(yintercept = 0) +
    scale_color_viridis(discrete = TRUE, guide = "legend") +
    scale_x_discrete(expand = expansion(mult = 0.02), na.translate = FALSE) +
    theme_cowplot() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + # axis.title.x=element_blank()
    NULL
  
  if(length(gene_symbol) == 1){
    p  <- p +
      guides(color = "none")
  } else {
    p
  }
  return(p)
}

#data_table <- achilles
#expression_table <- expression_join
#gene_symbol <- "SDHA"
#gene_symbol <- c("SDHA", "SDHB", "SDHC", "SDHD")
#mean <- mean_virtual_achilles

