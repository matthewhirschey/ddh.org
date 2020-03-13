library(tidyverse)
library(cowplot)
library(viridis)


make_cellbins <- function(gene_symbol) {
  achilles %>% #plot setup
    select(X1, gene_symbol) %>%
    left_join(expression_join, by = "X1") %>%
    rename(dep_score = gene_symbol) %>%
    select(cell_line, lineage, dep_score) %>%
    arrange(dep_score) %>%
    ggplot() +
    geom_vline(xintercept = 1, color = "lightgray") +
    geom_vline(xintercept = -1, color = "lightgray") +
    geom_vline(xintercept = 0) +
    geom_histogram(aes(x = dep_score), binwidth = 0.25, color = "gray", fill = "#02224C") +
    labs(x = "Dependency Score (binned)") +
    theme_cowplot()
}

make_celldeps <- function(gene_symbol) {
  achilles %>% #plot setup
    select(X1, gene_symbol) %>%
    left_join(expression_join, by = "X1") %>%
    rename(dep_score = gene_symbol) %>%
    select(cell_line, lineage, dep_score) %>%
    arrange(dep_score) %>%
    ggplot() +
    geom_point(aes(x = fct_reorder(cell_line, dep_score, .desc = FALSE), 
                   y = dep_score, 
                   text = paste0("Cell Line: ", cell_line)), 
               alpha = 0.2, color = "#02224C") +
    labs(x = "Cell Lines", y = "Dependency Score") +
    geom_hline(yintercept = mean_virtual_achilles) +
    geom_hline(yintercept = 1, color = "lightgray") +
    geom_hline(yintercept = -1, color = "lightgray") +
    geom_hline(yintercept = 0) +
    scale_x_discrete(expand = expansion(mult = 0.02), na.translate = FALSE) +
    theme_cowplot() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + # axis.title.x=element_blank()
    NULL
}