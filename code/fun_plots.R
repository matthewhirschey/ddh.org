library(tidyverse)
library(cowplot)
library(viridis)


make_cellbins <- function(data_table, expression_table, gene_symbol) {
  data_table %>% #plot setup
    select(X1, all_of(gene_symbol)) %>%
    left_join(expression_table, by = "X1") %>%
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

make_celldeps <- function(data_table, expression_table, gene_symbol, mean) {
  data_table %>% #plot setup
    select(X1, all_of(gene_symbol)) %>%
    left_join(expression_table, by = "X1") %>%
    pivot_longer(cols = !any_of(c("X1", "cell_line", "lineage")), names_to = "gene_symbol", values_to = "dep_score") %>% 
    select(-X1) %>%
    arrange(dep_score) %>%
    ggplot(aes(x = fct_reorder(cell_line, dep_score, .desc = FALSE), 
               y = dep_score, 
               text = paste0("Cell Line: ", cell_line), 
               color = gene_symbol)) +
    geom_point(alpha = 0.2) + #color = "#02224C"
    labs(x = "Cell Lines", y = "Dependency Score") +
    geom_hline(yintercept = mean) +
    geom_hline(yintercept = 1, color = "lightgray") +
    geom_hline(yintercept = -1, color = "lightgray") +
    geom_hline(yintercept = 0) +
    scale_color_viridis(discrete = TRUE) +
    scale_x_discrete(expand = expansion(mult = 0.02), na.translate = FALSE) +
    theme_cowplot() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + # axis.title.x=element_blank()
    labs(color = "") +
    NULL
}

#data_table <- achilles
#expression_table <- expression_join
#gene_symbol <- "SDHA"
#gene_symbol <- c("SDHA", "SDHB", "SDHC", "SDHD")
#mean <- mean_virtual_achilles

