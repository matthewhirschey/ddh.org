library(tidyverse)
library(cowplot)
library(viridis)
library(plotly)
library(gganatogram)

make_cellbins <- function(data_table, expression_table, gene_symbol) {
  p <- data_table %>% #plot setup
    select(X1, any_of(gene_symbol)) %>%
    left_join(expression_table, by = "X1") %>%
    select(-X1) %>%
    pivot_longer(cols = !any_of(c("cell_line", "lineage")), names_to = "gene_symbol", values_to = "dep_score") %>% 
    ggplot() +
    geom_vline(xintercept = 1, color = "lightgray") +
    geom_vline(xintercept = -1, color = "lightgray") +
    geom_vline(xintercept = 0) +
    geom_density(aes(x = dep_score, 
                     fill = fct_reorder(gene_symbol, dep_score, .fun = median), 
                     text = paste0("Gene: ", gene_symbol)), alpha = 0.6, color = "gray") +
    labs(x = "Dependency Score (binned)", fill = "Query \nGene", y = "Density") +
    scale_fill_viridis(discrete = TRUE, option = "D", direction = 1, guide = "legend") +
    theme_cowplot()
  
  if(length(gene_symbol) == 1){
    p  <- p +
      guides(fill = "none")
  } else {
    p
  }
  return(p)
}

#figure legend
plot_cellbins_title <- "Kernel density estimate."
plot_cellbins_legend <- "A smoothed version of the histogram of Dependency Scores. Dependency scores across all cell lines for queried genes, revealing overall influence of a gene on cellular fitness"

make_celldeps <- function(data_table, expression_table, gene_symbol, mean) {
  p <- data_table %>% #plot setup
    select(X1, any_of(gene_symbol)) %>%
    left_join(expression_table, by = "X1") %>%
    select(-X1) %>%
    pivot_longer(cols = !any_of(c("cell_line", "lineage")), names_to = "gene_symbol", values_to = "dep_score") %>% 
    ggplot(aes(x = fct_reorder(cell_line, dep_score, .fun = max, .desc = FALSE), 
               y = dep_score, 
               text = paste0("Cell Line: ", cell_line), 
               color = fct_reorder(gene_symbol, dep_score, .fun = max, .desc = TRUE)
               )) +
    geom_point(alpha = 0.4) + #color = "#02224C"
    #geom_boxplot() +
    labs(x = "Cell Lines", y = "Dependency Score", color = "Query \nGene") +
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

#figure legend
plot_celldeps_title <- "Cell Line Dependency Curve."
plot_celldeps_legend <- "Each point shows the ranked dependency score for a given cell line. Cells with dependency scores less than -1 indicate a cell that the query gene is essential within. Cells with dependency scores close to 0 show no changes in fitness when the query gene is knocked out. Cells with dependency scores greater than 1 have a gain in fitness when the query gene is knocked-out"

# make cell anatogram
make_cellanatogram <- function(data_table = subcell, gene_symbol) {
  p <- data_table %>% 
    filter_all(any_vars(gene_name %in% gene_symbol)) %>% 
    filter(!is.na(type)) %>% 
    select(-value) %>% 
    add_count(main_location) %>% 
    mutate(value = as_factor(n)) %>% 
    gganatogram(outline = T, fillOutline='grey95', organism="cell", fill = "value") +
    theme_void() +  
    coord_fixed() +
    scale_fill_viridis(discrete = TRUE) +
    labs(fill = "Count")
  
  if(length(gene_symbol) == 1){
    p  <- p +
      guides(fill = "none")
  } else {
    p
  }
  return(p)
}
