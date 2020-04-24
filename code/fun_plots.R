library(tidyverse)
library(cowplot)
library(viridis)
library(plotly)

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
  return(ggplotly(p, tooltip = c("text")))
}

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
