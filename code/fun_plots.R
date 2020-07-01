library(tidyverse)
library(cowplot)
library(viridis)
library(plotly)
library(gganatogram)

make_cellbins <- function(cellbins_data = achilles, expression_data = expression_join, gene_symbol) {
  plot_complete <- cellbins_data %>% #plot setup
    select(X1, any_of(gene_symbol)) %>%
    left_join(expression_data, by = "X1") %>%
    select(-X1) %>%
    pivot_longer(cols = !any_of(c("cell_line", "lineage", "lineage_subtype")), names_to = "gene_symbol", values_to = "dep_score") %>% #replace any_of with where(is.numeric)
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
    plot_complete  <- plot_complete +
      guides(fill = "none")
  } else {
    plot_complete
  }
  return(plot_complete)
}

#figure legend
plot_cellbins_title <- "Kernel density estimate."
plot_cellbins_legend <- "A smoothed version of the histogram of Dependency Scores. Dependency scores across all cell lines for queried genes, revealing overall influence of a gene on cellular fitness"

make_celldeps <- function(celldeps_data = achilles, expression_data = expression_join, gene_symbol, mean) {
  plot_complete <- celldeps_data %>% #plot setup
    select(X1, any_of(gene_symbol)) %>%
    left_join(expression_data, by = "X1") %>%
    select(-X1) %>%
    pivot_longer(cols = !any_of(c("cell_line", "lineage", "lineage_subtype")), names_to = "gene_symbol", values_to = "dep_score") %>% 
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
    plot_complete  <- plot_complete +
      guides(color = "none")
  } else {
    plot_complete
  }
  return(plot_complete)
}

#figure legend
plot_celldeps_title <- "Cell Line Dependency Curve."
plot_celldeps_legend <- "Each point shows the ranked dependency score for a given cell line. Cells with dependency scores less than -1 indicate a cell that the query gene is essential within. Cells with dependency scores close to 0 show no changes in fitness when the query gene is knocked out. Cells with dependency scores greater than 1 have a gain in fitness when the query gene is knocked-out"

# make cell anatogram
make_cellanatogram <- function(cellanatogram_data = subcell, gene_symbol) {
  plot_complete <- cellanatogram_data %>% 
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
    plot_complete  <- plot_complete +
      guides(fill = "none")
  } else {
    plot_complete
  }
  return(plot_complete)
}

# make lineage plot
make_lineage <- function(celldeps_data = achilles, expression_data = expression_join, gene_symbol) {
  data_full <- celldeps_data %>% #plot setup
    select(X1, any_of(gene_symbol)) %>%
    left_join(expression_data, by = "X1") %>%
    select(-X1) %>%
    pivot_longer(cols = !any_of(c("cell_line", "lineage", "lineage_subtype")), names_to = "gene_symbol", values_to = "dep_score") %>%
    dplyr::mutate_at("lineage", function(str) {
      str <- str_replace_all(str, "\\_", " ")
      str <- str_to_title(str)
      return(str)
    }) %>% 
    drop_na(lineage) %>% 
    drop_na(dep_score)
  
  data_mean <- data_full %>% 
    group_by(lineage) %>% 
    summarize(dep_score = mean(dep_score))
  
  plot_complete <- ggplot() +
    geom_boxplot(data = data_full, aes(x = fct_reorder(lineage, dep_score, .fun = mean, .desc = TRUE), 
                     y = dep_score
                     )) +
    geom_point(data = data_mean, aes(x = lineage, y = dep_score), color = "blue", alpha = 0.5) +
    coord_flip() +
    labs(y = "Dependency Score", x = "Lineage") +
    theme_minimal_vgrid()
  return(plot_complete)
}

#figure legend
plot_celllin_title <- "Cell Line Lineage Dependencies"
plot_celllin_legend <- "Each point shows the mean dependency score for a given cell lineage, with box plots showing median, interquartile ranges, and outliers."

# make sublineage plot
make_sublineage <- function(celldeps_data = achilles, expression_data = expression_join, gene_symbol) {
  data_full <- celldeps_data %>% #plot setup
    select(X1, any_of(gene_symbol)) %>%
    left_join(expression_data, by = "X1") %>%
    select(-X1) %>%
    pivot_longer(cols = is.numeric, names_to = "gene_symbol", values_to = "dep_score") %>% 
    dplyr::mutate_at("lineage_subtype", function(str) {
      str <- str_replace_all(str, "\\_", " ")
      str <- if_else(str_detect(str, "^[:lower:]"), str_to_title(str), str)
      return(str)
    })  %>% 
    drop_na(lineage_subtype) %>% 
    drop_na(dep_score)
  
  data_mean <- data_full %>% 
    group_by(lineage_subtype) %>% 
    summarize(dep_score = mean(dep_score))
  
  plot_complete <- ggplot() +
    geom_boxplot(data = data_full, aes(x = fct_reorder(lineage_subtype, dep_score, .fun = mean, .desc = TRUE), 
                               y = dep_score
    )) +
    geom_point(data = data_mean, aes(x = lineage_subtype, y = dep_score), color = "lightblue", alpha = 0.5) +
    coord_flip() +
    labs(y = "Dependency Score", x = "Sublineage") +
    theme_minimal_vgrid()
  return(plot_complete)
}

#figure legend
plot_cellsublin_title <- "Cell Line Sub-Lineage Dependencies"
plot_cellsublin_legend <- "Each point shows the mean dependency score for a given cell sublineage, with box plots showing median, interquartile ranges, and outliers."

