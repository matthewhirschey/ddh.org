library(tidyverse)
library(cowplot)
library(colorspace)
library(ggdist)
library(plotly)
library(gganatogram)

#library(extrafont)

#extrafont::loadfonts()


## COLORS ----------------------------------------------------------------------

## colors are hard-coded for each of the three groups:
## genes, cells, drugs
source(here::here("code", "generate_colors.R"))

## three color gradient: light version, main color, dark version
## main color
## color_set_gene[2] 
color_set <- color_set_gene
main_color <- color_set_gene[2] 

## for n colors use:
## color_pal_gene(n)
## which returns n color scaled between color_set[1] and color_set[3],
## centered around color_set[2]

## DENSITY PLOT ----------------------------------------------------------------
make_cellbins <- function(cellbins_data = achilles, expression_data = expression_names, gene_symbol) {
  plot_complete <- cellbins_data %>% #plot setup
    select(X1, any_of(gene_symbol)) %>%
    left_join(expression_data, by = "X1") %>%
    select(-X1) %>%
    pivot_longer(
      cols = where(is.numeric), 
      names_to = "gene_symbol", 
      values_to = "dep_score"
    ) %>% 
    group_by(gene_symbol) %>% 
    arrange(dep_score) %>% 
    mutate(med = median(dep_score, na.rm= TRUE)) %>% 
    ungroup() %>% 
    filter(!is.na(dep_score)) %>% 
    ggplot() +
    ## annotation range -1 to 1
    geom_rect(
      xmin = -1, xmax = 1,
      ymin = -Inf, ymax = Inf,
      fill = "grey95"
    ) +
    ## indicator line y axis
    geom_linerange(
      aes(xmin = -Inf, xmax = med, 
          y = fct_reorder(gene_symbol, -med), 
          color = med < -1),
      linetype = "dotted",
      size = .2
    ) +
    ## density curves via {ggdist}
    stat_halfeye(aes(x = dep_score, y = fct_reorder(gene_symbol, -med),
                     fill = stat(abs(x) > 1),
                     point_fill = after_scale(fill)),
                 .width = c(.025, .975),
                 color = "black",
                 shape = 21,
                 stroke = .7,
                 point_size = 2) +
    ## zero line
    geom_vline(
      xintercept = 0, 
      color = "grey80", 
      linetype = "dashed"
    ) +
    ## titles
    labs(
      x = "Dependency Score (binned)", y = NULL, 
      color = "Query \nGene", fill = "Query \nGene"
    ) +
    ## scales + legends
    scale_y_discrete(expand = c(.03, .03)) +
    scale_color_manual(values = c("grey70", main_color)) +
    scale_fill_manual(values = c("grey70", main_color)) +
    guides(
      color = guide_legend(size = 1, reverse = TRUE),
      fill = guide_legend(size = 1, reverse = TRUE)
    ) +
    ## theme changes
    theme_cowplot(font_size = 16) +
    theme(
      text = element_text(family = "Nunito Sans"),
      legend.position = "none", 
      axis.line.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      axis.text = element_text(family = "Roboto Slab"),
      axis.text.x = element_text(size = 12, color = "grey30"),
      axis.title = element_text(size = 15)
    ) +
    NULL
    return(plot_complete)
}

#figure legend
plot_cellbins_title <- "Computed density curves."
plot_cellbins_legend <- "Kernel density estimate of dependency scores. Dependency scores across all cell lines for queried genes, revealing overall influence of a gene on cellular fitness. The interval indicates the 95% quantile of the data, the dot indicates the median dependency score. The gray background highlights weak dependency values between -1 and 1."


## DOT PLOT --------------------------------------------------------------------
make_celldeps <- function(celldeps_data = achilles, expression_data = expression_names, gene_symbol, mean) {
  
  ## use main color in case of 1 selected gene, otherwise use palette function
  if(length(gene_symbol) == 1) { 
    cols <- main_color 
  }else{
    cols <- color_pal_gene(length(gene_symbol))
  }
  
  plot_complete <- celldeps_data %>% #plot setup
    select(X1, any_of(gene_symbol)) %>%
    left_join(expression_data, by = "X1") %>%
    select(-X1) %>%
    pivot_longer(
      cols = where(is.numeric), 
      names_to = "gene_symbol", 
      values_to = "dep_score"
    ) %>% 
    group_by(gene_symbol) %>% 
    arrange(dep_score) %>% 
    mutate(
      rank = 1:n(),
      med = median(dep_score, na.rm= TRUE)
    ) %>% 
    ungroup() %>% 
    ggplot(aes(x = rank, 
               y = dep_score, 
               text = paste0("Gene: ", gene_symbol, "\nCell Line: ", cell_line), 
               color = fct_reorder(gene_symbol, med),
               fill = fct_reorder(gene_symbol, med)
    )) +
    ## indicator lines dep. score
    geom_hline(yintercept = mean) +
    geom_hline(yintercept = 1, size = .2, color = "grey70", linetype = "dashed") +
    geom_hline(yintercept = -1, size = .2, color = "grey70", linetype = "dashed") +
    geom_hline(yintercept = 0, size = .2, color = "grey50") +
    ## dot plot
    geom_point(size = 1.1, stroke = .25, alpha = 0.6) +
    ## scales + legends
    scale_x_discrete(expand = expansion(mult = 0.02), na.translate = FALSE) +
    scale_color_manual(values = cols, guide = "legend") +
    scale_fill_manual(values = cols, guide = "legend") +
    guides(
      color = guide_legend(reverse = TRUE, override.aes = list(size = 4, stroke = .8)),
      fill = guide_legend(reverse = TRUE, override.aes = list(size = 3.8, stroke = .8))
    ) +
    ## titles
    labs(
      x = NULL, y = paste0(gene_symbol, " Dependency Score"), 
      color = "Query Gene", fill = "Query Gene"
    ) +
    ## theme changes
    theme_cowplot(font_size = 12) +
    theme(
      text = element_text(family = "Nunito Sans"),
      axis.text = element_text(family = "Roboto Slab"),
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(), 
      axis.line.x = element_blank()
    ) + 
    NULL
  
  ##only plot legend in case of more than 1 gene selected
  if(length(gene_symbol) == 1){
    plot_complete  <- plot_complete +
      theme(legend.position = "none")
    plot_complete
  } else {
    plot_complete
  }
}

#figure legend
plot_celldeps_title <- "Cell Line Dependency Curve."
plot_celldeps_legend <- "Each point shows the ranked dependency score ordered from low to high scores. Cells with dependency scores less than -1 indicate a cell that the query gene is essential within. Cells with dependency scores close to 0 show no changes in fitness when the query gene is knocked out. Cells with dependency scores greater than 1 have a gain in fitness when the query gene is knocked-out."

## ANATOGRAM -------------------------------------------------------------------
# make cell anatogram
make_cellanatogram <- function(cellanatogram_data = subcell, gene_symbol) {
  plot_complete <- cellanatogram_data %>% 
    filter_all(any_vars(gene_name %in% gene_symbol)) %>% 
    filter(!is.na(type)) %>% 
    select(-value) %>% 
    add_count(main_location) %>% 
    mutate(value = as_factor(n)) %>% 
    gganatogram(outline= TRUE, fillOutline='grey95', organism="cell", fill = "value") +
    theme_void(base_size = 14) +  
    theme(
      text = element_text(family = "Nunito Sans"),
      plot.margin = margin(5, 10, 5, 5)
    ) +
    coord_fixed() +
    scale_fill_viridis_d() +
    labs(fill = "Count") +
    NULL
  
  if(length(gene_symbol) == 1){
    plot_complete  <- plot_complete +
      guides(fill = "none")
  } 
  return(plot_complete)
}


## LINEAGE LINERANGE PLOT ------------------------------------------------------
# make lineage plot
make_lineage <- function(celldeps_data = achilles, expression_data = expression_names, gene_symbol) {
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
    drop_na(dep_score) %>% 
    group_by(lineage) %>% 
    mutate(mean = mean(dep_score)) %>% 
    ungroup() %>% 
    mutate(lineage = fct_reorder(lineage, -mean))
  
  data_mean <- data_full %>% 
    dplyr::group_by(lineage) %>% 
    dplyr::summarize(dep_score = mean(dep_score))
  
  plot_complete <- 
    ggplot(data_mean, aes(dep_score, lineage)) +
    ## annotation range -1 to 1
    geom_rect(
      xmin = -1, xmax = 1,
      ymin = -Inf, ymax = Inf,
      fill = "grey95"
    ) +
    ## zero line
    geom_vline(
      xintercept = 0, 
      color = "grey80", 
      linetype = "dashed"
    ) +
    ## indicator lines lineages
    geom_linerange(
      aes(xmin = -Inf, xmax = dep_score),
      color = "grey60",
      linetype = "dotted"
    ) +
    ## lineranges as "boxplots"
    stat_interval(
      data = data_full,
      orientation = "horizontal",
      .width = c(.05, .5, .95)
    ) +
    ## dot indicating mean
    geom_point(
      color = "black", fill = "white", 
      shape = 21, stroke = .5, 
      size = 1.8
    ) +
    ## scales + legends
    scale_x_continuous(
      sec.axis = dup_axis()
    ) +
    scale_color_manual(
      values = color_set, 
      labels = c("95% of the data fall in these ranges", "50%", "5%"), 
      name = ""
    ) +
    guides(color = guide_legend(reverse = TRUE)) +
    ## titles
    labs(
      x = "Dependency Score", 
      y = NULL, 
      title = paste0("Cell lineage dependencies for ", gene_symbol)
    ) +
    ## theme changes
    theme_cowplot(font_size = 16) +
    theme(
      #text = element_text(family = "Nunito Sans"),
      legend.position = "top", 
      axis.line.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      axis.text = element_text(family = "Roboto Slab"),
      axis.text.x = element_text(size = 12, color = "grey30"),
      axis.title = element_text(size = 15),
      axis.title.x.bottom = element_blank(),
      plot.title = element_text(size = 17), 
      plot.title.position = "plot"
    ) +
    NULL
  return(plot_complete)
}

#figure legend
plot_celllin_title <- "Cell Line Lineage Dependencies."
plot_celllin_legend <- "Each point shows the mean dependency score for the gene query within a given cell lineage. The intervals show the 5% quantiles centered on the median, the interquartile ranges, and the 95% quantiles. The gray background highlights weak dependency values between -1 and 1."


## SUBLINE RANGE PLOT ---------------------------------------------------

# make sublineage plot
make_sublineage <- function(celldeps_data = achilles, expression_data = expression_names, gene_symbol) {
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
    drop_na(dep_score) %>% 
    group_by(lineage_subtype) %>% 
    mutate(mean = mean(dep_score)) %>% 
    ungroup() %>% 
    mutate(lineage_subtype = fct_reorder(lineage_subtype, -mean))
  
  data_mean <- data_full %>% 
    dplyr::group_by(lineage_subtype) %>% 
    dplyr::summarize(dep_score = mean(dep_score))
  
  plot_complete <- 
    ggplot(data_mean, aes(dep_score, lineage_subtype)) +
    ## annotation range -1 to 1
    geom_rect(
      xmin = -1, xmax = 1,
      ymin = -Inf, ymax = Inf,
      fill = "grey95"
    ) +
    ## zero line
    geom_vline(
      xintercept = 0, 
      color = "grey80", 
      linetype = "dashed"
    ) +
    ## indicator lines sublineages
    geom_linerange(
      aes(xmin = -Inf, xmax = dep_score),
      color = "grey60",
      linetype = "dotted"
    ) +
    ## lineranges as "boxplots"
    stat_interval(
      data = data_full,
      orientation = "horizontal",
      .width = c(.05, .5, .95)
    ) +
    ## dot indicating mean
    geom_point(
      color = "black", fill = "white", 
      shape = 21, stroke = .5, 
      size = 1.8
    ) +
    ## scales + legends
    scale_x_continuous(
      sec.axis = dup_axis()
    ) +
    scale_color_manual(
      values = color_set, 
      labels = c("95% of the data fall in these ranges", "50%", "5%"), 
      name = ""
    ) +
    guides(color = guide_legend(reverse = TRUE)) +
    ## titles
    labs(
      x = "Dependency Score", 
      y = NULL, 
      title = paste0("Cell sub-lineage dependencies for ", gene_symbol)
    ) +
    ## theme changes
    theme_cowplot(font_size = 16) +
    theme(
      #text = element_text(family = "Nunito Sans"),
      legend.position = "top", 
      axis.line.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      axis.text = element_text(family = "Roboto Slab"),
      axis.text.x = element_text(size = 12, color = "grey30"),
      axis.title = element_text(size = 15),
      axis.title.x.bottom = element_blank(),
      plot.title = element_text(size = 17), 
      plot.title.position = "plot"
    ) +
    NULL
  return(plot_complete)
}

#figure legend
plot_cellsublin_title <- "Cell Line Sub-Lineage Dependencies"
plot_cellsublin_legend <- "Each point shows the mean dependency score for the gene query within a given cell lineage. The intervals show the 5% quantiles centered on the median, the interquartile ranges, and the 95% quantiles. The gray background highlights weak dependency values between -1 and 1."

## CELL EXPRESSION PLOT --------------------------------------------------------
make_cellexpression <- function(expression_data = expression, expression_join = expression_names, gene_symbol, mean = mean_virtual_expression, upper_limit = expression_upper, lower_limit = expression_lower) {
  plot_complete <- expression_data %>% #plot setup
    dplyr::select(where(is.character), any_of(gene_symbol)) %>% 
    dplyr::left_join(expression_join, by = "X1") %>%
    dplyr::select(-X1) %>%
    dplyr::select(cell_line, lineage, lineage_subtype, everything()) %>% 
    dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>% 
    pivot_longer(cols = where(is.numeric), names_to = "gene_symbol", values_to = "cell_exp") %>% 
    ggplot(aes(x = fct_reorder(cell_line, cell_exp, .fun = max, .desc = FALSE), 
               y = cell_exp, 
               text = paste0("Cell Line: ", cell_line), 
               color = fct_reorder(gene_symbol, cell_exp, .fun = max, .desc = TRUE)
    )) +
    geom_point(alpha = 0.4) + #color = "#02224C"
    labs(x = "Cell Lines", y = paste0(gene_symbol, " Expression"), color = "Query \nGene") +
    geom_hline(yintercept = mean) +
    geom_hline(yintercept = upper_limit, color = "lightgray") + #3SD
    geom_hline(yintercept = 0, color = "lightgray") + #3SD is below zero
    scale_color_viridis(discrete = TRUE, guide = "legend") +
    scale_x_discrete(expand = expansion(mult = 0.02), na.translate = FALSE) +
    theme_cowplot() +
    theme(
      text = element_text(family = "Nunito Sans"),
      axis.text = element_text(family = "Roboto Slab"),
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank()
    ) + 
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
plot_cellexp_title <- "Gene Expression in CCLE Cell Lines."
plot_cellexp_legend <- paste0("Each point shows the ranked expression value for a given cell line. Black line indicates resampled mean expression value (", round(mean_virtual_expression, 2), "). Gray line indicates 3 standard deviations away from the resampled mean (", round(expression_upper, 2), ").")

## PROTEIN SIZE PLOT --------------------------------------------------------
# make_proteinsize <- function(protein_data = proteins, gene_symbol) {
# 
# proteins %>% 
#   ggplot() +
#   geom_histogram(aes(x = mass), bins = 50)
# 
#   return(plot_complete)
# }
# 
# #figure legend
# plot_proteinsize_title <- "Protein size"
# plot_proteinsize_legend <- paste0("Mass of ", gene_symbol, " compared to all protein masses. Density")

## CORRELATION PLOT --------------------------------------------------------
make_correlation <- function(table_data = achilles_cor_nest, gene_symbol, mean = mean_virtual_achilles, upper_limit = achilles_upper, lower_limit = achilles_lower) {
  ## use main color in case of 1 selected gene, otherwise use palette function
  if(length(gene_symbol) == 1) { 
    cols <- main_color 
  }else{
    cols <- color_pal_gene(length(gene_symbol))
  }

  plot_complete <- 
    table_data %>%
    dplyr::filter(fav_gene %in% gene_symbol) %>%
    tidyr::unnest(data) %>% 
    dplyr::group_by(fav_gene) %>% 
    arrange(desc(r2)) %>% 
    mutate(
      rank = 1:n(),
      med = median(r2, na.rm= TRUE)
    ) %>% 
    ungroup() %>% 
    ggplot() +
    ## sd square
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = upper_limit, ymax = lower_limit), fill = "gray90") +
    geom_hline(yintercept = 0, color = "gray80") +
    annotate("text", x = Inf, y = 0.005, label = "Gene Rank", color = "gray50", hjust = 1.15 ,vjust = 0) +
    ## dot plot
    geom_point(aes(x = rank, 
                   y = r2, 
                   text = paste0("Gene: ", gene), 
                   color = fct_reorder(fav_gene, med),
                   fill = fct_reorder(fav_gene, med)
                   ), 
               size = 1.1, stroke = .1, alpha = 0.4) +
    ## scales + legends
    scale_x_discrete(expand = expansion(mult = 0.02), na.translate = FALSE) +
    scale_color_manual(values = cols, guide = "legend") +
    scale_fill_manual(values = cols, guide = "legend") +
    guides(
      color = guide_legend(reverse = TRUE, override.aes = list(size = 5)),
      fill = guide_legend(reverse = TRUE, override.aes = list(size = 5))
    ) +
    ## theme changes
    theme_cowplot(font_size = 16) +
    theme(
      text = element_text(family = "Nunito Sans"),
      axis.text = element_text(family = "Roboto Slab"),
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(), 
      axis.line.x = element_blank()
    ) + 
    NULL
  
  ##change axis in case of more than 1 gene selected
  if(length(gene_symbol) == 1){
    plot_complete  <- plot_complete +
      ## titles
      labs(
        x = NULL, 
        y = paste0("Gene Correlations with ", gene_symbol),
        color = "Query Gene", 
        fill = "Query Gene"
      ) +
      guides(color = "none", fill = "none")
    plot_complete
  } else {
    plot_complete  <- plot_complete +
      ## titles
      labs(
        x = NULL, 
        y = "Gene Correlations",
        color = "Query Gene", 
        fill = "Query Gene" 
      )
  }
  return(plot_complete)
}

#figure legend
plot_genecorrelations_title <- "Gene Correlation Curve."
plot_genecorrelations_legend <- paste0("Each point shows the ranked gene correlation value ordered from high to low for each gene correlation set. Genes with correlation values outside the gray box and greather than ", round(achilles_upper, 2)," indicates the gene has a correlation value greater than ", sd_threshold, " standard deviations away from the mean. Conversely, genes with correlation values outside the gray box and less than ", round(achilles_lower, 2)," indicates the gene has a correlation value lower than ", sd_threshold, " standard deviations away from the mean.")

## EXPvDEP PLOT --------------------------------------------------------
make_expdep <- function(expression_data = expression, celldeps_data = achilles, expression_join = expression_names, gene_symbol) {
  ## use main color in case of 1 selected gene, otherwise use palette function
  if(length(gene_symbol) == 1) { 
    cols <- main_color 
  }else{
    cols <- color_pal_gene(length(gene_symbol))
  }
  exp_data <- expression_data %>% #plot setup
    dplyr::select(where(is.character), any_of(gene_symbol)) %>% 
    pivot_longer(cols = where(is.numeric), names_to = "gene_symbol", values_to = "exp")
  
  dep_data <- celldeps_data %>% #plot setup
    dplyr::select(where(is.character), any_of(gene_symbol)) %>% 
    pivot_longer(cols = where(is.numeric), names_to = "gene_symbol", values_to = "dep")  
  
  combined_data <- exp_data %>% 
    left_join(dep_data, by = c("X1", "gene_symbol")) %>% 
    filter(!is.na(dep)) %>% 
    dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>% 
    dplyr::mutate(med = median(dep, na.rm= TRUE)) %>% 
    dplyr::left_join(expression_join, by = "X1") %>% 
    dplyr::select(gene_symbol, exp, dep, med, cell_line, lineage)

  plot_complete <- 
    ggplot() +
    ## dot plot
    geom_rect(aes(xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf), fill = "gray90") +
    geom_point(data = combined_data, 
               aes(x = dep, 
                   y = exp, 
                   #text = paste0("Gene: ", gene_symbol, "\nCell Line: ", cell_line), 
                   color = fct_reorder(gene_symbol, med),
                   fill = fct_reorder(gene_symbol, med)),
               size = 2, stroke = .1, alpha = 0.4) +
    ## scales + legends
    scale_color_manual(values = cols, guide = "legend") +
    scale_fill_manual(values = cols, guide = "legend") +
    guides(
      color = guide_legend(reverse = TRUE, override.aes = list(size = 5)),
      fill = guide_legend(reverse = TRUE, override.aes = list(size = 5))
    ) +
    ## titles
    labs(
      x = "Dependency", 
      y = "Expression", 
      color = "Query Gene", 
      fill = "Query Gene"
    ) +
    ## theme changes
    theme_cowplot(font_size = 16) +
    theme(
      text = element_text(family = "Nunito Sans"),
      axis.text = element_text(family = "Roboto Slab")
    ) + 
    NULL
  
  if(length(gene_symbol) == 1){
    plot_complete  <- plot_complete +
      guides(color = "none", 
             fill = "none") +
       labs(x = paste0(gene_symbol, " Dependency"), 
            y = paste0(gene_symbol, " Expression"))
  } else {
    plot_complete
  }
  return(plot_complete)
}

#figure legend
plot_expdep_title <- "Gene Dependency versus Expression"
plot_expdep_legend <- paste0("Each point shows the dependency value compared to the expression value for gene within a given cell line. Gray area indicates dependency values that are between -1 and 1.")

