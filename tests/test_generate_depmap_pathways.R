library(here)
library(tidyverse)

context("generate_depmap_pathways")

source(here::here("code/generate_depmap_pathways.R"), local=TRUE)

test_that("enrichr_loop with empty gene_list can be arranged by Adjusted.P.value", {
  # test for error: "arrange() failed at implicit mutate() step.
  #   Could not create a temporary column for `Adjusted.P.value`."
  result = enrichr_loop(c(), focused_lib) %>%
    arrange(Adjusted.P.value)
  expect_equal(nrow(result), 0)
})

test_that("cast_enrichr_data casts Term, Overlap and Genes", {
  result = cast_enrichr_data(tibble(
    Term=1:3,
    Overlap=2:4,
    Genes=3:5,
  ))
  expect_equal(result$Term, c("1","2","3"))
  expect_equal(result$Overlap, c("2","3","4"))
  expect_equal(result$Genes, c("3","4","5"))
})
