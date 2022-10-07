library(here)
library(tidyverse)
source(here::here("code/fun_tables.R"))

test_that("sort_dedup_and_limit removes lower ranked duplicate items", {
  df <- data.frame(
    key=c('A','B','A', 'B'),
    rank=c(0.5, 0.8, 0.9, 0.6))
  expected_result <- data.frame(
    key=c('A', 'B'),
    rank=c(0.9, 0.8))
  result <- sort_dedup_and_limit(df, limit_rows=4)
  expect_equal(result, expected_result)
})

test_that("search_gene_data filters out duplicate items for ADSS", {
  gene_summary <- tibble(
    approved_symbol=c('ADSS1','ADSS2'),
    aka=c('ADSSL1','ADSS'),
    approved_name=c('Adenylo1','Adenylo2')
  )
  pathways <- tibble(
    pathway=character(),
    go=character()
  )
  query_str <- 'ADSS'
  result <- search_gene_data(gene_summary, pathways, query_str, limit_rows=4)
  expect_equal(result %>% pull(rank), c(1.0,0.8))
  expect_equal(result %>% pull(key), c('ADSS2','ADSS1'))
})
