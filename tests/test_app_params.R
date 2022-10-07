library(here)

test_that("app_data_dir is set to 'data'", {
  source(here::here("code/app_params.R")) 
  expect_equal(app_data_dir, "data")
})

test_that("privateMode is set to TRUE", {
  source(here::here("code/app_params.R")) 
  expect_equal(privateMode, TRUE)
})
