source(here::here("code", "fun_structures.R"), local = TRUE)
source(here::here("code", "fun_text.R"), local = TRUE)
source(here::here("code", "shiny_structures.R"), local = TRUE)

testShinyModule("compoundStructureServer", args = compound_data_args, {
        c(output$compound_structure)
})

testShinyModule("cellImageServer", args = cell_data_args, {
      c(output$cell_image)
})

