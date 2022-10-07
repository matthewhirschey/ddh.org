#drug repurposing

repurposing_url <- "https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_drugs_20200324.txt"

repurposing <- 
  read_tsv(repurposing_url, skip = 9)

repurposing %>% filter(pert_iname %in% prism_names$name)
repurposing %>% filter(!pert_iname %in% prism_names$name)

prism_meta %>% filter(!name %in% repurposing$pert_iname)
