expression_long %>% 
  select(-protein_expression) %>% 
  filter(!is.na(gene_expression),
         #X1 == "ACH-001113",
         gene_expression > 0) %>% 
  group_by(X1) %>% 
  summarize(n = n(), min = min(gene_expression), max = max(gene_expression)) %>% 
  arrange(n) %>% 
  slice(1) %>% 
  left_join(expression_names, by = "X1")

#ACH-001097 has the fewest

#18567 is most
12876/18567

(18567-12876)/18567

