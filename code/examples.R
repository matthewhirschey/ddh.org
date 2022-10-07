source(here::here("code", "app_params.R"), local = TRUE)

#need to dynamically switch depending on privateMode
tld <- dplyr::if_else(privateMode == TRUE, "com", "org")
url <- glue::glue('http://www.datadrivenhypothesis.{tld}/')

#need full urls so works in methods, instead of reference url
examples <- glue::glue('<h4>Search for genes</h5>
              <ul>
                <li>A single gene, such as <a href="{url}?show=gene&query=TP53">TP53</a> or <a href="{url}?show=gene&query=BRCA1">BRCA1</a></li>
                <li>A pathway name, such as <a href="{url}?show=search&query=cholesterol">cholesterol</a>, which will lead you to <a href="{url}?show=pathway&query=0006695">Cholesterol Biosynthetic Process</a></li>
                <li>The Gene Ontology biological process identifier, such as <a href="{url}?show=search&query=1901989">1901989</a>, which will find <a href="{url}?show=pathway&query=1901989">Pathway: Positive Regulation Of Cell Cycle Phase Transition (GO:1901989)</a></li>
                <li>A custom list of genes (separated by commas), such as <a href="{url}?show=search&query=BRCA1,%20BRCA2">BRCA1, BRCA2</a>, which will search <a href="{url}?show=gene_list&query=BRCA1,BRCA2">a custom gene list</a></li>
              </ul>') 
