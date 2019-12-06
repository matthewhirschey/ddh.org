#### Data file generation
To generate the data files, run:  
1. code/create_gene_summary.R  
2. code/generate_depmap_data.R  
3. code/generate_depmap_stats.R  
4. code/generate_depmap_tables.R  
5. code/generate_depmap_pathways.R  
  
The files generated in steps 1-3 are required for steps 4 and 5. Step 4 takes about 60' to run locally. Step 5 requires some parallization, and you'll see objects dec(ile)1-10 that could be run in parallel. The code for step 5 has `gene_group <- sample` so it can be tested. 
