# CMBE-Bacteroides-ECM


## Overview

This repository contains scripts, files, and data used in the genomic analysis and annotation of Bacteroides genomes to identify carbohydrate-active enzymes and proteases specific to extracellular matrix components. The authors of this manuscript are 



### Annotations

| File Name           | Input                | Output                | Type                 | Description                        |
|---------------------|----------------------|-----------------------|----------------------|----------------------------------|
| scripts_Prokka_eggNOG_run-dbCAN_MEROPS     | genome.fasta         | prokka_output<br> eggNOG_results<br> dbcan_results.csv<br> merops_results.csv     | Shell script (.sh)<br> Ran on HPC    | Runs Prokka genome annotation    |



### ECM_database

| File Name           | Input                | Output                | Type                 | Description                        |
|---------------------|----------------------|-----------------------|----------------------|----------------------------------|
| build_ecm_database.R   | None                 | UniProt_PLUS_BRENDA_CollapsedByEC.csv<br> (supplementary_table1.csv) | R script (.R)                 | ECM component sequences database  |
| map_ec_cazyme_protease_annotations.R   | UniProt_PLUS_BRENDA_CollapsedByEC.csv<br> (supplementary_table1.csv)        | filtered_ecm.csv<br> supplementary_table2.csv      | R script (.R)    | Script to parse ECM sequence data |



### Statistical_analysis

| File Name           | Input                | Output                | Type                 | Description                        |
|---------------------|----------------------|-----------------------|----------------------|----------------------------------|
| mantel_test.R       | distance_matrix.csv  | mantel_results.txt    | R script (.R)        | Performs Mantel test              |
| spearman_corr.R    | abundance_table.csv  | spearman_results.csv  | R script (.R)   | Calculates Spearman's correlation|



### Figures

| File Name           | Input                | Output                | Type                 | Description                        |
|---------------------|----------------------|-----------------------|----------------------|----------------------------------|
| plot_dotplot.R      | genome_comparison.csv| dotplot.png           | R script (.R)        | Creates dotplot of protease/CAZymes families<br> (Figure 3)|
| heatmap_plot.R     | enzyme_abundance.csv | heatmap.png           | R script (.R)  | Generates heatmap of enzyme data<br> (Figure 4) |



### Supplementary_files


| File Name           | Type                 | Description                        |
|---------------------|----------------------|----------------------------------|
| supplementary_table1.csv |   file (.csv)    | Full ECM database including all BRENDA results from individual ECM terms and the UniProtKB results            |
| supplementary_table2.csv|  file (.csv)  | *Bacteroides* annotations with mappings to ECM substrates    |
