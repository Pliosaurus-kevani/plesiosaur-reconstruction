# Body reconstruction and size estimation of plesiosaurs
Code repository for Zhao 2025, run in R version 4.3.3.  
This repository contains commented code to enable readers to reproduce the analyses conducted in Zhao 2025.

The folder **"data"** contains the files required in the R scripts. The folder **"raw"** consists of the .xlsx files
containing the raw data. The detailed references of the data sources were submitted as part of the supplementary material of the article. Note that all species names in the dataset have been standardized to match those in the phylogeny, if they are present in the tree. Any discrepancies with recent taxonomic revisions do not represent the author's opinion on taxonomy.

The folder **"scripts"** contains the R scripts required to reproduce the analyses and figures.

**restoration_PGLS & OLS.R**
>This script fits regression models using ordinary least squares (OLS), phylogenetic generalized least squares (PGLS)
>or nonlinear functions to reconstruct the missing puzzles in plesiosaur reconstruction.

**find proxy_OLS & PGLS.R**
>This script fits regression models using ordinary least squares (OLS) and phylogenetic generalized least squares (PGLS)
>to evaluate the performance of various skeletal structures as size proxies.
