# Body reconstruction and size estimation of plesiosaurs
Code repository for Zhao 2025, run in R version 4.3.3.  
This repository contains commented code to enable readers to reproduce the analyses conducted in Zhao, 2025.

The folder **"data"** contains the files required in the R scripts. The folder **"raw"** consists of the .xlsx files
containing the raw data and a PDF file containing the references. The detailed references of the data sources were submitted as part of the supplementary material of the article. Note that all species names in the dataset have been standardized to match those in the phylogeny, if they are present in the tree. Any discrepancies with recent taxonomic revisions do not represent the author's opinion on taxonomy.

The file **"consensus.nex"** is the strict consensus tree from [Sachs & Madzia, 2025](https://doi.org/10.7717/peerj.19665) (using implied weighting, K = 12).

The folder **"scripts"** contains the R scripts required to reproduce the analyses and figures.

**time-calibration.R**
>This script performs time‑calibration of the strict consensus tree under the minimum branch length (mbl) method. Polytomies are randomly resolved.

**restoration_PGLS & OLS.R**
>This script fits regression models using ordinary least squares (OLS), phylogenetic generalized least squares (PGLS)
>or a nonlinear function (log-logistic) to reconstruct the missing puzzles in plesiosaur reconstruction.
>It also generates the figure 8 in the manuscript.

**find proxy_OLS & PGLS.R**
>This script fits regression models using ordinary least squares (OLS) and phylogenetic generalized least squares (PGLS)
>to evaluate the performance of various skeletal structures as mass proxies.
>It also generates the figure 6 in the manuscript.
