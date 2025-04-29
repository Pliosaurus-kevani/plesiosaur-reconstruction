# Body reconstruction and size estimation of plesiosaurs
Code repository for Zhao 2025, run in R version 4.3.3.  
This repository contains commented code to enable readers to reproduce the analyses conducted in Zhao 2025.

The folder **"data"** contains the .csv files required in the R scripts.
The variable-rate model parameters computed by [BayesTraits 4.0.0](https://www.evolution.reading.ac.uk/BayesTraitsV4.0.0/BayesTraitsV4.0.0.html)
are also included here as .zip packages.

The folder **"phylogeny"** contains the tree files constructed using a Bayesian method (implemented in [RevBayes](https://revbayes.github.io/)),
which will be described in detail in another recent publication of mine.
It contains the maximum clade credibility tree and 100 trees randomly selected from the post burn-in Bayesian samples.

The folder **"scripts"** contains the R scripts required to reproduce the analyses and figures.

**restoration_PGLS & OLS.R**
>This script fits regression models using ordinary least squares (OLS), phylogenetic generalized least squares (PGLS)
>or nonlinear functions to reconstruct the missing puzzles in plesiosaur reconstruction.

**find proxy_OLS & PGLS.R**
>This script fits regression models using ordinary least squares (OLS) and phylogenetic generalized least squares (PGLS)
>to evaluate the performance of various skeletal structures as size proxies.

**generate files for BayesTraits.R**
>This script generates the command files and data files required by [BayesTraits 4.0.0](https://www.evolution.reading.ac.uk/BayesTraitsV4.0.0/BayesTraitsV4.0.0.html).

**convergence and log Bayes Factors.R**
>This script checks the convergence of parameters estimated using [BayesTraits 4.0.0](https://www.evolution.reading.ac.uk/BayesTraitsV4.0.0/BayesTraitsV4.0.0.html), and compute the log Bayes Factors from marginal likelihoods.

**BayesTraits_plots.R**
>This script summarizes the variable-rate model and maps the rates onto the branches of the maximum clade credibility tree.
