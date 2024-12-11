# MultiRFM
High-dimensional multi-study robust factor model

=========================================================================

To robustly extract meaningful features from data derived from multiple heterogeneous sources, we introduce a high-dimensional multi-study robust factor model, called MultiRFM, which learns latent features and accounts for the heterogeneity among sources.


# Installation
"MultiRFM" depends on the 'Rcpp' and 'RcppArmadillo' package, which requires appropriate setup of computer. For the users that have set up system properly for compiling C++ files, the following installation command will work.
```{Rmd}
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("feiyoung/MultiRFM")


```


## Simulated codes
For the codes in simulation study, check the `simuCodes` directory of the repo.


## News

MultiRFM version 1.1 released! (2024-12-11) 


