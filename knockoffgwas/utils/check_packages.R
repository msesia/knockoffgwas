#!/usr/bin/env Rscript
# 

########################
## Packages from CRAN ##
########################

check.packages.cran <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    sapply(pkg, require, character.only = TRUE)
}
packages.cran <- c("tidyverse", "fastcluster", "latex2exp", "grid", "gridExtra")
check.packages.cran(packages.cran)

################################
## Packages from Bioconductor ##
################################

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
check.packages.bio <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        BiocManager::install(new.pkg)
    sapply(pkg, require, character.only = TRUE)
}
packages.bio <- c()
check.packages.bio(packages.bio)

#####################
## Github packages ##
#####################
check.packages.dev <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        devtools::install_github(sprintf("privefl/%s", new.pkg))
    sapply(pkg, require, character.only = TRUE)
}
packages.dev <- c("bigstatsr", "bigsnpr")
check.packages.dev(packages.dev)
