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
packages.cran <- c("shiny", "tidyverse", "gridExtra", "dqshiny", "tibbletime", "scales",
                   "cowplot", "latex2exp", "ggrepel", "egg", "grid")
check.packages.cran(packages.cran)
