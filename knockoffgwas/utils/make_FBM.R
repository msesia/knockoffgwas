#!/usr/bin/env Rscript

# Install package bigsnpr
# devtools::install_github("privefl/bigsnpr")

# Documentation here:
# https://privefl.github.io/bigsnpr/reference/index.html

# Default arguments (for debugging)
geno.basename <- "../../tmp/knockoffs_full/example_res7"
out.basename  <- "../../tmp/knockoffs_full/example_res7"

# Input arguments
args <- commandArgs(trailingOnly=TRUE)
geno.basename <- as.character(args[1])
out.basename  <- as.character(args[2])

cat("Converting augmented genotypes into FBM...\n")

# Load packages
suppressMessages(library(devtools))
suppressMessages(library(bigsnpr))

# Delete temporary file if it exists
tmp.file <- sprintf("%s.bk", out.basename)
if(file.exists(tmp.file)) {
    file.remove(tmp.file)
}
# Load and convert genotypes
bed.file <- sprintf("%s.bed", geno.basename)
cat("Reading genotype file and creating FBM object... ")
dat <- snp_readBed(bed.file, backingfile = out.basename)
cat("done.\n")

cat(sprintf("Written %s.rds\n", out.basename))
