#!/usr/bin/env Rscript

# Input arguments
resolution <- "Radj1"

args <- commandArgs(trailingOnly=TRUE)
resolution <- as.character(args[1])

suppressMessages(library(tidyverse))

scratch <- "/scratch/groups/candes/ukbiobank_tmp"

# Load list of variants
Variants <- tibble()
chr.list <- seq(1,22)
for(chr in chr.list) {
    cat(sprintf("Loading list of variants on chromosome %d... ", chr))
    key.file <- sprintf("%s/knockoffs/%s_K50/ukb_gen_chr%d.key", scratch, resolution, chr)
    Variants.chr <- read_delim(key.file, delim=" ", col_types=cols())
    Variants.chr <- Variants.chr %>% mutate(CHR=Chr) %>% select(CHR, Variant, Position, Group, Knockoff)
    colnames(Variants.chr) <- c("CHR", "SNP", "BP", "Group", "Knockoff")
    Variants <- rbind(Variants, Variants.chr)
    cat("done.\n")
}

# Extract list of original variables
original.file <- sprintf("%s/augmented_data_big/ukb_gen_%s_original.txt", scratch, resolution)
Variants %>%
    filter(Knockoff==FALSE) %>%
    select(SNP) %>%
    write_delim(original.file, delim=" ", col_names=F)
cat(sprintf("List of original variants written on:\n %s\n", original.file))

knockoff.file <- sprintf("%s/augmented_data_big/ukb_gen_%s_knockoff.txt", scratch, resolution)
Variants %>%
    filter(Knockoff==TRUE) %>%
    select(SNP) %>%
    write_delim(knockoff.file, delim=" ", col_names=F)
cat(sprintf("List of knockoff variants written on:\n %s\n", knockoff.file))
