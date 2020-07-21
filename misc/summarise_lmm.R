#!/usr/bin/env Rscript

# Default parameters
clump.file <- "../data/lmm/example_lmm_clumped.tab"
out.file   <- "../data/lmm/example_lmm_clumped.txt"

# Load libraries
suppressMessages(library(tidyverse))

# Load clumped results
Clumped.raw <- read.table(clump.file, header=TRUE) %>% as_tibble() %>% select(CHR, SNP, BP, SP2)

# Load list of variants
chr.list <- unique(Clumped.raw$CHR)
Variants <- lapply(chr.list, function(chr) {
    cat(sprintf("Loading list of variants for chromosome %d...\n", chr))
    bim.file <- sprintf("../data/genotypes/example_chr%d.bim", chr)
    read_tsv(bim.file, col_names=c("CHR", "SNP", "X1", "BP", "A0", "A1"), col_types=cols())
})
Variants <- do.call("rbind", Variants)

# Parse clumped results
if(nrow(Clumped.raw)>0) {
    # Extract lead SNPs
    Clumped.raw.lead <- Clumped.raw %>% mutate(SP2=SNP) %>% select(CHR, SNP, BP, SP2)
    # Extract secondary SNPs
    Clumped.raw.secondary <- Clumped.raw %>%
        mutate(SP2 = gsub("\\(1\\)", "", SP2)) %>%
        separate_rows(SP2, sep=",") %>%
        select(CHR, SNP, BP, SP2) %>%
        filter(SP2!="NONE")
    # Cross-reference with complete list of variants
    Clumped <- rbind(Clumped.raw.lead, Clumped.raw.secondary) %>%
        as_tibble() %>% arrange(CHR, BP, SP2) %>%
        mutate(SNP.lead=as.character(SNP), BP.lead=BP, SNP=as.character(SP2)) %>%
        select(CHR, SNP.lead, BP.lead, SNP) %>%
        left_join(select(Variants, CHR, SNP, BP), by=c('CHR', 'SNP')) %>%
        select(CHR, SNP.lead, BP.lead, SNP, BP)
} else {
    Clumped <- tibble()
}

# Load the LMM p-values and add them to the list of clumped results
lmm.file <- sprintf("../data/lmm/example_lmm.txt")
LMM <- read_tsv(lmm.file, col_types=cols())
if(! "P" %in% colnames(LMM)) {
    if("P_BOLT_LMM" %in% colnames(LMM)) {
        LMM <- LMM %>% mutate(P=P_BOLT_LMM)
    } else {
        LMM <- LMM %>% mutate(P=P_BOLT_LMM_INF)
    }
}
Clumped <- Clumped %>% inner_join(LMM %>% select(CHR,SNP,BP,P), by = c("CHR", "SNP", "BP"))

# Summarise results by clump
Discoveries <- Clumped %>% group_by(CHR, SNP.lead, BP.lead) %>%
    summarise(P=min(P), BP.min=min(BP), BP.max=max(BP), BP.width=BP.max-BP.min, Size=n())

# Save results
Clumped %>% write_delim(out.file, delim=" ")
cat(sprintf("List of discovered clumps written on: %s\n", out.file))
