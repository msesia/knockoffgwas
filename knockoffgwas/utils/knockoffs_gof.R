#!/usr/bin/env Rscript

# Load packages
suppressMessages(library(tidyverse))
suppressMessages(library(latex2exp))

# Default arguments (for debugging)
chr.name <- 21
res.name <- 1
stats.basename  <- "../../tmp/knockoffs_gof/knockoffs_chr21_res1"
groups.file <- "../../tmp/knockoffs/knockoffs_chr21_res1_grp.txt"
out.basename  <- "../../tmp/knockoffs_gof_plots/knockoffs_chr21_res1"

# Input arguments
args <- commandArgs(trailingOnly=TRUE)
chr.name <- as.character(args[1])
res.name <- as.character(args[2])
stats.basename  <- as.character(args[3])
groups.file <- as.character(args[4])
out.basename  <- as.character(args[5])


plot.knockoff.diagnostics <- function(chr.name, res.name, stats.basename, groups.file) {

    ## Load variable grouping
    Variants <- read_delim(groups.file, delim=" ", progress = FALSE, col_types = cols()) %>%
        mutate(Group = as.integer(Group))

    ## Load variant frequency table
    frq.file <- sprintf("%s_self.frq", stats.basename)
    Frq <- read_table(frq.file, col_types=cols())
    ## Compute diagnostics
    Frq <- Frq %>%
        separate(SNP, ".k", into=c("SNP", "Knockoff")) %>%
        mutate(Knockoff = ifelse(is.na(Knockoff), FALSE, TRUE))
    Diagnostics <- Frq %>%
        pivot_wider(names_from="Knockoff", values_from=c("MAF")) %>%
        mutate(x=`FALSE`, xk=`TRUE`) %>%
        select(CHR, SNP, x, xk)
    ## Plot frequency diagnostics
    p.frq <- Diagnostics %>%
        ggplot(aes(x=x, y=xk)) +
        geom_point(alpha=0.2) +
        geom_abline(intercept = 0, slope = 1, color="red", linetype=2) +
        xlab(TeX("MAF ($X$)")) + ylab(TeX("MAF ($\\tilde{X}$)")) +
        xlim(0,1) +
        ylim(0,1) +
        theme_bw()
    
    ## Load LD table
    ld.file <- sprintf("%s.ld", stats.basename)
    LD <- suppressWarnings(read_table2(ld.file, col_types=cols())) %>%
        mutate(CHR=CHR_A) %>% select(-CHR_A, -CHR_B)

    ## Add grouping information
    LD <- LD %>%
        mutate(SNP=SNP_A) %>%
        mutate(SNP = gsub(".k", "", SNP)) %>%
        inner_join(Variants %>% select(SNP, Group), by = c("SNP")) %>%
        mutate(Group_A=Group) %>% select(-SNP, -Group) %>%
        mutate(SNP=SNP_B) %>%
        mutate(SNP = gsub(".k", "", SNP)) %>%
        inner_join(Variants %>% select(SNP, Group), by = c("SNP")) %>%
        mutate(Group_B=Group) %>% select(-SNP, -Group)

    ## Add knockoff key information
    LD <- LD %>%
        separate(SNP_A, ".k", into=c("SNP_A", "Knockoff_A")) %>%
        mutate(Knockoff_A = ifelse(is.na(Knockoff_A), FALSE, TRUE)) %>%
        separate(SNP_B, ".k", into=c("SNP_B", "Knockoff_B")) %>%
        mutate(Knockoff_B = ifelse(is.na(Knockoff_B), FALSE, TRUE))

    # Create correlation tables between different groups
    group.range <- seq(0,10)
    LD.XX <- LD %>%
        filter(abs(Group_B-Group_A) %in% group.range, Knockoff_A==FALSE, Knockoff_B==FALSE) %>%
        mutate(R.XX=R2) %>%
        mutate(SNP_A=str_replace(SNP_A,".A",""), SNP_A=str_replace(SNP_A,".B","")) %>%
        mutate(SNP_B=str_replace(SNP_B,".A",""), SNP_B=str_replace(SNP_B,".B","")) %>%
        select(Group_A, Group_B, SNP_A, SNP_B, R.XX) %>%
        distinct(Group_A, Group_B, SNP_A, SNP_B, R.XX)
    LD.XkXk <- LD %>%
        filter(abs(Group_B-Group_A) %in% group.range, Knockoff_A==TRUE, Knockoff_B==TRUE) %>%
        mutate(R.XkXk=R2) %>%
        mutate(SNP_A=str_replace(SNP_A,".A",""), SNP_A=str_replace(SNP_A,".B","")) %>%
        mutate(SNP_B=str_replace(SNP_B,".A",""), SNP_B=str_replace(SNP_B,".B","")) %>%
        select(Group_A, Group_B, SNP_A, SNP_B, R.XkXk) %>%
        distinct(Group_A, Group_B, SNP_A, SNP_B, R.XkXk)
    LD.XXk <- LD %>%
        filter((Group_B-Group_A) %in% seq(1,10), Knockoff_A*Knockoff_B==FALSE) %>%
        mutate(R.XXk=R2) %>%
        mutate(SNP_A=str_replace(SNP_A,".A",""), SNP_A=str_replace(SNP_A,".B","")) %>%
        mutate(SNP_B=str_replace(SNP_B,".A",""), SNP_B=str_replace(SNP_B,".B","")) %>%
        select(Group_A, Group_B, SNP_A, SNP_B, R.XXk) %>%
        distinct(Group_A, Group_B, SNP_A, SNP_B, R.XXk)

    ## Plot originality
    LD.cross <- inner_join(LD.XX, LD.XkXk, by = c("Group_A", "Group_B", "SNP_A", "SNP_B"))
    p.orig <- LD.cross %>%
        mutate(Distance = as.factor(abs(Group_A-Group_B))) %>%
        ggplot(aes(x=abs(R.XX), y=abs(R.XkXk))) +
        geom_abline(color="red") +
        geom_point(alpha=0.1) +
        xlim(0,1) + ylim(0,1) +
        xlab(TeX("|corr($X_{j},X_{k}$)|")) + ylab(TeX("|corr($\\tilde{X}_{j},\\tilde{X}_{k}$)|")) +
        theme_bw()
    
    ## Plot exchangeability
    options(repr.plot.width=4, repr.plot.height=3)
    LD.cross <- inner_join(LD.XX, LD.XXk, by = c("Group_A", "Group_B", "SNP_A", "SNP_B")) %>%
        filter(Group_A!=Group_B)
    p.exch <- LD.cross %>%
        mutate(Distance = as.factor(abs(Group_A-Group_B))) %>%
        ggplot(aes(x=abs(R.XX), y=abs(R.XXk))) +
        geom_abline(color="red") +
        geom_point(alpha=0.1) +
        xlim(0,1) + ylim(0,1) +
        xlab(TeX("|corr($X_{j},X_{k}$)|")) + ylab(TeX("|corr($X_{j},\\tilde{X}_{k}$)|")) +
        theme_bw()
    
    ## Plot histogram of self-correlations
    p.self <- LD %>% filter(BP_A==BP_B) %>%
        ggplot(aes(x=R2)) +
        geom_histogram(bins=30) +
        xlab(TeX("|corr($X_{j},\\tilde{X}_{j}$|)")) +
        theme_bw()

    ## Combine plots
    plot.title <- sprintf("Knockoff GOF for chromosome %s, resolution %s", chr.name, res.name)
    p.combined <- gridExtra::grid.arrange(p.frq, p.orig, p.exch, p.self, nrow = 2,
                                          top = grid::textGrob(plot.title,gp=grid::gpar(fontsize=15,font=1)))
    return(p.combined)
}

## Make plot
pp <- plot.knockoff.diagnostics(chr.name, res.name, stats.basename, groups.file)

## Save plot
out.file <- sprintf("%s.png", out.basename)
ggsave(out.file, plot=pp)
cat(sprintf("GOF plots saved on %s\n", out.file))
