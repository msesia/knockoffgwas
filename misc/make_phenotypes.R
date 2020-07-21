library(tidyverse)
library(SNPknockG)
library(snpStats)

out.dir <- "../data"

amplitude <- 8
chr.list <- seq(21,22)
set.seed(2019)

# Load list of variants
Variants <- lapply(chr.list, function(chr) {
    bim   <- sprintf("../data/genotypes/example_chr%d.bim", chr)
    Variants.chr <- read_delim(bim, delim="\t",
                               col_names=c("CHR", "SNP", "X", "BP", "A1", "A2"), col_types=cols())
    return(Variants.chr)
})
Variants <- do.call("rbind", Variants)

# Choose the causal variants
choose_causal <- function(n.loci, n.signals) {

    # Check that the signals can be split exactly into the desired number of loci
    locus.size <- n.signals/n.loci
    if(locus.size!=round(locus.size)) {
        cat(sprintf("Error: %d is not divisible by %d.\n", n.signals, n.loci))
        return(NULL)
    }

    # Allocate causal SNPs to chromosomes
    Variants.num <- Variants %>% group_by(CHR) %>% summarise(N=n())
    Variants.num <- Variants.num %>% mutate(Share = N/sum(Variants.num$N), Loci = round(Share*n.loci))
    # Add signals if rounding errors occurred
    if(sum(Variants.num$Loci) != n.loci) {
        Variants.num$Loci[1] = Variants.num$Loci[1]+(n.loci-sum(Variants.num$Loci))
    }
    n.loci.chr <- Variants.num$Loci
    names(n.loci.chr) <- chr.list
    cat(sprintf("Divided %d causal loci into %d chromosomes.\n", sum(n.loci.chr),nrow(Variants.num)))

    # Distance between causal SNPs (in kB)
    causal.dist <- round(sort(runif(locus.size, 0, 100)),2)

    # Pick causal SNPs
    Variants.causal <- tibble()
    for(chr in chr.list) {
        chr.idx <- as.character(chr)
        cat(sprintf("Assembling %d clumps of size %d for chromosome %d...\n",
                    n.loci.chr[chr.idx], locus.size, chr))
        Variants.chr <- Variants %>% filter(CHR==chr)
        causal.idx <- round(seq(round(0.02*nrow(Variants.chr)),
                                round(0.98*nrow(Variants.chr)), length.out=n.loci.chr[chr.idx]))
        bp.first <- Variants.chr$BP[causal.idx]
        bp.next <- as.numeric(sapply(causal.dist[-1], function(dkb) {bp.first + 1000 * dkb}))
        bp.next <- matrix(bp.next, ncol=length(causal.dist)-1)
        locus.idx <- seq(length(causal.idx))
        for(j in 1:length(bp.first)) {
            if(length(causal.dist)>1) {
                for(k in 1:ncol(bp.next)) {
                    candidates <- order(abs(Variants.chr$BP-bp.next[j,k]))
                    candidates <- candidates[which(candidates>causal.idx[j])]
                    candidates <- setdiff(candidates, causal.idx)
                    causal.idx <- c(causal.idx, candidates[1])
                    locus.idx <- c(locus.idx, j)
                }
            }
        }
        Variants.causal.chr <- Variants.chr[causal.idx,]
        Variants.causal.chr$Locus <- locus.idx
        Variants.causal.chr <- Variants.causal.chr %>% arrange(BP) %>% distinct()
        Variants.causal <- rbind(Variants.causal, Variants.causal.chr)
    }

    # Assign effect sizes
    Variants.causal <- Variants.causal %>%
        mutate(Locus=rep(seq(1,n.loci),each=locus.size)) %>%
        mutate(Sign=rep(2*rbinom(nrow(Variants.causal)/locus.size,1,0.5)-1,each=locus.size)) %>%
        mutate(Scale=round(runif(nrow(Variants.causal),0.1,1.9),3))

    return(Variants.causal)
}

Variants.causal <- choose_causal(20, 40)

linear.model <- function(X,beta) {
    y <- X %*% beta + rnorm(nrow(X))
    return(y)
}

# Load the causal genotypes
cat(sprintf("Loading genotypes for causal loci... "))
G <- lapply(chr.list, function(chr) {
    bed   <- sprintf("../data/genotypes/example_chr%d.bed", chr)
    bim   <- sprintf("../data/genotypes/example_chr%d.bim", chr)
    fam   <- sprintf("../data/genotypes/example_chr%d.fam", chr)
    dat <- read.plink(bed, bim, fam)
    G.chr <- dat$genotypes
    SNP <- dat$map$snp.name
    col.idx <- which(SNP %in% Variants.causal$SNP)
    X.chr <- as(G.chr[,col.idx], "numeric")
    return(X.chr)
})
G <- do.call("cbind", G)
cat("done.\n")

# Load sample information
fam <- sprintf("../data/genotypes/example_chr%d.fam", chr.list[1])
Subjects <- read_delim(fam, delim="\t", col_names=c("FID", "IID", "X0", "X1", "sex", "X2"),
                         col_types=cols()) %>%
    select(FID, IID, sex)

# Relevant covariates
Z <- 2-Subjects$sex

# Compute phenotype mean vector
cat(sprintf("Generating phenotypes... "))
beta <- amplitude * Variants.causal$Sign * Variants.causal$Scale / sqrt(nrow(G))
y.mean <- 10*scale(Z)/sqrt(nrow(G)) + scale(G) %*% beta
y <- as.numeric(y.mean + rnorm(nrow(G)))
cat("done.\n")

# Store phenotypes
Phenotypes <- Subjects %>% mutate(y=round(y,3))
pheno.file <- sprintf("%s/phenotypes/phenotypes.tab", out.dir)
Phenotypes %>% write_delim(pheno.file, delim="\t")
cat(sprintf("Phenotypes written on: %s\n", pheno.file))

# Store causal variants
causal.file <- sprintf("%s/phenotypes/causal.txt", out.dir)
Variants.causal %>% write_delim(causal.file, delim=" ")
cat(sprintf("List of causal variants written on: %s\n", causal.file))
