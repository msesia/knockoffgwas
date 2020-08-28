#!/usr/bin/env Rscript

# Install packages bigsnpr and bigstatsr
# devtools::install_github("privefl/bigstatsr")
# devtools::install_github("privefl/bigsnpr")

# Documentation here:
# https://privefl.github.io/bigsnpr/reference/index.html
# https://privefl.github.io/bigstatsr/reference/big_spLinReg.html

# Load packages
suppressMessages(library(tidyverse))
suppressMessages(library(bigsnpr))

# Default arguments (for debugging)
basename       <- "../../tmp/knockoffs_full/example_res0"
part.file     <- "../../tmp/partitions/example_chr21.txt"
pheno.file     <- "../../data/phenotypes/phenotypes.tab"
pheno.name     <- "y"
out.basename   <- "../../tmp/example_res0"

# Input arguments
args <- commandArgs(trailingOnly=TRUE)
basename       <- as.character(args[1])
pheno.file     <- as.character(args[2])
pheno.name     <- as.character(args[3])
out.basename   <- as.character(args[4])

# Other parameters
ncores <- 1
dfmax  <- 10000

# Set random seed for cross-validated lasso statistics
set.seed(2020)

####################
## Load genotypes ##
####################
fbm.file <- sprintf("%s.rds", basename)

# Attach the "bigSNP" object in R session
cat("Attaching bigSNP object... ")
obj.bigSNP <- snp_attach(fbm.file)
cat("done.\n")

# Extract list of variants
map <- obj.bigSNP$map %>% as_tibble()
colnames(map) <- c("CHR", "SNP", "gd", "BP", "a1", "a2")
map <- map %>% select(CHR, SNP, BP)

# Extract list of subjects
Subjects <- obj.bigSNP$fam %>% as_tibble()
colnames(Subjects) <- c("FID", "IID", "X1", "X2", "sex", "X3")
Subjects <- Subjects %>% select(FID, IID) %>%
    mutate(FID=as.character(FID), IID=as.character(IID))

# Get aliases for useful slots
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
SNP <- obj.bigSNP$map$marker.ID

#############################
## Load variant partitions ##
#############################

## Load list of variants and partitions
bim.file <- sprintf("%s.bim", basename)
Variants <- read_delim(bim.file, delim="\t", col_names=c("CHR", "SNP", "X0", "BP", "X1", "X2"),
                       col_types=cols(.default=col_integer(), CHR=col_character(), SNP=col_character())) %>%
    separate(SNP, c("SNP", "Knockoff"), fill="right") %>%
    mutate(Knockoff = ifelse(is.na(Knockoff), FALSE, TRUE)) %>%
    select(CHR, SNP, BP, Knockoff)
grp.file <- sprintf("%s_grp.txt", basename)
df <- read_delim(grp.file, delim=" ", col_names=c("SNP", "Group"),
                 col_types=cols(.default=col_integer(), SNP=col_character()))
Variants <- Variants %>% left_join(df, by = "SNP")

# Compute scaling factor for the genotypes
cat("Computing scaling factors for all variants... ")
scaler <- big_scale()
G.scale <- scaler(G)
scaling.factors <- G.scale$scale
cat("done.\n")

#####################
## Load phenotypes ##
#####################

cat("Reading phenotype file... ")
# Load phenotype table
Phenotypes <- read_delim(pheno.file, delim="\t", col_types=cols()) %>%
    mutate(FID=as.character(FID), IID=as.character(IID))
# Make sure that the rows of the genotypes match the rows of the phenotypes
Phenotypes <- Phenotypes %>% right_join(Subjects, by=c("FID", "IID"))
cat("done.\n")

###################
## Fit the lasso ##
###################

# Extract response variable
ind.train <- which(!is.na(Phenotypes[[pheno.name]]))
cat(sprintf("%d individuals out of %d have missing phenotype.\n",
            nrow(Subjects)-length(ind.train), nrow(Subjects)))
y <- Phenotypes[[pheno.name]][ind.train]

# Extract covariates (sex)
Covariates <- Phenotypes %>% select(sex)
covar.train <- as.matrix(Covariates)[ind.train,,drop=F]

# Find the class of response (numeric or binary factor)
y.unique <- unique(y[!is.na(y)])
if(length(y.unique)==2) {
    phenotype.class <- "binary"
    y <- factor(y, levels=c(1,2), labels=c(0,1))
    y <- as.numeric(levels(y))[y]
} else {
    phenotype.class <- "continuous"
}

# Fit the lasso
if(phenotype.class=="binary") {
    cat(sprintf("Fitting sparse logistic regression with %d observations, %d variants and %d covariates... ",
                length(y), ncol(G), ncol(covar.train)))
    lasso.fit <- big_spLogReg(G, y01.train=y, ind.train=ind.train, covar.train=covar.train,
                              dfmax=dfmax, ncores=ncores)
} else {
    cat(sprintf("Fitting sparse linear regression with %d observations, %d variants and %d covariates... ",
                length(y), ncol(G), ncol(covar.train)))
    lasso.fit <- big_spLinReg(G, y.train=y, ind.train=ind.train, covar.train=covar.train,
                              dfmax=dfmax, ncores=ncores)
}
cat("done.\n")

# Extract beta from each fold and combine them
cat("Extracting regression coefficients... ")
beta <- sapply(1:10, function(k) lasso.fit[[1]][k][[1]]$beta)

# Separate the coefficients of the genetic variants from the coefficients of the covariates
beta.variants <- beta[1:ncol(G),]
beta.covariates <- beta[(ncol(G)+1):nrow(beta),]

# Undo scaling of lasso coefficients
beta.variants <- beta.variants * scaling.factors
Beta <- cbind(tibble(CHR=Variants$CHR,
                     SNP=Variants$SNP, BP=Variants$BP, Knockoff=Variants$Knockoff),
                     as_tibble(beta.variants)) %>% as_tibble()
colnames(Beta) <- c("CHR", "SNP", "BP", "Knockoff", paste("K", seq(ncol(beta.variants)),sep=""))
Beta <- Beta %>%
    mutate(Z=(K1+K2+K3+K4+K5+K6+K7+K8+K9+K10)/10,
           Nonzero=(K1!=0)+(K2!=0)+(K3!=0)+(K4!=0)+(K5!=0)+(K6!=0)+(K7!=0)+(K8!=0)+(K9!=0)+(K10!=0)) %>%
    select(CHR, SNP, BP, Knockoff, Z)

# Extract the estimated coefficients
Lasso.res <- Beta %>%
    inner_join(Variants, by = c("CHR", "SNP", "BP", "Knockoff")) %>%
    filter(Z!=0) %>%
    arrange(desc(Z)) %>%
    select(CHR, SNP, BP, Z, Group, Knockoff)
cat("done.\n")

##################################
## Compute the test stastistics ##
##################################
cat("Computing test statistics... ")

# Compute the knockoff statistics
W.stats <- function(Z, knockoff) {
    importance <- abs(Z)
    z  <- sum(importance[which(knockoff==FALSE)], na.rm=T)
    zk <- sum(importance[which(knockoff==TRUE)], na.rm=T)
    w <- z-zk
}

Stats <- Lasso.res %>%
    select("CHR", "Group", "SNP", "BP", "Knockoff", "Z") %>%
    filter(Z!=0) %>%
    group_by(CHR, Group) %>%
    summarize(W = W.stats(abs(Z),Knockoff),
              Lead=which.max(abs(Z)), SNP.lead=SNP[Lead], BP.lead=BP[Lead],
              Size=n()) %>%
    ungroup() %>%
    arrange(desc(abs(W))) %>%
    select(CHR, Group, SNP.lead, BP.lead, Size, W) %>%
    filter(W!=0)
cat("done.\n")

# Give preview
Stats %>% print(n=20)

# Save results
stats.file <- sprintf("%s_stats.txt", out.basename)
Stats %>% write_delim(stats.file, delim=" ")
cat(sprintf("Test statistics written on:\n%s\n", stats.file))
