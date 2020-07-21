library(tidyverse)
library(SNPknockG)
library(snpStats)

scratch <- "/scratch/groups/candes/ukbiobank_public/data"

chr.list <- seq(21,22)
p.keep <- 2000

###################
## Download data ##
###################

# Download the archived HAPMAP data
url.hap <- "https://mathgen.stats.ox.ac.uk/wtccc-software/rel21_poly/HapMap_rel21_b35_CEU.tgz"
url.map <- "https://mathgen.stats.ox.ac.uk/wtccc-software/recombination_rates/genetic_map_b35_CEU.tgz"
if(!file.exists(sprintf("%s/tmp/tmp.hap.tgz", scratch))) {
    download.file(url.hap,destfile=sprintf("%s/tmp/tmp.hap.tgz", scratch))
}
if(!file.exists(sprintf("%s/tmp/tmp.map.tgz", scratch))) {
    download.file(url.map,destfile=sprintf("%s/tmp/tmp.map.tgz", scratch))
}

dictionary <- NULL
variants <- NULL
for(chr in chr.list){
   hap.file <- sprintf("genotypes_chr%d_CEU_r21_nr_fwd_phased_all_by_snp_no_mono", chr)
   leg.file <- sprintf("genotypes_chr%d_CEU_r21_nr_fwd_legend_all_no_mono", chr)
   map.file <- sprintf("genetic_map_CEU_chr%d.txt", chr)
   # Unpack the archives
   untar(sprintf("%s/tmp/tmp.hap.tgz", scratch),files=c(hap.file,leg.file))
   untar(sprintf("%s/tmp/tmp.map.tgz", scratch),files=c(map.file))
   # Load the tables
   dictionary.raw <- read_table(hap.file,col_names=F)
   genmap <- read_delim(map.file, delim=" ", col_names=c("BP", "Rate", "Genetic"), skip=1) %>%
       mutate(CHR=chr) %>% select(CHR, everything())
   positions.raw <- read_delim(leg.file, delim=" ", col_names=c("SNP", "BP", "X0", "X1"), skip=1) %>%
       mutate(CHR=chr) %>% select(CHR, everything())
   # Make the problem smaller
   dictionary.raw <- dictionary.raw[1:(p.keep*2),]
   positions.raw <- positions.raw[1:(p.keep*2),]   
   # Cross-reference with positions
   variants.chr <- left_join(positions.raw, genmap, by = c("CHR","BP")) %>%
       mutate(Distance=c(1,diff(Genetic)))
   # Store results
   variants <- rbind(variants, variants.chr)
   dictionary <- rbind(dictionary, dictionary.raw)
}

# Fix reference alleles
dictionary[rowMeans(dictionary)>0.5,] = 1-dictionary[rowMeans(dictionary)>0.5,]

# Compute MAF
variants$MAF <- rowMeans(dictionary)
    
# Keep only as subset of common variants
variants.clean <- lapply(chr.list, function(chr) {
    variants.clean.chr <- variants %>%
        filter(CHR==chr, MAF>=0.2, BP>14e6) %>%
        arrange(BP) %>%
        head(n=round(p.keep/length(chr.list)))
    return(variants.clean.chr)
})
variants.clean <- do.call("rbind", variants.clean)
keep.idx <- which(variants$SNP %in% variants.clean$SNP)    
dictionary.clean <- dictionary[keep.idx,]

#############
## Fit HMM ##
#############

fp_path <- "/home/groups/candes/bin/fastphase"

# Fit an HMM, chromosome-by-chromosome
hmm.list <- lapply(chr.list, function(chr) {   
    # Extract indices for this chromosome
    chr.idx <- which(variants.clean$CHR==chr)
    # Store the dictionary as a matrix
    H <- as.matrix(dictionary.clean[chr.idx,])
    storage.mode(H) <- "integer"
    colnames(H) <- NULL
    # Write H as inp
    Hinp_file <- SNPknock.fp.writeX(t(H), phased=TRUE)
    # Call fastPhase and return the path to the parameter estimate files
    fp_out_path <- sprintf("%s/hmm/example_chr%d", scratch, chr)
    SNPknock.fp.runFastPhase(fp_path, Hinp_file, K=10, numit=15, phased=TRUE, out_path=fp_out_path)
    # Load the HMM
    r_file <- paste(fp_out_path, "_rhat.txt", sep="")
    theta_file <- paste(fp_out_path, "_thetahat.txt", sep="")
    alpha_file <- paste(fp_out_path, "_alphahat.txt", sep="")
    char_file <- paste(fp_out_path, "_origchars", sep="")
    hmm.chr <- SNPknock.fp.loadFit_hmm(r_file, alpha_file, theta_file, char_file, phased=TRUE)
    return(hmm.chr)
})
names(hmm.list) <- chr.list

#######################
## Generate new data ##
#######################

# Generate haplotype sequences, chromosome-by-chromosome
n.samples <- 1000
n <- 2 * n.samples
cat(sprintf("Generating haplotypes for %d subjects... ", n.samples))
H.list <- lapply(chr.list, function(chr) {
    # Extract hmm for this chromosome
    hmm <- hmm.list[[as.character(chr)]]
    # Extract indices for this chromosome
    H.chr <- SNPknock.models.sampleHMM(hmm$pInit, hmm$Q, hmm$pEmit, n=n)
    return(H.chr)
})
H <- do.call("cbind", H.list)
cat("done. \n")

# Convert haplotypes to genotypes
X <- H[seq(1,2*n.samples,2),] + H[seq(2,2*n.samples,2),]

##########################
## Save genotypes (BED) ##
##########################

# Choose sex of individuals
X.sex <- 1+rbinom(nrow(X),1,0.5)

# Variant names
colnames(X) <- variants.clean$SNP

# Sample IDs (make sure they all have the same length)
rownames(X) <- nrow(X) + seq(1, nrow(X)) - 1

# Save genotypes into PLINK format, chromosome-by-chromosome
for(chr in chr.list){
    cat(sprintf("Saving genotypes on chromosome %d.\n", chr))
    # Extract indices for this chromosome
    chr.idx <- which(variants.clean$CHR==chr)
    # Create SnpMatrix object
    X.SnpMatrix <- as(X[,chr.idx], "SnpMatrix")
    # Save SnpMatrix
    file.base <- sprintf("%s/genotypes/example_chr%d", scratch, chr)
    write.plink(file.base, snp.major = TRUE, X.SnpMatrix, id=rownames(X), sex=X.sex,
                chromosome=variants.clean$CHR[chr.idx], position=variants.clean$BP[chr.idx])
}

############################
## Save haplotypes (HAPS) ##
############################

for(chr in chr.list){
    cat(sprintf("Saving haplotypes on chromosome %d.\n", chr))
    # Extract indices for this chromosome
    chr.idx <- which(variants.clean$CHR==chr)
    # Write haplotypes in HAPS format
    haps.file <- sprintf("%s/tmp/example_chr%d.haps", scratch, chr)
    H[,chr.idx] %>% t() %>% as_tibble() %>% write_delim(haps.file, delim=" ", col_names=F)
    cat(sprintf("Results written on: %s\n", haps.file))
    # Write legend
    legend.chr <- variants.clean %>% filter(CHR==chr) %>% select(SNP, BP, X0, X1)
    colnames(legend.chr) <- c("id", "position", "a0", "a1")
    leg.file <- sprintf("%s/tmp/example_chr%d.leg", scratch, chr)
    legend.chr %>% write_delim(leg.file, delim=" ")
    cat(sprintf("Results written on: %s\n", leg.file))
    # Write sample file
    sample.chr <- tibble(ID_1=c(0,rownames(X)), ID_2=c(0,rownames(X)), missing=0, sex=c("D",X.sex))
    sample.file <- sprintf("%s/tmp/example_chr%d.sample", scratch, chr)
    sample.chr %>% write_delim(sample.file, delim=" ")
    cat(sprintf("Results written on: %s\n", sample.file))
}


Sigma <- cov2cor(cov(X))

plot(abs(Sigma[400,]))
