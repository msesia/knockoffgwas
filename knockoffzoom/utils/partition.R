#!/usr/bin/env Rscript

LocationOfThisScript = function() {
    # Solution from:
    # https://github.com/molgenis/molgenis-pipelines/wiki/How-to-source-another_file.R-from-within-your-R-script
    
    this.file = NULL
    # This file may be 'sourced'
    for (i in -(1:sys.nframe())) {
        if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
    }

    if (!is.null(this.file)) return(dirname(this.file))

    # But it may also be called from the command line
    cmd.args = commandArgs(trailingOnly = FALSE)
    cmd.args.trailing = commandArgs(trailingOnly = TRUE)
    cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
    res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)

    # If multiple --file arguments are given, R uses the last one
    res = tail(res[res != ""], 1)
    if (0 < length(res)) return(dirname(res))
    
    # Both are not the case. Maybe we are in an R GUI?
    return(NULL)
}
local.source <- function(script.name) {
    current.dir <- LocationOfThisScript()
    source(sprintf("%s/%s", current.dir, script.name))
}

# Default arguments (for debugging)
map.file    <- "../data/maps/genetic_map_chr22.txt"
bim.file   <- "../data/genotypes/example_chr22.bim"
qc.file    <- "../data/qc/variants_qc.txt"
out.file   <- "../data_processed/partitions/example_chr22.txt"

# Input arguments
args <- commandArgs(trailingOnly=TRUE)
map.file  <- as.character(args[1])
bim.file <- as.character(args[2])
qc.file  <- as.character(args[3])
out.file <- as.character(args[4])

## Values of resolution (measured in cM)
resolution.list <- c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 1)

# Load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(fastcluster))
#local.source("utils.R")

# Load list of variants
Variants <- read_delim(bim.file, delim="\t", col_names=c("CHR", "SNP", "X", "BP", "A1", "A2"), col_types=cols()) %>%
    select(CHR, SNP, BP)

# Remove variants that did not pass QC
variants.qc <- read_delim(qc.file, delim=" ", col_names="SNP", col_types=cols())
Variants <- Variants %>% inner_join(variants.qc, by="SNP")

# Load genetic map
Map <- read_delim(map.file, delim="\t", col_types=cols()) %>%
    transmute(CHR=parse_number(Chromosome), BP=`Position(bp)`, Rate=`Rate(cM/Mb)`, Map=`Map(cM)`)

## Cross-reference genetic map and list of variants
Map <- Map %>% inner_join(Variants, by = c("CHR", "BP"))

## Compute distances
Map$Map.diff <- c(0,pmax(diff(Map$Map),1e-6))
Map$Map.diff[is.na(Map$Map.diff)] <- 1e-6
Map$Map <- cumsum(Map$Map.diff)

## K-means clustering
n.clusters <- 5
clustering <- kmeans(Map$Map, n.clusters, nstart = 20)
clusters <- factor(clustering$cluster, levels=unique(clustering$cluster), labels=1:n.clusters) %>% as.integer()

if(FALSE) {
  tibble(position=Map$Map, cluster=clusters) %>% rownames_to_column() %>%
    mutate(rowname=as.integer(rowname), cluster=as.factor(cluster)) %>%
    ggplot(aes(x=rowname, y=position, color=cluster)) +
    geom_point() +
    theme_bw()
}

## Hierarchical clustering within each block
hc.list <- lapply(1:n.clusters, function(c) {
  c.idx <- which(clusters==c)
  map <- Map$Map[c.idx]
  D <- dist(map)
  hc <- fastcluster::hclust(D, method="complete")
  return(hc)
})

cut.trees <- function(hc.list, h) {
  ## Cut each tree at a certain height
  hc.clusters.list <- lapply(1:n.clusters, function(c) {
    hc.clusters <- cutree(hc.list[[c]], h=h)
    return(hc.clusters)
  })
  ## Combine list of clusters
  partial.sum <- 0
  hc.clusters <- c()
  for(c in 1:n.clusters) {
    hc.clusters <- c(hc.clusters, hc.clusters.list[[c]] + partial.sum)
    partial.sum <- partial.sum + max(hc.clusters.list[[c]])
  }
  return(hc.clusters)
}
 
## Cut dendrogram to define groups
Partitions <- Map %>% select(SNP, Map)
for(resolution.idx in 1:length(resolution.list)) {
  resolution <- resolution.list[resolution.idx]
  col.name <- sprintf("res_%d", length(resolution.list)-resolution.idx+1)
  hc.clusters <- cut.trees(hc.list, resolution)    
  Partitions[[col.name]] <- hc.clusters
}

## Compute average group sizes
group.sizes <- apply(Partitions[,-c(1,2)], 2, function(x) mean(table(x)))
cat(sprintf("Mean group sizes: \n"))
print(group.sizes)

if(FALSE) {
  col.name <- sprintf("res_%d", 1)
  hc.clusters <- Partitions[[col.name]]
  tibble(position=Map$Map, cluster=hc.clusters) %>% rownames_to_column() %>%
    mutate(rowname=as.integer(rowname), cluster=as.factor(cluster)) %>%
                                        #  filter(rowname>1100, rowname<1500) %>%
    ggplot(aes(x=rowname, y=position, color=cluster)) +
    geom_point() +
    theme_bw()
}

## Save results
Partitions %>% select(-Map) %>% write_delim(out.file, delim=" ")
cat(sprintf("Partitions written to: %s\n", out.file))
