suppressMessages(library(tidyverse))
suppressMessages(library(gridExtra))
suppressMessages(library(shiny))
source("utils_clumping.R")
source("utils_manhattan.R")
source("utils_plotting.R")
source("utils_shiny.R")

data_dir <- "../data"
res_dir <- "../results"
chromosomes <- 21:22
phenotype <- "example"

annotations <- load_annotations(data_dir)

# Load associations
source("utils_shiny.R")
association_results <- load_association_results(res_dir, phenotype)

source("utils_manhattan.R")
plot_manhattan_knockoffs(association_results$LMM, association_results$Stats, ytrans="identity")

#Make annotation plot
source("utils_plotting.R")
window.chr <- 21
window.center <- 41.0
window.width <- 2
window.boundaries <- find_chr_boundaries(association_results, window.chr)
window.left <- window.boundaries$min.BP
window.right <- window.boundaries$max.BP
plot_combined(window.chr, window.left, window.right,
              association_results$Discoveries,
              association_results$LMM,
              association_results$LMM.clumped,
              highlight.gene="KIF1B", max.gene.rows=10)

Discoveries <- association_results$Discoveries
LMM <- association_results$LMM
Clumped <- association_results$LMM.clumped
Annotations.func <- annotations$Annotations.func
Exons.canonical <- annotations$Exons.canonical


#Make annotation plot
source("utils_plotting.R")
window.chr <- 22
window.center <- 41.0
window.width <- 0.1
window.left <- 1e6*max(0, window.center - window.width)
window.right <- 1e6*min(window.center + window.width,
                    1e-6*max(filter(association_results$LMM, CHR==window.chr)$BP))
plot_combined(window.chr, window.left, window.right,
                 association_results$Discoveries,
                 association_results$LMM,
                 association_results$LMM.clumped,
                 Annotations.func=annotations$Annotations.func,
                 Exons.canonical=annotations$Exons.canonical,
                 highlight.gene="KIF1B", max.gene.rows=10)



window.chr <- 1
source("utils_manhattan.R")
plot_manhattan_knockoffs(association_results$LMM, association_results$Stats, ytrans="identity", chr=window.chr)
