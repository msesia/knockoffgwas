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
local.source("utils_fdp.R")

load_annotations <- function(data_dir) {
  # Functional annotations
  annotations.func.file <- sprintf("%s/annotations/wgEncodeBroadHmmGm12878HMM.txt", data_dir)
  Annotations.func.raw <- read_tsv(annotations.func.file, 
                                   col_names = FALSE,
                                   col_types=cols(col_integer(), 
                                                  col_character(),
                                                  col_integer(),
                                                  col_integer(),
                                                  col_character(),
                                                  col_integer(),
                                                  col_character(),
                                                  col_integer(),
                                                  col_integer(),
                                                  col_integer()))
  names(Annotations.func.raw) = c("#bin", "chrom", "chromStart", "chromEnd",
                                  "name", "score", "strand", "thickStart", 
                                  "thickEnd", "itemRgb")
  
  Annotations.func.raw <- Annotations.func.raw %>% mutate(name=as.factor(name))
  
  # Remove weird chromosomes and convert colors
  valid.chrom <- paste("chr", seq(1,22), sep="")
  Annotations.func <- Annotations.func.raw %>%
    mutate(name=as.factor(name)) %>%
    filter(chrom %in% valid.chrom) %>% mutate(chrom=parse_number(as.character(chrom))) %>%
    mutate(itemR = as.integer(floor(itemRgb/256^2) %% 256), 
           itemG = as.integer(floor(itemRgb/256) %% 256), 
           itemB = as.integer(itemRgb %% 256)) %>%
    mutate(itemColor = rgb(red=itemR, blue=itemB, green=itemG, maxColorValue=255)) %>%
    select(c("#bin", "chrom", "chromStart", "chromEnd", "name", "score",
             "strand", "thickStart", "thickEnd", "itemR", "itemG", "itemB", "itemColor")) %>%
    arrange(chrom, chromStart)
  
  # Extract color map
  annotation.color.map <- Annotations.func %>% group_by(name, itemColor) %>% summarise() %>% 
    ungroup() %>%
    mutate(name.num=parse_number(as.character(name))) %>% 
    mutate(label=gsub("\\d+_", "",name), label=gsub(fixed("_"), " ",label)) %>%
    arrange(name.num)
  
  # Convert names to factors according to color maps
  Annotations.func <- Annotations.func %>%
    mutate(name=factor(name, levels=annotation.color.map$name, labels=annotation.color.map$name))
  
  # Gene annotations
  annotations.genes.file <- sprintf("%s/annotations/ncbiRefSeq.txt", data_dir)
  Annotations.genes.raw <- read_tsv(annotations.genes.file, 
                                    col_names = FALSE,
                                    col_types=cols(col_integer(),
                                                   col_character(),
                                                   col_character(),
                                                   col_character(),
                                                   col_integer(),
                                                   col_integer(),
                                                   col_integer(),
                                                   col_integer(),
                                                   col_integer(),
                                                   col_character(),
                                                   col_character(),
                                                   col_integer(),
                                                   col_character(),
                                                   col_character(),
                                                   col_character(),
                                                   col_character()))
  names(Annotations.genes.raw) = c("#bin", "name", "chrom", "strand", "txStart",
                                   "txEnd", "cdsStart", "cdsEnd", "exonCount",
                                   "exonStarts", "exonEnds", "score", "name2", 
                                   "cdsStartStat", "cdsEndStat", "exonFrames")
  
  # Remove weird chromosomes
  valid.chrom <- paste("chr", seq(1,22), sep="")
  Annotations.genes <- Annotations.genes.raw %>%
    filter(chrom %in% valid.chrom) %>% mutate(chrom=parse_number(as.character(chrom))) %>%
    arrange(chrom, `#bin`)
  
  # Split rows corresponding to same gene but different exons
  Exons <- Annotations.genes %>% 
    separate_rows(exonStarts, exonEnds, exonFrames, sep=",", convert=TRUE) %>%
    drop_na()
  
  # Pick the canonical transcripts
  # That is, for each unique "name2", keep only the rows corresponding to the "name" 
  # with the largest sum of exon lengths
  Exons.canonical <- Exons %>%
    mutate(exonLength=exonEnds-exonStarts) %>%
    group_by(name, name2) %>% summarise(Length=sum(exonLength)) %>%
    ungroup() %>% group_by(name2) %>% top_n(1, Length) %>%
    inner_join(Exons, by=c("name", "name2")) 
  
  annotations = c()
  annotations$Annotations.func = Annotations.func
  annotations$Exons.canonical = Exons.canonical
  
  return(annotations)
}

load_association_results <- function(data_dir, lmm_dir, phenotype){

    # Load knockoffs discovereries
    resolution.list <- paste("res", 0:6, sep="")
    phenotype.list <- c(phenotype)
    Params <- expand.grid(Resolution=resolution.list, Phenotype=phenotype.list) %>% as_tibble()
    Discoveries <- lapply(1:nrow(Params), function(idx) {
        resolution <- Params$Resolution[idx]
        phenotype <- Params$Phenotype[idx]
        knockoffs.file <- sprintf("%s/%s_%s_discoveries.txt", data_dir, phenotype, resolution)
        if(file.exists(knockoffs.file)) {
            Discoveries.knockoffs <- read_delim(knockoffs.file, delim=" ", col_types=cols()) %>%
                mutate(Phenotype=phenotype, Method="Knockoffs", Importance=W, Resolution=resolution) %>%
                select(-c("W", "Group"))
        } else {
            cat(sprintf("File %s not found!\n", knockoffs.file))
            Discoveries.knockoffs <- tibble()
        }
        return(Discoveries.knockoffs)
    })
    Discoveries <- do.call("rbind", Discoveries) %>% mutate(Resolution=as.character(Resolution))

    # Load knockoffs stats
    Stats <- lapply(1:nrow(Params), function(idx) {
        resolution <- Params$Resolution[idx]
        phenotype <- Params$Phenotype[idx]
        stats.file <- sprintf("%s/%s_%s_stats.txt", data_dir, phenotype, resolution)
        Stats <- read_delim(stats.file, delim=" ", col_types=cols()) %>% mutate(Resolution=resolution)
    })
    Stats <- do.call("rbind", Stats) %>% mutate(Resolution=as.character(Resolution))

    # Compute local FDP
    Discoveries <- local_fdp(Discoveries, Stats)

    # Load LMM p-values
    lmm.file <- sprintf("%s/%s_lmm.txt", lmm_dir, phenotype)
    if(file.exists(lmm.file)) {
        LMM <- read_delim(lmm.file, delim="\t", col_types=cols()) %>% as_tibble()
        if("P_BOLT_LMM" %in% colnames(LMM)) {
            LMM <- LMM %>% mutate(P=P_BOLT_LMM)
        } else {
            LMM <- LMM %>% mutate(P=P_BOLT_LMM_INF)
        }
    } else {
        LMM <- tibble()
    }

    # Load clumped LMM results
    lmm.file <- sprintf("%s/%s_lmm_clumped.txt", lmm_dir, phenotype)
    if(file.exists(lmm.file)) {
        LMM.clumped <- read_delim(lmm.file, delim=" ", col_types=cols()) %>%
            mutate(Phenotype=phenotype, Method="LMM", Importance=-log10(P), Resolution="GWAS") %>%
            select(-c("P"))
    } else {
        LMM.clumped <- tibble()
    }

    # Load low-res knockoff stats
    Pvalues <- Stats %>% filter(Resolution=="res2")

    association_results <- c()
    association_results$Discoveries <- Discoveries
    association_results$Stats <- Stats
    association_results$Pvalues <- Pvalues
    association_results$LMM <- LMM
    association_results$LMM.clumped <- LMM.clumped

    return(association_results)
}

plot_combined_state <- function(state, annotations){
    plot_combined(state$chr, state$window.left, state$window.right,
                  state$association_results$Discoveries,
                  state$association_results$LMM,
                  state$association_results$LMM.clumped,
                  Annotations.func=annotations$Annotations.func,
                  Exons.canonical=annotations$Exons.canonical,
                  highlight.gene=state$highlight.gene)
}

# From https://github.com/daqana/dqshiny
autocomplete_input <- function(
  id, label, options, value = "", width = NULL, placeholder = NULL,
  max_options = 0, hide_values = FALSE, create = FALSE, contains = FALSE
) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is needed to convert list of options into json!")
  }
  value <- shiny::restoreInput(id = id, default = value)
  js_opts <- jsonlite::toJSON(as.list(options), auto_unbox = TRUE)
  width <- shiny::validateCssUnit(width)
  if (length(value) == 0L) value <- ""
  shiny::div(
    class = "form-group shiny-input-container autocomplete",
    style = if (!is.null(width)) paste0("width: ", width, ";"),
    if (!is.null(label)) shiny::tags$label(label, `for` = id),
    shiny::tags$input(
      id = id, type = "text", class = "form-control", result = value,
      value = value, placeholder = placeholder, "data-options" = js_opts,
      "data-max" = max_options, "data-hide" = logical_js(hide_values),
      "data-create" = logical_js(create), "data-contains" = logical_js(contains)
    ),
    htmltools::htmlDependency(
      "autocomplete", "0.0.1", c(href = "dqshinyRes"),
      script = "js/autocomplete-binding.js", stylesheet = "css/autocomplete.css"
    )
  )
}

logical_js <- function(b) {
  tolower(isTRUE(b))
}