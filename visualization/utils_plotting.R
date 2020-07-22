suppressMessages(library(ggrepel))
suppressMessages(library(latex2exp))
suppressMessages(library(egg))
suppressMessages(library(grid))
suppressMessages(library(cowplot))

font.size <- 18
title.font.size <- 18
axis.font.size <- 18
legend.font.size <- 15
dot.size <- 2.5

bp.labeler <- function(x) {
    round(x*1e-6,3)
}

find_chr_boundaries <- function(association_results, chr) {
    # Initialize boundaries
    min.BP <- NULL
    max.BP <- NULL

    # Look at the range of LMM p-values
    if(nrow(association_results$LMM)>0) {
        if(nrow(filter(association_results$LMM, CHR==chr))>0) {
            min.BP <- min(filter(association_results$LMM, CHR==chr)$BP)
            max.BP <- max(filter(association_results$LMM, CHR==chr)$BP)
        } else {
            if(nrow(association_results$Stats)>0) {
                if(nrow(filter(association_results$Stats, CHR==chr))>0) { 
                    min.BP <- min(filter(association_results$Stats, CHR==chr)$BP.min)
                    max.BP <- max(filter(association_results$Stats, CHR==chr)$BP.max)
                } else {
                    min.BP <- 0
                    max.BP <- 0.1
                }
            } else {
                min.BP <- 0
                max.BP <- 0.1
            }
        }
    } else {
        if(nrow(association_results$Stats)>0) {
            if(nrow(filter(association_results$Stats, CHR==chr))>0) { 
                min.BP <- min(filter(association_results$Stats, CHR==chr)$BP.min)
                max.BP <- max(filter(association_results$Stats, CHR==chr)$BP.max)
            } else {
                min.BP <- 0
                max.BP <- 0.1
            }
        } else {
            min.BP <- 0
            max.BP <- 0.1
        }
    }

    # Return boundaries
    boundaries <- c()
    boundaries$min.BP <- min.BP
    boundaries$max.BP <- max.BP
    return(boundaries)
}

place_segments <- function(Segments, gap=5e4, verbose=FALSE) {

  # Define function that checks whether a new segment would fit in an existing row
  segment_fits <- function(start, end, row, gap) {
    if(length(row)==0) {
      return(TRUE)
    }
    for(segment in row) {
      if(start<=segment[2]+gap && end>=segment[1]-gap){
        return(FALSE)
      }
    }
    return(TRUE)
  }

  # Count number of segments
  n.segments <- nrow(Segments)
  segment.heights <- rep(NA, n.segments)
  segment.rows <- list()

  n.rows <- 0

  for(j in 1:n.segments) {
    start <- Segments$start[j]
    end <- Segments$end[j]

    if(length(segment.rows)==0) {
      # Add the first segment to the first row
      segment.rows[[1]] <- list()
      segment.rows[[1]][[1]] <- c(start,end)
      segment.heights[j] <- 1
      n.rows <- n.rows + 1
    } else {
      would.fit <- sapply(1:length(segment.rows), function(g) {
        segment_fits(start, end, segment.rows[[g]], gap=gap)
      })
      if(any(would.fit)) {
        # Add segment to existing row
        k <- min(which(would.fit))
        #print(length(segment.rows[[k]]))
        row.length <- length(segment.rows[[k]])
        if(verbose) cat(sprintf("Appending (%d,%d) to row %d \n",start,end,k))
        segment.rows[[k]][[row.length+1]] <- c(start,end)
        segment.heights[j] <- k
      } else {
        # Add segment to new row
        n.rows <- n.rows + 1
        if(verbose) cat(sprintf("No room to append (%d,%d). Creating new row %d. \n",start,end,n.rows))
        segment.rows[[n.rows]] <- list()
        segment.rows[[n.rows]][[1]] <- c(start,end)
        segment.heights[j] <- n.rows
      }
    }
  }

  # Assign heights to each segment
  Segments$Height <- segment.heights

  # Return modified data
  return(Segments)
}

# Minimal theme
theme_minimal <- theme_bw() +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank(),
          panel.border=element_blank(),
          panel.grid.major.x = element_line(size = 0.2, colour = "darkgray"),
          panel.grid.minor.x = element_line(size = 0.1, colour = "darkgray")
          )

plot_pvalues <- function(window.chr, window.left, window.right, LMM, LMM.clumped,
                         p.significant=5e-8, p.max=1e-0) {

    # Extract LMM p-values within this window
    LMM.window <- LMM %>% filter(CHR==window.chr, BP<=window.right, BP>=window.left) %>%
        left_join(LMM.clumped, by = c("SNP", "CHR", "BP")) %>%
        mutate(BP.lead=factor(BP.lead))

    #cat(sprintf("There are %d LMM pvalues within this window, %d of which are significant.\n",
    #            nrow(LMM.window), sum(LMM.window$P<p.significant)))

    # Extract clumps that overlap with this window
    LMM.clumped.window <- LMM.clumped %>% filter(CHR==window.chr) %>%
        group_by(CHR, SNP.lead, BP.lead) %>%
        summarise(BP.min=min(BP), BP.max=max(BP)) %>%
        filter(BP.min<=window.right, BP.max>=window.left) %>%
        select(CHR, SNP.lead) %>%
        inner_join(LMM.clumped, by = c("CHR", "SNP.lead"))

    #cat(sprintf("There are %d LMM clumps that overlap with this window.\n",
    #            length(unique(LMM.clumped.window$SNP.lead))))

    # Manhattan plot
    if(nrow(LMM.window)>0) {

        # Significance level for pvalues
        Window.nominal <- LMM.window %>% mutate(Importance = -log10(p.significant))

        LMM.window$BP.lead <- round((LMM.window$BP.lead %>% as.character %>% parse_number)/1e6,3) %>% factor
        clump.lead.snps <- unique(LMM.clumped.window$SNP.lead)
        
        if(all(is.na(LMM.window$BP.lead))) {
            p.manhattan <- LMM.window %>%
                filter(P<p.max) %>%
                mutate(P=pmax(1e-300,P)) %>%
                ggplot(aes(x=BP, y=-log10(P))) +
                geom_point(color="black", alpha=0.25, size=dot.size)
        } else {
            p.manhattan <- LMM.window %>%
                filter(P<p.max) %>%
                mutate(SNP.lead=factor(SNP.lead, levels=clump.lead.snps, labels=clump.lead.snps)) %>%
                mutate(P=pmax(1e-300,P)) %>%
                ggplot(aes(x=BP, y=-log10(P), color=SNP.lead, alpha=is.na(SNP.lead))) +
                geom_point(size=dot.size) +
                scale_colour_discrete(na.value = "black", name="Lead SNP (LMM)", guide=FALSE)
        }
    } else {
        p.manhattan <- ggplot(tibble()) + geom_blank()
    }

    # Vertical limits for Manhattan plots
    manhattan.y.max <- NA
    manhattan.y.breaks <- c(0,7.3,NA)

    p.manhattan <- p.manhattan +
        geom_hline(yintercept=-log10(p.significant), linetype="dashed", color = "black") +
        scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 1), guide=FALSE) +
        coord_cartesian(xlim = c(window.left,window.right)) +
        scale_x_continuous(expand=c(0.01,0.01), labels=bp.labeler) +
        scale_y_continuous(limits=c(1.8,manhattan.y.max), breaks=manhattan.y.breaks, trans="sqrt") +
        theme_bw() +
        xlab(sprintf("Chromosome %d (Mb)", window.chr)) + ylab(TeX("$-\\log_{10}(p)$")) +
        theme(text = element_text(size=font.size),
              plot.title = element_text(size=title.font.size),
              axis.title = element_text(size=axis.font.size),
              legend.text = element_text(size=legend.font.size),
              legend.title = element_text(size=legend.font.size),
              panel.border = element_blank(),
              panel.grid.major.x = element_line(size = 0.2, colour = "darkgray"),
              panel.grid.minor.x = element_line(size = 0.1, colour = "darkgray")
              ) +
        ggtitle(sprintf("Manhattan plot (BOLT-LMM)"))

    # Plot clumped LMM discoveries
    if(nrow(LMM.clumped.window)>0) {
        LMM.clumps.window <- LMM.clumped.window %>% group_by(CHR, SNP.lead, BP.lead) %>%
            summarise(BP.min=min(BP), BP.max=max(BP)) %>%
            ungroup() %>%
            mutate(BP.lead=factor(BP.lead)) %>%
            mutate(SNP.lead=factor(SNP.lead, levels=clump.lead.snps, labels=clump.lead.snps)) %>%
            mutate(start=BP.min, end=BP.max) %>% mutate(width=end-start) %>% arrange(desc(width)) %>%
            place_segments(gap=1e4) %>%
            mutate(Height=-Height)

        p.clumped <- LMM.clumps.window %>%
            ggplot() +
            geom_segment(aes(x=BP.min, y=Height,
                             xend=BP.max, yend=Height, color=SNP.lead)) +
            geom_segment(aes(x=BP.min, y=Height-0.4, xend=BP.min, yend=Height+0.4, color=SNP.lead)) +
            geom_segment(aes(x=BP.max, y=Height-0.4, xend=BP.max, yend=Height+0.4, color=SNP.lead)) +
            scale_colour_discrete(na.value = "black", name="Lead SNP (LMM)") +
                    coord_cartesian(xlim = c(window.left,window.right)) +
            scale_x_continuous(expand=c(0.01,0.01), labels=bp.labeler) +
            scale_y_continuous(breaks = NULL) +
            ylab("") + xlab("") +
            theme_void() +
            theme(text = element_text(size=font.size),
                  plot.title = element_text(size=title.font.size),
                  axis.title = element_text(size=axis.font.size),
                  legend.text = element_text(size=legend.font.size),
                  legend.title = element_text(size=legend.font.size),
                  legend.key.size = unit(0.5,"line")) +
            guides(color=guide_legend(ncol=2))

      # Determine whether we would list all lead SNPs in the legend
      if(length(unique(LMM.clumps.window$SNP.lead))>20) {
        p.clumped <- p.clumped +
          theme(legend.position = "none")
      }
    } else {
      LMM.clumps.window <- tibble()
      p.clumped <- ggplot(tibble()) + geom_blank()
    }

    # Return plot objects
    plots <- c()
    plots$manhattan <- p.manhattan
    plots$clumped <- p.clumped
    return(plots)
}

plot_chicago <- function(window.chr, window.left, window.right, Discoveries) {
    
    # Extract knockoff discoveries within this window
    if(!is.null(Discoveries)) {
        Knockoffs.window <- Discoveries %>% filter(Method=="Knockoffs") %>%
            filter(CHR==window.chr, BP.min<=window.right, BP.max>=window.left)
        #cat(sprintf("There are %d knockoff discoveries within this window.\n", nrow(Knockoffs.window)))
    } else {
        Knockoffs.window <- tibble()
    }
    
    # Plot knockoff discoveries
    resolution.list <- c("res0", "res1", "res2", "res3", "res4", "res5", "res6") %>% rev
    resolution.heights <- seq(length(resolution.list))
    names(resolution.heights) <- resolution.list
    #resolution.labels <- paste(parse_number(resolution.list), "\\%", sep="")
    resolution.labels <- c("single-SNP", "6", "22", "38", "75", "149", "358") %>% rev

    if(nrow(Knockoffs.window)>0) {
        p.knockoffs <- Knockoffs.window %>%
            mutate(Resolution=as.character(Resolution)) %>%
            mutate(Height=resolution.heights[Resolution]) %>%
            ggplot() +
            geom_rect(aes(xmin=BP.min, xmax=BP.max, ymin=Height-0.5, ymax=Height+0.5, fill=1-FDP.local),
                      color="black") +
            scale_fill_gradient(name="Local FDP (estimated)",
                                #low="red1", high="dodgerblue1",
                                low="gray98", high="gray30",
                                limits=c(0.5, 1),
                                breaks=c(0.5,0.75,1),
                                space="Lab", na.value="gray", guide="colourbar",
                                labels=function(x){1-x})
    } else {
        p.knockoffs <- ggplot(tibble()) + geom_blank()
    }

    p.knockoffs <- p.knockoffs +
        ylab("Resolution (Mb)") + xlab("") +
        coord_cartesian(xlim = c(window.left,window.right)) +
        scale_x_continuous(expand=c(0.01,0.01), labels=bp.labeler) +
        scale_y_continuous(limits=c(0.5,max(resolution.heights)+0.5),
                           labels=resolution.labels, breaks=resolution.heights) +
        ggtitle("Chicago plot (KnockoffZoom)") +
        theme_bw() +
        theme(panel.grid.minor.y = element_blank(),
              axis.line=element_blank(),
              axis.title.x=element_blank(),
              panel.border=element_blank(),
              panel.grid.major.x = element_line(size = 0.2, colour = "darkgray"),
              panel.grid.minor.x = element_line(size = 0.1, colour = "darkgray"),
              text = element_text(size=font.size),
              axis.title.y = element_text(size=title.font.size),
              plot.title = element_text(size=title.font.size),
              legend.text = element_text(size=legend.font.size),
              legend.title = element_text(size=legend.font.size),
              legend.key.height = unit(0.75,"line")
              )

    return(p.knockoffs)
}

plot_annotations <- function(window.chr, window.left, window.right, Annotations.func) {
    # Extract color map
    annotation.color.map <- Annotations.func %>% group_by(name, itemColor) %>% summarise() %>%
        ungroup() %>%
        mutate(name.num=parse_number(as.character(name))) %>%
        mutate(label=gsub("\\d+_", "",name), label=gsub(fixed("_"), " ",label)) %>%
        mutate(label=as.factor(label)) %>%
        arrange(name.num)

    # Convert names to factors according to color maps
    Annotations.func <- Annotations.func %>%
        mutate(label=gsub("\\d+_", "",name), label=gsub(fixed("_"), " ",label)) %>%
        mutate(label=factor(label, levels=annotation.color.map$label, labels=annotation.color.map$label))

    # Plot functional annotations
    myColors <- annotation.color.map$itemColor
    names(myColors) <- annotation.color.map$label

    # Select functional annotations within this window
    Functional.window <- Annotations.func %>%
        filter(chrom==window.chr, chromStart<=window.right, chromEnd>=window.left)
    #cat(sprintf("There are %d functional annotations within this window.\n",
    #            nrow(Functional.window)))

    # Make plot
    if(nrow(Functional.window)>0) {
        p.functional <- Functional.window %>%
            mutate(chromStart=chromStart, chromEnd=chromEnd) %>%
            ggplot() +
            geom_rect(aes(xmin=chromStart, xmax=chromEnd, ymin=0.5, ymax=1.5, fill=label)) +
            ylab("") + xlab("") +
            coord_cartesian(xlim = c(window.left,window.right)) +
            scale_x_continuous(expand=c(0.01,0.01), labels=bp.labeler) +
            scale_color_manual(values=myColors, guide=FALSE) +
            scale_fill_manual(values=myColors, name="Variant annotation") +
            ggtitle("Functional annotations") +
            theme_void() +
            theme(legend.key.size = unit(0.5,"line"),
                  text = element_text(size=font.size),
                  plot.title = element_text(size=title.font.size),
                  axis.title = element_text(size=axis.font.size),
                  legend.text = element_text(size=legend.font.size),
                  legend.title = element_text(size=legend.font.size),
                  )+
            guides(fill=guide_legend(ncol=2))
    } else {
        p.functional <- ggplot(tibble()) + geom_blank()
    }

    # Return plot
    return(p.functional)
}

plot_genes <- function(window.chr, window.left, window.right, Exons.canonical,
                       highlight.gene=NULL, max.gene.rows=10) {

    # Select exons within this windows
    Exons.window <- Exons.canonical %>%
        filter(chrom==window.chr, txStart<=window.right, txEnd>=window.left)
    #cat(sprintf("There are %d exons within this window, divided into %d genes.\n",
    #            nrow(Exons.window), length(unique(Exons.window$name2))))


    # Find out how many genes there are and determine whether we would plot all of them
    Genes.window <- Exons.window %>% group_by(name, name2, strand) %>%
        summarise(txStart=min(txStart), txEnd=max(txEnd)) %>%
        mutate(txStart=max(txStart, window.left), txEnd=min(txEnd, window.right)) %>%
        mutate(start=txStart, end=txEnd)

    # Highlight special gene, if available
    if(nrow(Genes.window) > 0) {
        if(is.null(highlight.gene)) {
            Genes.window$highlight <- FALSE
        } else {
            Genes.window <- Genes.window %>% mutate(highlight = (name2==highlight.gene))
        }

        # Do not attempt to place more than a max number of genes
        n.genes <- length(unique(Genes.window$name))
        max.genes <- 50
        if(n.genes<=max.genes) {
            plot.genes <- TRUE
        } else {
            plot.genes <- FALSE
            genes.toshow <- unique(Genes.window$name)[1:max.genes]
            Genes.window <- Genes.window %>% filter((name %in% genes.toshow) || (highlight))
        }

        # Sort the genes by length, making sure that the highlighted gene is first;
        # then, place them on different rows so they don't overlap
        Genes.window <- Genes.window %>% mutate(width = end-start) %>% arrange(desc(highlight), desc(width)) %>%
            place_segments()

        # Remove genes that do not fit in 3 rows
        if(length(unique(Genes.window$Height))>max.gene.rows) {
            plot.genes <- FALSE
            Genes.window <- Genes.window %>% filter(Height<=max.gene.rows)
        }
    } else{
        plot.genes <- FALSE
        n.genes <- 0
    }        
    
    # Plot exons and genes, if they fit
    if(plot.genes) {
        # Rescale the gene heights
        Genes.window <- Genes.window %>% mutate(Height=-(Height-0.5)) %>%
            mutate(txCenter=(txStart+txEnd)/2, Height=Height/2)

        # Extract genes within this window
        strand.levels <- c("+", "-")
        strand.labels <- c("->", "<-")
        Genes.window <- Genes.window %>%
            mutate(strand.latex=factor(strand, labels=strand.labels, levels=strand.levels)) %>%
            mutate(strand.latex=as.character(strand.latex))

        p.genes <- Exons.window %>%
            filter(exonStarts>=window.left, exonEnds<=window.right) %>%
            inner_join(Genes.window %>% select(name, name2, Height), by = c("name", "name2")) %>%
            ggplot() +
            geom_rect(aes(xmin=exonStarts, xmax=exonEnds, ymin=Height-0.1, ymax=Height+0.1),
                      alpha=1, color="black", fill="black") +
            geom_segment(data=Genes.window, aes(x=txStart, y=Height, xend=txEnd, yend=Height, group=name2),
                         color="black") +
            geom_label_repel(data=Genes.window, aes(x=txCenter, y=Height+0.1, label="NA"),
                             size=1, alpha=0, max.iter=0) +
            geom_label_repel(data=Genes.window, aes(x=txCenter, y=Height-0.1, label="NA"),
                             size=1, alpha=0, max.iter=0) +
            geom_label_repel(data=Genes.window, aes(x=txStart, y=Height+0.1, label="NA"),
                             size=1, alpha=0, max.iter=0) +
            geom_label_repel(data=Genes.window, aes(x=txStart, y=Height-0.1, label="NA"),
                             size=1, alpha=0, max.iter=0) +
            geom_label_repel(data=Genes.window, aes(x=txEnd, y=Height+0.1, label="NA"),
                             size=1, alpha=0, max.iter=0) +
            geom_label_repel(data=Genes.window, aes(x=txEnd, y=Height-0.1, label="NA"),
                             size=1, alpha=0, max.iter=0) +
            geom_label_repel(data=Genes.window, aes(x=txCenter, y=Height,
                                                    label=paste(name2,strand.latex,sep=" "),
                                                    fill=highlight),
                             size=5, direction="both", force=1, max.iter=2000,
                             box.padding=0.1,
                             point.padding=1,
                             label.padding=0.15,
                             segment.color = 'grey50', segment.alpha=0.5, seed=2019) +
            ylab("") + xlab("") +
            coord_cartesian(xlim = c(window.left,window.right)) +
            scale_x_continuous(expand=c(0.01,0.01), labels=bp.labeler) +
            scale_y_continuous(expand=c(0.01,0.01), breaks = NULL,
                               limits=c(min(Genes.window$Height)-0.25,max(Genes.window$Height)+0.25)) +
            scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "yellow"), guide=FALSE) +
            theme_void() +
            ggtitle("Genes") +
            theme(text = element_text(size=font.size),
                  plot.title = element_text(size=title.font.size),
                  axis.title = element_text(size=axis.font.size),
                  panel.grid.major.x = element_line(size = 0.2, colour = "darkgray"),
                  panel.grid.minor.x = element_line(size = 0.1, colour = "darkgray")
                  )
    } else {
        if(n.genes==0) {
            gene.title <- sprintf("There are no genes in this region.")
        } else {
            gene.title <- sprintf("There are %d genes in this region. Zoom in to see them.", n.genes)
        }
        p.genes <- ggplot(tibble()) + geom_blank() + ggtitle(gene.title)
    }

    # Return plot
    return(p.genes)
}

plot_combined <- function(window.chr, window.left, window.right, Discoveries, LMM, LMM.clumped,
                          Annotations.func=NULL, Exons.canonical=NULL,
                          highlight.gene=NULL, max.gene.rows=10) {

    # Make sure that the window is not empty
    if(window.right<=window.left) {
        return(ggplot(tibble()) + geom_blank())
    }

    # Make Chicago plot with KnockoffZoom discoveries
    p.knockoffs <- plot_chicago(window.chr, window.left, window.right, Discoveries)

    # Make plots with LMM p-values
    if(nrow(LMM)>0) {
        p.lmm <- plot_pvalues(window.chr, window.left, window.right, LMM, LMM.clumped)
    } else {
        p.lmm <- c()
        p.lmm$manhattan <- ggplot(tibble()) + geom_blank()
        p.lmm$clumped <- ggplot(tibble()) + geom_blank()
    }

    # Plot functional annotations
    if(!is.null(Annotations.func)) {
        p.functional <- plot_annotations(window.chr, window.left, window.right, Annotations.func)
    } else {
        p.functional <- ggplot(tibble()) + geom_blank()
    }

    # Plot genes
    if(!is.null(Exons.canonical)) {
        p.genes <- plot_genes(window.chr, window.left, window.right, Exons.canonical,
                              highlight.gene=highlight.gene, max.gene.rows=max.gene.rows)
    } else {
        p.genes <- ggplot(tibble()) + geom_blank()
    }

    # Determine relative heights of each subplot
    height.manhattan <- 0.75
    height.clumps <- 0.3
    height.knockoffs <- 1.25
    height.functional <- 0.3
    height.genes <- 0.25*max.gene.rows
    heights <- c(height.manhattan, height.clumps, height.knockoffs, height.functional, height.genes)

    # Convert the plot objects for placement
    debug.lines <- FALSE
    g1 <-  ggplotGrob(p.lmm$manhattan)
    g2 <-  ggplotGrob(p.lmm$clumped + theme(legend.position="none"))
    g3 <-  ggplotGrob(p.knockoffs + theme(legend.position="none"))
    g4 <-  ggplotGrob(p.functional + theme(legend.position="none"))
    g5 <-  ggplotGrob(p.genes)
    fg1 <- gtable_frame(g1, width = unit(1, "null"), height = unit(heights[1], "null"), debug = debug.lines)
    fg2 <- gtable_frame(g2, width = unit(1, "null"), height = unit(heights[2], "null"), debug = debug.lines)
    fg3 <- gtable_frame(g3, width = unit(1, "null"), height = unit(heights[3], "null"), debug = debug.lines)
    fg4 <- gtable_frame(g4, width = unit(1, "null"), height = unit(heights[4], "null"), debug = debug.lines)
    fg5 <- gtable_frame(g5, width = unit(1, "null"), height = unit(heights[5], "null"), debug = debug.lines)

    # Combine the main plots
    fg.l <- gtable_frame(gtable_rbind(fg1, fg2, fg3, fg4, fg5),
                         width = unit(4, "null"), height = unit(1, "null"))

    # Extract the legends
    g6 <- ggplotGrob(ggplot())
    try(g6 <- ggplotGrob(ggplotify::as.ggplot(get_legend(p.lmm$clumped))+
                         theme(text = element_text(size=legend.font.size))
                         ), silent=TRUE)
    g7 <- ggplotGrob(ggplot())
    try(g7 <- ggplotGrob(ggplotify::as.ggplot(get_legend(p.knockoffs))+
                         theme(text = element_text(size=legend.font.size))
                         ), silent=TRUE)
    g8 <- ggplotGrob(ggplot())
    try(g8 <- ggplotGrob(ggplotify::as.ggplot(get_legend(p.functional))+
                         theme(text = element_text(size=legend.font.size))
                         ), silent=TRUE)

    # Combine the legends
    g0 <- ggplotGrob(ggplot())
    fg6 <- gtable_frame(g6, width = unit(1, "null"), height = unit(0.75, "null"), debug = debug.lines)
    fg7 <- gtable_frame(g7, width = unit(1, "null"), height = unit(1, "null"), debug = debug.lines)
    fg8 <- gtable_frame(g8, width = unit(1, "null"), height = unit(1, "null"), debug = debug.lines)
    fg00 <- gtable_frame(g0, width = unit(1, "null"), height = unit(1, "null"), debug = debug.lines)
    fg.r <- gtable_frame(gtable_rbind(fg6,fg7,fg8,fg00),
                         width = unit(1, "null"), height = unit(1, "null"))

    # Combine main plots and legends
    grid.newpage()
    combined <- gtable_frame(gtable_cbind(fg.l, fg.r),
                         width = unit(1, "null"),
                         height = unit(1, "null"))
    p.final <- ggplotify::as.ggplot(combined)

    # Return complete plot
    return(p.final)
}
