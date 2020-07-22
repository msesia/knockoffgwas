suppressMessages(library(cowplot))
suppressMessages(library(latex2exp))

font.size <- 20
dot.size <- 2.5

knockoff.threshold <- function(W, fdr=0.10, offset=1) {
  if(offset>1 | offset<0) {
    stop('Input offset must be between 0 or 1')
  }
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}
                 
pre_process <- function(gwas) {
    gwas.don <- gwas %>%
        # Compute chromosome size
        group_by(CHR) %>% 
        summarise(chr_len=max(BP)) %>% 

        # Calculate cumulative position of each chromosome
        mutate(tot=cumsum(chr_len)-chr_len) %>%
        select(-chr_len) %>%

        # Add this info to the initial dataset
        left_join(gwas, ., by=c("CHR"="CHR")) %>%

        # Add a cumulative position of each SNP
        arrange(CHR, BP) %>%
        mutate( BPcum=BP+tot)

    return(gwas.don)
}

manhattan_simple <- function(gwas.don, axisdf, limit.left, limit.right,
                             ytrans="log10", yintercept=-log10(5e-8)) {
    # Rescale things to Mb
    limit.left <- limit.left/1e6
    limit.right <- limit.right/1e6
    axisdf$center <- axisdf$center/1e6
    gwas.don$BP <- gwas.don$BP/1e6
    gwas.don$BPcum <- gwas.don$BPcum/1e6

    # Make LMM Manhattan plot
    y.max <- 1.1 * max(-log10(pmax(1e-300,gwas.don$P)))
    p.manhattan <- gwas.don %>%
        mutate(P=pmax(1e-300,P)) %>%
        filter(P<1e-3) %>%
        ggplot(aes(x=BPcum, y=-log10(P))) +
    
        # Show all points
        geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=dot.size) +
        scale_color_manual(values = rep(c("darkgrey", "black"), 22 )) +

        # Show significance threshold
        geom_hline(yintercept=yintercept, linetype="dashed", color = "red") +

        # Custom axes:
        scale_x_continuous(label=axisdf$CHR, breaks=axisdf$center, limits=c(limit.left,limit.right),
                           expand=c(0.01,0.01)) +
        scale_y_continuous(trans=ytrans, expand = c(0,0), limits=c(3,y.max)) + 
        xlab("Chromosome") + ylab("Importance") +
    
        # Custom the theme:
        theme_bw() +
        theme(legend.position="none", panel.border = element_blank(), 
              panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
              panel.grid.major = element_line(size = 0.2, colour = "darkgray"),
              panel.grid.minor = element_line(size = 0.1, colour = "lightgray"),
              text = element_text(size=font.size)
              )

    return(p.manhattan)
}

manhattan_knock <- function(gwas.don, axisdf, limit.left, limit.right,
                            ytrans="log10", yintercept=0, use.Mb=FALSE) {

    # Rescale things to Mb
    limit.left <- limit.left/1e6
    limit.right <- limit.right/1e6
    axisdf$center <- axisdf$center/1e6
    gwas.don$BP <- gwas.don$BP/1e6
    gwas.don$BPcum <- gwas.don$BPcum/1e6
    
    # Make LMM Manhattan plot
    y.max <- 1.1 * max(gwas.don$W)
    p.manhattan <- gwas.don %>%        
        ggplot(aes(x=BPcum, y=W)) +
    
        # Show all points
        geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=dot.size) +
        scale_color_manual(values = rep(c("darkgrey", "black"), 22 )) +

        # Show significance threshold
        geom_hline(yintercept=yintercept, linetype="dashed", color = "red") +

        # Custom axes:
        scale_x_continuous(label=axisdf$CHR, breaks=axisdf$center, limits=c(limit.left,limit.right),
                           expand=c(0.01,0.01)) +
        scale_y_continuous(trans=ytrans, expand = c(0,0), limits=c(0,y.max)) + 
        xlab("Chromosome") + ylab("Importance") +
    
        # Custom the theme:
        theme_bw() +
        theme(legend.position="none", panel.border = element_blank(), 
              panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
              panel.grid.major = element_line(size = 0.2, colour = "darkgray"),
              panel.grid.minor = element_line(size = 0.1, colour = "lightgray"),
              text = element_text(size=font.size)
              )

    return(p.manhattan)
}

plot_manhattan <- function(gwas, ytrans="log10") {
    
    don <- pre_process(gwas)
    
    # Compute plot limits
    limit.left <- min(don$BPcum)
    limit.right <- max(don$BPcum)
    
    # Find centers of each chromosome
    axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

    # Make plot
    p.manhattan <- manhattan_simple(don, axisdf, limit.left, limit.right, ytrans=ytrans)
    
    return(p.manhattan)
}

plot_manhattan_knockoffs <- function(LMM, Knockoffs, ytrans="identity", chr=NULL) {
    
    if(!is.null(chr)) {
        if(nrow(LMM)>0) LMM <- LMM %>% filter(CHR==chr)
        if(nrow(Knockoffs)>0) Knockoffs <- Knockoffs %>% filter(CHR==chr)
    }

    if(nrow(LMM)>0) {
        # Create chromosome blocks
        don <- pre_process(LMM)
        # Compute plot limits
        limit.left <- min(don$BPcum)
        limit.right <- max(don$BPcum)
        # Find centers of each chromosome
        axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
        # Make LMM Manhattan plot
        p.manhattan <- manhattan_simple(don, axisdf, limit.left, limit.right, ytrans=ytrans)
        p.manhattan <- p.manhattan +
            ylab(TeX("$-\\log_{10}(p)$")) +
            theme(axis.title.x=element_blank())        
    } else {
        p.manhattan <- ggplot(tibble())
    }
    
    # Transform knockoffs results into pseudo p-values
    W.thresh <- knockoff.threshold(Knockoffs$W)
    Knockoffs <- Knockoffs %>% 
        filter(W>0) %>% 
        mutate(SNP=SNP.lead, BP=BP.lead) %>%
        select(CHR, SNP, BP, W)
    
    # Create chromosome blocks
    if(nrow(LMM)>0) {
        Knockoffs.don <- Knockoffs %>%
            select(CHR, BP, W) %>%
            left_join(don, by = c("CHR", "BP"))
    } else {
        Knockoffs.don <- pre_process(Knockoffs)
        don <- Knockoffs.don
        # Compute plot limits
        limit.left <- min(don$BPcum)
        limit.right <- max(don$BPcum)
        # Find centers of each chromosome
        axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
    }
       
    # Make Knockoffs Manhattan plot
    p.knockoffs <- manhattan_knock(Knockoffs.don, axisdf, limit.left, limit.right, ytrans=ytrans,
                                   yintercept=W.thresh)
    p.knockoffs <- p.knockoffs +
        scale_y_continuous(name="Test statistics", labels=function(x) round(x,3))
    
    # Modify plots if we are focusing on a single chromosome
    if(!is.null(chr)) {
      p.knockoffs <- p.knockoffs +
        xlab(sprintf("Mb (Chromosome %d)", chr)) +
        scale_x_continuous(expand=c(0.01,0.01)) +
        scale_color_manual(values = "black")
      p.manhattan <- p.manhattan +
        xlab(sprintf("Mb (Chromosome %d)", chr)) +
        scale_x_continuous(expand=c(0.01,0.01)) +
        scale_color_manual(values = "black")
    }
    
    plot_grid(p.manhattan, p.knockoffs, ncol=1, align="v", axis="tblr") 
}
