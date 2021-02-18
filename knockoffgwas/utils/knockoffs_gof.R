#!/usr/bin/env Rscript

plot.knockoff.diagnostics <- function(stats.basename, key.file, groups.file, out.basename) {

    # Load variable grouping
    Variants <- read_delim(groups.file, delim=" ", progress = FALSE, col_types = cols()) %>%
        mutate(Group = as.integer(Group))

    # Load knockoff keys
    Key <- read_delim(key.file, delim=" ", col_types=cols())

    # Load variant frequency table
    frq.file <- sprintf("%s.frq", stats.basename)
    Frq <- read_table(frq.file, col_types=cols())
    # Compute diagnostics
    maf <- Frq$MAF
    maf.x  <- maf[seq(1,length(maf), by=2)]
    maf.xk <- maf[seq(2,length(maf), by=2)]
    Diagnostics <- tibble(index=seq(length(maf.x)), x=maf.x, xk=maf.xk) %>%
        mutate(error=(x-xk)^2)

    # Plot frequency diagnostics
    p.frq <- Diagnostics %>%
        ggplot(aes(x=x, y=xk)) +
        geom_point(alpha=0.2) +
        geom_abline(intercept = 0, slope = 1, color="red", linetype=2) +
        xlab(TeX("MAF ($X$)")) + ylab(TeX("MAF ($\\tilde{X}$)")) +
        theme_bw()

    # Load LD table
    ld.file <- sprintf("%s.ld", stats.basename)
    LD <- suppressWarnings(read_table2(ld.file, col_types=cols())) %>%
        mutate(CHR=CHR_A) %>% select(-CHR_A, -CHR_B)

    # Add grouping information
    LD <- LD %>%
        mutate(SNP=SNP_A) %>%
        mutate(SNP = gsub(".A", "", SNP), SNP = gsub(".B", "", SNP)) %>%
        inner_join(Variants %>% select(CHR, SNP, Group), by = c("CHR", "SNP")) %>%
        mutate(Group_A=Group) %>% select(-SNP, -Group) %>%
        mutate(SNP=SNP_B) %>%
        mutate(SNP = gsub(".A", "", SNP), SNP = gsub(".B", "", SNP)) %>%
        inner_join(Variants %>% select(CHR, SNP, Group), by = c("CHR", "SNP")) %>%
        mutate(Group_B=Group) %>% select(-SNP, -Group)

    # Add knockoff key information
    LD <- LD %>%
        mutate(SNP=SNP_A) %>%
        inner_join(Key %>% select(CHR, SNP, Knockoff), by = c("CHR", "SNP")) %>%
        mutate(Knockoff_A=Knockoff) %>% select(-SNP, -Knockoff) %>%
        mutate(SNP=SNP_B) %>%
        inner_join(Key %>% select(CHR, SNP, Knockoff), by = c("CHR", "SNP")) %>%
        mutate(Knockoff_B=Knockoff) %>% select(-SNP, -Knockoff)

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

    # Plot originality
    LD.cross <- inner_join(LD.XX, LD.XkXk, by = c("Group_A", "Group_B", "SNP_A", "SNP_B"))
    p.orig <- LD.cross %>%
        mutate(Distance = as.factor(abs(Group_A-Group_B))) %>%
        ggplot(aes(x=abs(R.XX), y=abs(R.XkXk))) +
        geom_abline(color="red") +
        geom_point(alpha=0.1) +
        xlim(0,1) + ylim(0,1) +
        xlab(TeX("|corr($X_{j},X_{k}$)|")) + ylab(TeX("|corr($\\tilde{X}_{j},\\tilde{X}_{k}$)|")) +
        theme_bw()

    # Plot exchangeability
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

    # Plot histogram of self-correlations
    p.self <- LD %>% filter(BP_A==BP_B) %>%
        ggplot(aes(x=R2)) +
        geom_histogram(bins=30) +
        xlab(TeX("|corr($X_{j},\\tilde{X}_{j}$|)")) +
        theme_bw()

    # Combine plots
    p.combined <- gridExtra::grid.arrange(p.frq, p.orig, p.exch, p.self, nrow = 2)
    return(p.combined)
}
