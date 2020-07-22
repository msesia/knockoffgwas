are.disjoint <- function(a1, b1, a2, b2, gap=1e5) {
    # Check whether two intervals are disjoint
    if(a1>a2) {
        a.tmp <- a1
        b.tmp <- b1
        a1 <- a2
        b1 <- b2
        a2 <- a.tmp
        b2 <- b.tmp
    }
    if(b1<a2-gap) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

clump.segments <- function(bp.min, bp.max, gap=1e5) {
    # Check whether any two segments are less than 100kB apart
    # If they are, assign them to the same clump
    if(length(bp.min)==1) {
        return(c(1))
    }
    
    # Compute adjacency matrix
    A <- sapply(1:length(bp.min), function(i) { 
        sapply(1:length(bp.min), function(j) {
            1-are.disjoint(bp.min[i],bp.max[i],bp.min[j],bp.max[j],gap=gap)
        })
    })
        
    # Find connected components
    A.graph <- igraph::graph.adjacency(A)
    A.clu <- igraph::components(A.graph)
    
    # Return group membership
    return(A.clu$membership)
}

my.which.min <- function(x) {
    if(all(is.na(x))) {
        return(NA)
    } else {
        return(which.min(x))
    }
}

consolidate_clumps <- function(Clumped, gap=1e5) {   
    # Store column names
    original.colnames <- colnames(Clumped)
    # Store groups
    original.groups <- group_vars(Clumped)

    if(! "P" %in% original.colnames) {
        Clumped <- Clumped %>%
            mutate(P=10^(-Importance)) %>%
            select(-Importance)
    }
    
    # Form consolidated clumps for each group
    Segments <- Clumped %>% 
        group_by(.dots = c(original.groups, "CHR")) %>% 
        mutate(Clump=clump.segments(BP.min, BP.max, gap=gap))
    
    # Update lead SNPs for each group
    if(! "Causal" %in% original.colnames) {
        Consolidated <- Segments %>%
            group_by(.dots = c(original.groups, "CHR", "Clump")) %>%
            summarise(Lead=my.which.min(P), P=min(P),
                      SNP.lead=SNP.lead[Lead], BP.lead=BP.lead[Lead],
                      BP.min=min(BP.min), BP.max=max(BP.max), BP.width=BP.max-BP.min,
                      Size=sum(Size)) %>%
            ungroup()
    } else {
        Consolidated <- Segments %>%
            group_by(.dots = c(original.groups, "CHR", "Clump")) %>%
            summarise(Lead=my.which.min(P), P=min(P),
                      SNP.lead=SNP.lead[Lead], BP.lead=BP.lead[Lead],
                      BP.min=min(BP.min), BP.max=max(BP.max), BP.width=BP.max-BP.min,
                      Causal=any(Causal),
                      Size=sum(Size)) %>%
            ungroup()
    }

    if(! "P" %in% original.colnames) {
        Consolidated <- Consolidated %>%
            mutate(Importance=-log10(P)) %>%
            select(-P)
    }

    # Return results
    new.colnames <- intersect(original.colnames, colnames(Consolidated))
    Consolidated %>% 
        select(new.colnames) %>%
        group_by(.dots = original.groups)
} 

consolidate_clumps_variants <- function(Clumped, threshold=0) {
    # Store column names
    original.colnames <- colnames(Clumped)
    # Store groups
    original.groups <- group_vars(Clumped)
    
    # Form consolidated clumps for each group
    Segments <- Clumped %>% 
        group_by(.dots = c(original.groups, "CHR", "SNP.lead", "BP.lead")) %>% 
        summarise(BP.min=min(BP[Importance>=threshold]), BP.max=max(BP[Importance>=threshold])) %>% 
        ungroup() %>% 
        group_by(.dots = c(original.groups, "CHR")) %>% 
        mutate(Clump=clump.segments(BP.min, BP.max))
    
    # Update lead SNPs for each group
    Consolidated <- Clumped %>%
        inner_join(Segments, by = c("CHR", "SNP.lead", "BP.lead", "Phenotype", "Method", "Resolution")) %>%
        group_by(.dots = c(original.groups, "CHR", "Clump")) %>% 
        mutate(Lead=which.max(Importance), SNP.lead=SNP[Lead], BP.lead=BP[Lead],
               BP.min=min(BP), BP.max=max(BP)) %>%
        ungroup()
    
    # Return results
    Consolidated %>% 
        select(original.colnames) %>%
        group_by(.dots = original.groups)
} 
