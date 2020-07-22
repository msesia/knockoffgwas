suppressMessages(library(tibbletime))
suppressMessages(library(scales))

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


knockoff.filter <- function(Stats, fdr=0.1, offset=1) {
    W.thres <- knockoff.threshold(Stats$W, fdr=fdr, offset=offset)
    Selected <- Stats %>% filter(W >= W.thres)
    return(Selected)
}

local_fdp <- function(Discoveries, Stats) {

    # Define function to estimate the local FDP
    local.fdp.helper <- function(W, offset=0) {
        t <- 0
        fdp.hat <- (offset + sum(W <= t)) / max(1, sum(W >= t))
        fdp.hat <- min(1,fdp.hat)
        return(fdp.hat)
    }
    if(nrow(Discoveries)>500) {
        window <- 100
    } else if(nrow(Discoveries)>=100) {
        window <- 10
    } else if(nrow(Discoveries)>=20) {
        window <- 5
    } else {
        cat("The number of discoveries is too small to estimate the local FDP.\n")
        return(Discoveries %>% mutate(FDP.local=NA))
    }
    estimate.local.fdp <- rollify(local.fdp.helper, window=window, na_value=0)

    # Compute local FDP at each resolution
    res.levels <- unique(Discoveries$Resolution)
    Discoveries <- lapply(res.levels, function(resolution) {
        # Extract values at this resolution
        Discoveries.res <- Discoveries %>% filter(Resolution==resolution)
        Stats.res <- Stats %>% filter(Resolution==resolution)
        # Sort the statistics
        Stats.res <- Stats.res %>%
            arrange(desc(abs(W)))
        # Estimate the FDP
        Stats.fdp <- Stats.res %>%
            mutate(FDP.local=estimate.local.fdp(W)) %>%
            select(Resolution, CHR, SNP.lead, FDP.local)
        # Cross-reference with selections  
        Discoveries.res <- Discoveries.res %>%
            left_join(Stats.fdp, by = c("CHR", "SNP.lead", "Resolution"))
        return(Discoveries.res)
    })
    Discoveries <- do.call("rbind", Discoveries)

    return(Discoveries)
}



