# Robust genome-wide pairwise Ne estimation using poolSeq
estimateNe_intervals_diag <- function(sync,
                                      gens = c(2,4,8,12,20,28,36,44,56),
                                      repl = 1,
                                      poolSize = 96,      # 48 diploid individuals -> 96 chromosomes
                                      minCov = 10,        # min coverage required at both timepoints
                                      minSites = 100,     # require at least this many SNPs to run estimateNe
                                      maxSites = NULL,    # if set, randomly sample up to this many sites (to limit mem)
                                      bootstrap = FALSE,
                                      B = 200,            # bootstrap replicates (only if bootstrap = TRUE)
                                      bootSampleSize = NULL, # number of sites per bootstrap (NULL => use all kept sites)
                                      verbose = TRUE) {
  
  stopifnot(is.numeric(poolSize), length(poolSize) == 1 || length(poolSize) == 2)
  
  # prepare output table
  out <- data.frame(gen0 = gens[-length(gens)],
                    gen1 = gens[-1],
                    Ne = NA_real_,
                    nSites = NA_integer_,
                    boot_median = NA_real_,
                    boot_lower = NA_real_,
                    boot_upper = NA_real_,
                    stringsAsFactors = FALSE)
  
  for (i in seq_len(nrow(out))) {
    g0 <- out$gen0[i]; g1 <- out$gen1[i]
    if (verbose) message(sprintf("\n---- Interval %d -> %d ----", g0, g1))
    
    # pull vectors from Sync
    p0_v   <- af(sync, gen = g0, repl = repl)
    pt_v   <- af(sync, gen = g1, repl = repl)
    cov0_v <- coverage(sync, gen = g0, repl = repl)
    covt_v <- coverage(sync, gen = g1, repl = repl)
    
    # align by common names (chr.pos keys)
    common <- intersect(names(p0_v), names(pt_v))
    if (length(common) == 0) {
      warning(sprintf("No common loci between gen %d and %d", g0, g1))
      next
    }
    p0_v   <- p0_v[common]; pt_v   <- pt_v[common]
    cov0_v <- cov0_v[common]; covt_v <- covt_v[common]
    
    # convert to numeric vectors (unnamed)
    p0n   <- unname(as.numeric(p0_v))
    ptn   <- unname(as.numeric(pt_v))
    cov0n <- unname(as.numeric(cov0_v))
    covtn <- unname(as.numeric(covt_v))
    
    # filtering:
    #  - remove NA entries
    #  - require coverage >= minCov at both timepoints
    #  - remove sites that are fixed at 0 or 1 at both timepoints (no information)
    keep <- !is.na(p0n) & !is.na(ptn) & !is.na(cov0n) & !is.na(covtn) &
      (cov0n >= minCov) & (covtn >= minCov) &
      !( (p0n == 0 & ptn == 0) | (p0n == 1 & ptn == 1) )
    
    n_total <- length(p0n)
    n_keep  <- sum(keep)
    
    if (verbose) {
      message(sprintf("Total loci available: %d", n_total))
      message(sprintf("Loci kept after filters (minCov=%d, exclude fixed): %d", minCov, n_keep))
    }
    
    out$nSites[i] <- n_keep
    
    if (n_keep < minSites) {
      warning(sprintf("Too few sites (%d) remain after filtering for interval %d->%d; returning NA", n_keep, g0, g1))
      next
    }
    
    # optionally subsample sites to limit memory/cpu
    keep_idx <- which(keep)
    if (!is.null(maxSites) && length(keep_idx) > maxSites) {
      set.seed(1) # reproducible sampling
      keep_idx <- sample(keep_idx, maxSites)
      if (verbose) message(sprintf("Subsampled to %d sites (maxSites=%d)", length(keep_idx), maxSites))
    }
    
    # final vectors for estimateNe
    p0_final   <- p0n[keep_idx]
    pt_final   <- ptn[keep_idx]
    cov0_final <- cov0n[keep_idx]
    covt_final <- covtn[keep_idx]
    
    # try running estimateNe
    est <- tryCatch({
      estimateNe(p0 = p0_final, pt = pt_final,
                 cov0 = cov0_final, covt = covt_final,
                 t = g1 - g0,
                 Ncensus = NA,
                 poolSize = rep(as.numeric(poolSize), 2))  # same pool size both times
    }, error = function(e) {
      message("estimateNe error: ", conditionMessage(e))
      return(NA_real_)
    })
    
    out$Ne[i] <- as.numeric(est)
    
    if (verbose) message(sprintf("Raw estimateNe for %d->%d: %s", g0, g1, ifelse(is.na(out$Ne[i]), "NA", format(out$Ne[i], digits=6))))
    
    # optional bootstrap for CI
    if (bootstrap && !is.na(out$Ne[i])) {
      if (verbose) message("Starting bootstrap ...")
      n_sites <- length(p0_final)
      if (is.null(bootSampleSize)) bootSampleSize <- n_sites
      
      boot_res <- numeric(B)
      for (b in seq_len(B)) {
        samp_idx <- sample(seq_len(n_sites), bootSampleSize, replace = TRUE)
        boot_res[b] <- tryCatch({
          estimateNe(p0 = p0_final[samp_idx],
                     pt = pt_final[samp_idx],
                     cov0 = cov0_final[samp_idx],
                     covt = covt_final[samp_idx],
                     t = g1 - g0,
                     Ncensus = NA,
                     poolSize = rep(as.numeric(poolSize), 2))
        }, error = function(e) NA_real_)
      }
      boot_res <- boot_res[!is.na(boot_res) & is.finite(boot_res)]
      if (length(boot_res) == 0) {
        warning("All bootstrap estimates returned NA/NaN.")
      } else {
        out$boot_median[i] <- median(boot_res)
        out$boot_lower[i]  <- quantile(boot_res, probs = 0.025)
        out$boot_upper[i]  <- quantile(boot_res, probs = 0.975)
        if (verbose) message(sprintf("Bootstrap median (95%% CI): %g (%g - %g)",
                                     out$boot_median[i], out$boot_lower[i], out$boot_upper[i]))
      }
    }
  } # end intervals loop
  
  return(out)
}
