# Calculates the effective population size (Ne) from allele frequency trajectories
# v.1.01 - Oct 2025

# ----LOAD LIBRARIES----
# Load library
library(poolSeq)
library(data.table)
library(foreach)
library(stringi)
library(matrixStats)
library(Rcpp)

# ----FUNCTIONS & DATASET----
# Load in sync file into variable 'Sync'
Sync <- read.sync(file = '/home/baron/Documents/PhD/Data/pop_size_analysis/PoolSeq/cages/cage_10.sync', gen = c(2,4,8,12,20,28,36,44,56), repl = 1)

# Load parameters from Sync into function
estimateNe_from_trajectories_A <- function(sync,
                                           gens = c(2,4,8,12,20,28,36,44,56),
                                           repl = 1,
                                           gen0 = gens[1],          # baseline generation (2)
                                           minCov = 20,             # Coverage threshold, change when required
                                           minSites = 1000,
                                           nboot = 200,             # bootstrap replicates for CI
                                           bootSampleSize = 20000,  # sample size per bootstrap (<= nSites)
                                           seed = 1,
                                           verbose = TRUE) {
  set.seed(seed)
  stopifnot(gen0 %in% gens)
  other_gens <- gens[gens != gen0]
  
  # pull p0 and cov0
  p0_all   <- af(sync, gen = gen0, repl = repl)
  cov0_all <- coverage(sync, gen = gen0, repl = repl)
  
  # for each later generation compute corrected F
  t_vals <- numeric(0)
  F_vals <- numeric(0)
  nSites_vec <- integer(0)
  
  # Store per-locus data keyed by names for bootstrap convenience
  # Vectors per later gen of the locus keys kept
  kept_keys_list <- list()
  
  for (g in other_gens) {
    p_g_all   <- af(sync, gen = g, repl = repl)
    cov_g_all <- coverage(sync, gen = g, repl = repl)
    
    # align keys present in both
    common <- intersect(names(p0_all), names(p_g_all))
    if (length(common) == 0) {
      if (verbose) message("No common loci for gen ", gen0, " and gen ", g)
      next
    }
    
    p0 <- as.numeric(p0_all[common])
    pg <- as.numeric(p_g_all[common])
    cov0 <- as.numeric(cov0_all[common])
    covg <- as.numeric(cov_g_all[common])
    
    # filters: no NA, cov >= minCov at both times, remove fixed in both!!!
    keep <- !is.na(p0) & !is.na(pg) & !is.na(cov0) & !is.na(covg) &
      cov0 >= minCov & covg >= minCov &
      !((p0 == 0 & pg == 0) | (p0 == 1 & pg == 1))
    
    n_sites <- sum(keep)
    if (verbose) message(sprintf("Interval %d->%d: %d loci kept", gen0, g, n_sites))
    if (n_sites < minSites) {
      if (verbose) message("Too few sites; skipping interval ", gen0, "->", g)
      next
    }
    
    p0k <- p0[keep]; pgk <- pg[keep]; cov0k <- cov0[keep]; covgk <- covg[keep]
    keys_kept <- common[keep]
    kept_keys_list[[as.character(g)]] <- keys_kept
    
    # Account for sampling variance
    # sampling variance per locus: approx p0(1-p0)/cov0 + p_g(1-p_g)/covg
    sampVar <- p0k*(1-p0k)/cov0k + pgk*(1-pgk)/covgk
    
    # corrected squared change
    sqchg_corr <- (pgk - p0k)^2 - sampVar
    
    # replace small negative corrected values by 0
    sqchg_corr[sqchg_corr < 0] <- 0
    
    # standardized F per locus: (corrected sq change) / [p0(1-p0)]
    denom <- p0k*(1-p0k)
    # avoid dividing by extremely small denom: drop loci with p0*(1-p0) < tiny
    tiny <- 1e-8
    good <- denom > tiny
    if (sum(good) < minSites) {
      if (verbose) message("After removing tiny-p loci, too few sites for interval ", gen0, "->", g)
      next
    }
    
    F_hat <- mean(sqchg_corr[good] / denom[good], na.rm = TRUE)
    
    t_vals <- c(t_vals, g - gen0)
    F_vals <- c(F_vals, F_hat)
    nSites_vec <- c(nSites_vec, sum(good))
  }
  
  if (length(t_vals) < 2) stop("Not enough intervals with data to fit regression.")
  
  # Fit linear regression through origin: F = s * t  => s = coef(lm(F ~ 0 + t))
  lm_fit <- lm(F_vals ~ 0 + t_vals, weights = nSites_vec)  # weight by number of loci per t
  s_hat <- as.numeric(coef(lm_fit))
  Ne_hat <- 1 / (2 * s_hat)
  
  # bootstrap to get CIs: resample loci (by index) per interval and refit regression
  boot_ne <- rep(NA_real_, nboot)
  all_keys_union <- unique(unlist(kept_keys_list))
  # sample indices within each interval independently (with replacement) and compute F, then regress
  for (b in seq_len(nboot)) {
    F_boot <- numeric(0)
    t_boot <- numeric(0)
    ns_boot <- numeric(0)
    for (j in seq_along(other_gens)) {
      g <- other_gens[j]
      keys <- kept_keys_list[[as.character(g)]]
      if (is.null(keys)) next
      # sample indices
      n_available <- length(keys)
      samp_n <- min(bootSampleSize, n_available)
      idx <- sample(seq_len(n_available), size = samp_n, replace = FALSE) #  SAMPLE WITHOUT REPLACEMENT - change to true if required
      # get vectors for sampled loci
      p0_all_sub   <- as.numeric(p0_all[keys])[idx]
      pg_all_sub   <- as.numeric(af(sync, gen = g, repl = repl)[keys])[idx]
      cov0_sub     <- as.numeric(cov0_all[keys])[idx]
      covg_sub     <- as.numeric(coverage(sync, gen = g, repl = repl)[keys])[idx]
      
      sampVar_sub <- p0_all_sub*(1-p0_all_sub)/cov0_sub + pg_all_sub*(1-pg_all_sub)/covg_sub
      sqcorr_sub <- (pg_all_sub - p0_all_sub)^2 - sampVar_sub
      sqcorr_sub[sqcorr_sub < 0] <- 0
      denom_sub <- p0_all_sub*(1-p0_all_sub)
      good_sub <- denom_sub > 1e-8
      if (sum(good_sub) < 10) {
        F_boot <- c(F_boot, NA); t_boot <- c(t_boot, g - gen0); ns_boot <- c(ns_boot, sum(good_sub))
      } else {
        F_boot <- c(F_boot, mean(sqcorr_sub[good_sub] / denom_sub[good_sub], na.rm = TRUE))
        t_boot <- c(t_boot, g - gen0)
        ns_boot <- c(ns_boot, sum(good_sub))
      }
    }
    # only fit if >1 point
    if (sum(!is.na(F_boot)) >= 2) {
      fitb <- tryCatch(lm(F_boot ~ 0 + t_boot, weights = ns_boot), error = function(e) NULL)
      if (!is.null(fitb)) {
        s_b <- as.numeric(coef(fitb))
        if (is.finite(s_b) && s_b > 0) boot_ne[b] <- 1 / (2 * s_b)
      }
    }
  }
  
  boot_ne <- boot_ne[!is.na(boot_ne)]
  out <- list(
    gen0 = gen0,
    t = t_vals,
    F = F_vals,
    nSites = nSites_vec,
    lm = lm_fit,
    Ne = Ne_hat,
    Ne_boot = boot_ne,
    Ne_CI = if (length(boot_ne) > 0) quantile(boot_ne, probs = c(0.025, 0.5, 0.975)) else NA
  )
  class(out) <- "NeTrajEstimate"
  return(out)
}

# ----  CREATE WRAPPER -----
summary.NeTrajEstimate <- function(object, ...) {
  if (!inherits(object, "NeTrajEstimate"))
    stop("Object must be of class 'NeTrajEstimate'")
  
  cat("Effective population size estimate from allele frequency trajectories\n")
  cat("---------------------------------------------------------------\n")
  cat("Baseline generation (gen0):", object$gen0, "\n")
  cat("Number of intervals used:   ", length(object$t), "\n")
  cat("Point estimate Ne:          ", round(object$Ne, 2), "\n")
  
  if (!is.null(object$Ne_CI) && !all(is.na(object$Ne_CI))) {
    cat("95% bootstrap CI:           ",
        paste(round(object$Ne_CI, 2), collapse = " – "), "\n")
  } else {
    cat("95% bootstrap CI:           Not available\n")
  }
  
  # R² from regression
  rsq <- summary(object$lm)$r.squared
  cat("Regression R²:              ", round(rsq, 3), "\n")
  
  cat("\nIntervals (generations since gen0):\n")
  df <- data.frame(
    t = object$t,
    F_hat = signif(object$F, 4),
    nSites = object$nSites
  )
  print(df, row.names = FALSE)
  
  invisible(df)
}

# ----OUTPUT RESULTS----
# Create table and print
resA <- estimateNe_from_trajectories_A(Sync)
summary(resA)
