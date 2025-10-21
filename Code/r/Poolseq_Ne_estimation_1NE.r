# Calculates the effective population size (Ne) from allele frequency trajectories - single Ne for 1 population
# v.1.02 - Oct 2025

# ---- LOAD LIBRARIES ----
# Load library
library(poolSeq)
library(data.table)
library(foreach)
library(stringi)
library(matrixStats)
library(Rcpp)

# ---- FUNCTIONS & DATASET ----
# Load in sync file into variable 'Sync'
Sync <- read.sync(file = '/home/baron/Documents/PhD/Data/pop_size_analysis/PoolSeq/cages/cage_3.sync', gen = c(2,4,8,12,20,28,36,44,56), repl = 1)

# Load parameters from Sync into function
# ---- FUNCTION TO ESTIMATE Ne PER INTERVAL ----
# ---- FUNCTION TO ESTIMATE Ne PER INTERVAL WITH BOOTSTRAP ----
# ----COMBINED FUNCTION: PER-INTERVAL + OVERALL Ne----
estimateNe_combined <- function(sync,
                                gens = c(2,4,8,12,20,28,36,44,56),
                                repl = 1,
                                minCov = 20,
                                minSites = 1000,
                                nboot = 200,
                                bootSampleSize = 20000,
                                seed = 1,
                                verbose = TRUE) {
  set.seed(seed)
  n_intervals <- length(gens) - 1
  
  # Store per-interval results
  res_list <- vector("list", n_intervals)
  
  # For regression (overall Ne)
  t_vals <- numeric(0)
  F_vals <- numeric(0)
  nSites_vec <- numeric(0)
  
  for (i in 1:n_intervals) {
    gen0 <- gens[i]
    gen1 <- gens[i+1]
    
    p0_all   <- af(sync, gen = gen0, repl = repl)
    cov0_all <- coverage(sync, gen = gen0, repl = repl)
    p1_all   <- af(sync, gen = gen1, repl = repl)
    cov1_all <- coverage(sync, gen = gen1, repl = repl)
    
    common <- intersect(names(p0_all), names(p1_all))
    if (length(common) == 0) next
    
    p0 <- as.numeric(p0_all[common])
    p1 <- as.numeric(p1_all[common])
    cov0 <- as.numeric(cov0_all[common])
    cov1 <- as.numeric(cov1_all[common])
    
    keep <- !is.na(p0) & !is.na(p1) & !is.na(cov0) & !is.na(cov1) &
      cov0 >= minCov & cov1 >= minCov &
      !((p0 == 0 & p1 == 0) | (p0 == 1 & p1 == 1))
    
    n_sites <- sum(keep)
    if (n_sites < minSites) {
      if (verbose) message("Too few sites for interval ", gen0, "->", gen1)
      next
    }
    
    p0k <- p0[keep]; p1k <- p1[keep]; cov0k <- cov0[keep]; cov1k <- cov1[keep]
    
    # Corrected F
    sampVar <- p0k*(1-p0k)/cov0k + p1k*(1-p1k)/cov1k
    sqchg_corr <- (p1k - p0k)^2 - sampVar
    sqchg_corr[sqchg_corr < 0] <- 0
    denom <- p0k*(1-p0k)
    tiny <- 1e-8
    good <- denom > tiny
    if (sum(good) < minSites) next
    
    F_hat <- mean(sqchg_corr[good] / denom[good], na.rm = TRUE)
    Ne_hat <- ifelse(F_hat > 0, 1/(2*F_hat), NA)
    
    # Bootstrap for CI per interval
    boot_ne <- numeric(nboot)
    keys <- which(good)
    n_keys <- length(keys)
    for (b in 1:nboot) {
      idx <- sample(keys, size = min(bootSampleSize, n_keys), replace = TRUE)
      p0b <- p0k[idx]; p1b <- p1k[idx]; cov0b <- cov0k[idx]; cov1b <- cov1k[idx]
      sampVar_b <- p0b*(1-p0b)/cov0b + p1b*(1-p1b)/cov1b
      sqb <- (p1b - p0b)^2 - sampVar_b
      sqb[sqb < 0] <- 0
      denom_b <- p0b*(1-p0b)
      good_b <- denom_b > tiny
      if (sum(good_b) < 10) {
        boot_ne[b] <- NA
      } else {
        F_b <- mean(sqb[good_b] / denom_b[good_b], na.rm = TRUE)
        boot_ne[b] <- ifelse(F_b > 0, 1/(2*F_b), NA)
      }
    }
    
    boot_ne <- boot_ne[!is.na(boot_ne)]
    Ne_CI <- if(length(boot_ne) > 0) quantile(boot_ne, probs = c(0.025, 0.5, 0.975)) else c(NA, NA, NA)
    
    res_list[[i]] <- data.frame(
      gen0 = gen0,
      gen1 = gen1,
      Ne = Ne_hat,
      Ne_lower = Ne_CI[1],
      Ne_median = Ne_CI[2],
      Ne_upper = Ne_CI[3],
      nSites = sum(good)
    )
    
    # Add to regression vectors
    t_vals <- c(t_vals, gen1 - gen0)
    F_vals <- c(F_vals, F_hat)
    nSites_vec <- c(nSites_vec, sum(good))
  }
  
  # Overall Ne by regression
  lm_fit <- lm(F_vals ~ 0 + t_vals, weights = nSites_vec)
  s_hat <- as.numeric(coef(lm_fit))
  Ne_overall <- 1 / (2 * s_hat)
  
  # Bootstrap overall Ne
  boot_ne_overall <- numeric(nboot)
  for (b in 1:nboot) {
    F_boot <- numeric(0); t_boot <- numeric(0); ns_boot <- numeric(0)
    for (i in 1:length(res_list)) {
      row <- res_list[[i]]
      if (is.null(row)) next
      boot_F <- 1 / (2 * row$Ne_median)  # use median Ne to compute F for resample
      F_boot <- c(F_boot, boot_F)
      t_boot <- c(t_boot, row$gen1 - row$gen0)
      ns_boot <- c(ns_boot, row$nSites)
    }
    if (length(F_boot) >= 2) {
      fitb <- tryCatch(lm(F_boot ~ 0 + t_boot, weights = ns_boot), error = function(e) NULL)
      if (!is.null(fitb)) {
        s_b <- as.numeric(coef(fitb))
        if (is.finite(s_b) && s_b > 0) boot_ne_overall[b] <- 1 / (2 * s_b)
      }
    }
  }
  boot_ne_overall <- boot_ne_overall[!is.na(boot_ne_overall)]
  Ne_overall_CI <- if(length(boot_ne_overall) > 0) quantile(boot_ne_overall, probs = c(0.025, 0.5, 0.975)) else c(NA, NA, NA)
  
  # Combine interval results
  res_intervals <- do.call(rbind, res_list)
  
  return(list(
    intervals = res_intervals,
    Ne_overall = Ne_overall,
    Ne_overall_CI = Ne_overall_CI,
    lm_overall = lm_fit
  ))
}

# ---- PRINT OUTPUT ----
res_combined <- estimateNe_combined(Sync)

# Interval-specific Ne
print(res_combined$intervals)

# Overall regression-based Ne
res_combined$Ne_overall
res_combined$Ne_overall_CI
