# Robust genome-wide pairwise Ne estimation using poolSeq
# Calculates the effective population size (Ne) between time point intervals using temporal allele frequency
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
Sync <- read.sync(file = '/home/baron/Documents/PhD/Data/pop_size_analysis/PoolSeq/cages/cage_2.sync', gen = c(2,4,8,12,20,28,36,44,56), repl = 1)

# ---- NE FROM BASELINE GEN0 USING TRAJECTORY METHOD (PER GENERATION) ----
estimateNe_from_baseline_F <- function(sync,
                                       gens = c(2,4,8,12,20,28,36,44,56),
                                       gen0 = 2,
                                       repl = 1,
                                       minCov = 20,
                                       minSites = 1000,
                                       nboot = 200,
                                       bootSampleSize = 20000,
                                       seed = 1,
                                       verbose = TRUE,
                                       bootstrap = TRUE) {
  set.seed(seed)
  
  # If bootstrap = FALSE, ignore nboot
  if (!bootstrap) nboot <- 0
  
  # Only consider generations after gen0
  later_gens <- gens[gens > gen0]
  n_intervals <- length(later_gens)
  
  # Prepare output table
  out <- data.frame(gen0 = rep(gen0, n_intervals),
                    gen1 = later_gens,
                    Ne = NA_real_,
                    nSites = NA_integer_,
                    Ne_lower = NA_real_,
                    Ne_median = NA_real_,
                    Ne_upper = NA_real_,
                    stringsAsFactors = FALSE)
  
  # Baseline allele frequencies and coverage
  p0_all <- af(sync, gen = gen0, repl = repl)
  cov0_all <- coverage(sync, gen = gen0, repl = repl)
  
  for (i in seq_len(n_intervals)) {
    g1 <- later_gens[i]
    p1_all <- af(sync, gen = g1, repl = repl)
    cov1_all <- coverage(sync, gen = g1, repl = repl)
    
    # Keep common loci
    common <- intersect(names(p0_all), names(p1_all))
    if (length(common) == 0) next
    
    p0 <- as.numeric(p0_all[common])
    p1 <- as.numeric(p1_all[common])
    cov0 <- as.numeric(cov0_all[common])
    cov1 <- as.numeric(cov1_all[common])
    
    # Filtering
    keep <- !is.na(p0) & !is.na(p1) & !is.na(cov0) & !is.na(cov1) &
      cov0 >= minCov & cov1 >= minCov &
      !((p0 == 0 & p1 == 0) | (p0 == 1 & p1 == 1))
    
    n_sites <- sum(keep)
    if (verbose) message(sprintf("Interval %d -> %d: loci kept = %d", gen0, g1, n_sites))
    out$nSites[i] <- n_sites
    if (n_sites < minSites) next
    
    p0k <- p0[keep]; p1k <- p1[keep]; cov0k <- cov0[keep]; cov1k <- cov1[keep]
    
    # Corrected squared allele frequency change (F)
    sampVar <- p0k*(1-p0k)/cov0k + p1k*(1-p1k)/cov1k
    sqchg_corr <- (p1k - p0k)^2 - sampVar
    sqchg_corr[sqchg_corr < 0] <- 0
    denom <- p0k*(1-p0k)
    tiny <- 1e-8
    good <- denom > tiny
    if (sum(good) < minSites) next
    
    # Compute F per generation
    F_hat <- mean(sqchg_corr[good] / denom[good], na.rm = TRUE)
    F_per_gen <- F_hat / (g1 - gen0)
    out$Ne[i] <- ifelse(F_per_gen > 0, 1/(2*F_per_gen), NA)
    
    # Bootstrap
    if (nboot > 0) {
      boot_ne <- numeric(nboot)
      keys <- which(good)
      n_keys <- length(keys)
      for (b in seq_len(nboot)) {
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
          F_b <- mean(sqb[good_b]/denom_b[good_b], na.rm=TRUE)
          F_b_per_gen <- F_b / (g1 - gen0)
          boot_ne[b] <- ifelse(F_b_per_gen > 0, 1/(2*F_b_per_gen), NA)
        }
      }
      boot_ne <- boot_ne[!is.na(boot_ne)]
      if (length(boot_ne) > 0) {
        out$Ne_lower[i] <- quantile(boot_ne, 0.025)
        out$Ne_median[i] <- median(boot_ne)
        out$Ne_upper[i] <- quantile(boot_ne, 0.975)
      }
    }
  }
  
  return(out)
}

# Ne relative to generation 2, with bootstrap
res_baseline_F <- estimateNe_from_baseline_F(Sync, nboot = 200, bootstrap = TRUE)

# View results
print(res_baseline_F)

