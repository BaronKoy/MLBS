#!/usr/bin/env Rscript

suppressMessages({
  library(data.table)
})

# --- User parameters ---
sync_file   <- "/home/baron/Documents/PhD/Data/pop_size_analysis/PoolSeq/cages/cage_2.sync"
output_file <- "Ne_estimates.csv"
allele_target <- "A"          # base to track (A,C,G,T)
min_cov <- 10
max_cov <- 200
min_maf <- 0.01
chunk_size <- 1e5             # process 100k loci per chunk

time_points <- c(2,4,8,12,20,28,36,44,56)
intervals <- data.frame(
  start_gen = time_points[-length(time_points)],
  end_gen   = time_points[-1],
  delta_t   = diff(time_points),
  sum_F     = 0,
  count     = 0
)

# --- helper functions ---
get_freq_cov <- function(count_str, allele="A"){
  counts <- as.numeric(unlist(strsplit(count_str, ":")))
  if (length(counts) < 4) return(c(NA,0))
  total <- sum(counts[1:4])
  if (total == 0) return(c(NA,0))
  allele_index <- switch(allele, "A"=1, "C"=2, "G"=3, "T"=4)
  freq <- counts[allele_index] / total
  return(c(freq, total))
}

estimate_chunk <- function(chunk, allele){
  pop_cols <- 4:ncol(chunk)
  res_list <- vector("list", length(intervals$start_gen))
  for (i in seq_along(res_list)){
    res_list[[i]] <- list(sumF=0, n=0)
  }
  
  for (r in 1:nrow(chunk)){
    vals <- sapply(chunk[r, pop_cols, with=FALSE], get_freq_cov, allele=allele)
    freqs <- vals[1,]
    covs  <- vals[2,]
    
    # coverage + maf filtering
    if (any(is.na(freqs)) || any(covs < min_cov) || any(covs > max_cov)) next
    mean_maf <- mean(freqs)
    if (mean_maf < min_maf || mean_maf > (1 - min_maf)) next
    
    for (i in seq_along(res_list)){
      p0 <- freqs[i]
      p1 <- freqs[i+1]
      if (is.na(p0) || is.na(p1) || p0==0 || p0==1 || p1==0 || p1==1) next
      F_i <- ((p1 - p0)^2) / (p0*(1 - p0))
      if (!is.finite(F_i)) next
      res_list[[i]]$sumF <- res_list[[i]]$sumF + F_i
      res_list[[i]]$n <- res_list[[i]]$n + 1
    }
  }
  return(res_list)
}

# --- streaming loop ---
cat("Streaming sync file by chunks...\n")

rows_read <- 0
repeat {
  chunk <- tryCatch(fread(sync_file, header=FALSE, sep="\t",
                          skip=rows_read, nrows=chunk_size, colClasses="character"),
                    error=function(e) NULL)
  if (is.null(chunk) || nrow(chunk)==0) break
  if (rows_read == 0) colnames(chunk)[1:3] <- c("chr","pos","ref")
  
  partial <- estimate_chunk(chunk, allele_target)
  
  for (i in seq_along(partial)){
    intervals$sum_F[i] <- intervals$sum_F[i] + partial[[i]]$sumF
    intervals$count[i] <- intervals$count[i] + partial[[i]]$n
  }
  
  rows_read <- rows_read + nrow(chunk)
  cat("Processed", rows_read, "lines...\n")
}

# --- final Ne calculation ---
intervals$mean_F <- intervals$sum_F / pmax(intervals$count,1)
intervals$Ne <- 1 / (2 * intervals$mean_F * intervals$delta_t)

# --- save & show ---
fwrite(intervals[, c("start_gen","end_gen","delta_t","count","Ne")], output_file)
cat("\n✅ Finished. Results written to", output_file, "\n")
print(intervals[, c("start_gen","end_gen","Ne")])