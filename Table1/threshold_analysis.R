#!/usr/bin/env Rscript
# ------------------------------------------------------------
#        CRC–IBD  Simulation  (sample-template, cleaned)
# ------------------------------------------------------------

.libPaths("/home/swang25/R/ubuntu/4.4.1")
suppressPackageStartupMessages({
  library(MASS)       # mvrnorm
  library(data.table)
  library(dplyr)      # cummean()
  library(csmGmm)
  library(here)
})


## ---------- user switches -------------------------------------
replication_mode <- FALSE   # FALSE = pleiotropy

## ---------- command-line ID -----------------------------------
args   <- commandArgs(trailingOnly = TRUE)
aID    <- if (length(args)) as.numeric(args[1]) else 1
set.seed(aID)

## ---------- constants -----------------------------------------
nSNPs <- 9234567
#nSNPs <- 1e6

Sigma <- diag(2)            # independent Z noise

# three scenario π-vectors (from sample code)
sProp <- list(
  s1 = c("0,0"=0.98,  "1,0"=0.005, "-1,0"=0.005, "0,1"=0.005, "0,-1"=0.005,
         "1,1"=0.0001,"-1,-1"=0.0001,"1,-1"=0.0001,"-1,1"=0.0001),
  s2 = c("0,0"=0.95,  "1,0"=0.01,  "-1,0"=0.01,  "0,1"=0.015, "0,-1"=0.015,
         "1,1"=0.0005,"-1,-1"=0.0005,"1,-1"=0.0005,"-1,1"=0.0005),
  s3 = c("0,0"=0.90,  "1,0"=0.02,  "-1,0"=0.02,  "0,1"=0.025, "0,-1"=0.025,
         "1,1"=0.001, "-1,-1"=0.001, "1,-1"=0.001, "-1,1"=0.001)
)


## ---------- empirical μ ---------------------------------------
library(data.table)

mu_file <- here::here("CRC_IBD","output","Fig4_data_aID1_muInfo.txt")
mu_emp  <- as.matrix(fread(mu_file, header = TRUE, data.table = FALSE)) #changed header to TRUE
storage.mode(mu_emp) <- "numeric"

first_nz <- function(x) { i <- which(is.finite(x) & x != 0); if (length(i)) i[1] else NA_integer_ }
last_nz  <- function(x) { i <- which(is.finite(x) & x != 0); if (length(i)) i[length(i)] else NA_integer_ }

c1_single <- first_nz(mu_emp[1, ])
c1_shared <- last_nz (mu_emp[1, ])
c2_single <- first_nz(mu_emp[2, ])
c2_shared <- last_nz (mu_emp[2, ])

stopifnot(all(is.finite(c(c1_single, c1_shared, c2_single, c2_shared))))

single_crc <- mu_emp[1, c1_single]  # from row1 first non-zero 
shared_crc <- mu_emp[1, c1_shared]  # from row1 last  non-zero
single_ibd <- mu_emp[2, c2_single]  # from row2 first non-zero
shared_ibd <- mu_emp[2, c2_shared]  # from row2 last  non-zero

print(list(single_crc, shared_crc, single_ibd, shared_ibd))

## ---------- results frame -------------------------------------
powerRes <- data.frame(
  scenario      = names(sProp),
  nCausal       = NA,  propCausal = NA,
  nRejNew       = NA,  nRejMeta   = NA,
  nRejThreshold1= NA,  nRejThreshold2 = NA, nRejThreshold3 = NA,
  powerNew      = NA,  powerMeta  = NA,
  powerThreshold1 = NA, powerThreshold2 = NA, powerThreshold3 = NA,
  fdpNew        = NA,  fdpMeta    = NA,
  fdpThreshold1 = NA,  fdpThreshold2 = NA,  fdpThreshold3 = NA,
  seed          = aID
)

## ---------- loop over the three π-scenarios -------------------
for (i in seq_along(sProp)) {
  
  ## ----- pattern catalogue {-1,0,+1}^2 ------------------------
  H <- CJ(h1 = -1:1, h2 = -1:1)
  H[, mu1 := fifelse(h1 == 0, 0,
                     sign(h1) * ifelse(h2 == 0, single_crc, shared_crc))]
  H[, mu2 := fifelse(h2 == 0, 0,
                     sign(h2) * ifelse(h1 == 0, single_ibd, shared_ibd))]
  
  ## map named sProp to this CJ order
  lab    <- paste(H$h1, H$h2, sep = ",")     # "-1,-1", "-1,0", ...
  pi_vec <- as.numeric(sProp[[i]][lab])
  pi_vec[is.na(pi_vec)] <- 0
  pi_vec <- pi_vec / sum(pi_vec)
  sc_name <- names(sProp)[i]

  ## ----- simulate Z ------------------------------------------
  pat <- sample.int(9, nSNPs, TRUE, pi_vec)        # *only this pi_vec*
  mu  <- cbind(H$mu1[pat], H$mu2[pat])
  
  Z   <- mu + mvrnorm(nSNPs, c(0,0), Sigma)
  Z[Z >  8.1] <-  8.1;  Z[Z < -8.1] <- -8.1
  P   <- 1 - pchisq(Z^2, 1)
  is_both <- (mu[,1] != 0 & mu[,2] != 0)

  
  ## ---------- quick sanity print (comment when happy) ---------
  # cat(sc_name, ":  π₀ =", round(mean(!causal), 3),
  #     "| shared =", round(mean(rowSums(abs(mu))==shared_crc+shared_ibd),4), "\n")
  
  ## ----- thresholds --------------------------------
  sameDir <- (Z[,1] > 0 & Z[,2] > 0) | (Z[,1] < 0 & Z[,2] < 0)
  pass <- function(cut) if (replication_mode)
    which(sameDir & P[,1] < cut & P[,2] < cut)
  else
    which(          P[,1] < cut & P[,2] < cut)
  
  idx <- list(t1 = pass(1e-8),
              t2 = pass(1e-6),
              t3 = pass(1e-5))

  ## ----- csmGmm -----------------------------------------------
  initPi <- list(c(0.82), c(0.02,0.02), c(0.02,0.02), c(0.10))
  initMu <- list(matrix(0,2,1),
                 matrix(c(0,3, 0,6),2),
                 matrix(c(3,0, 6,0),2),
                 matrix(c(8,8),2))
  
  lfdr <- tryCatch(
    symm_fit_ind_EM(Z, initMu, initPi,
                    sameDirAlt = replication_mode,
                    eps = 1e-5)$lfdrResults,
    error = function(e) rep(1, nSNPs))
  ord        <- order(lfdr)
  idx$EM     <- ord[cummean(lfdr[ord]) <= 0.10]
  
  ## ----- meta -------------------------------------------------
  w1 <- sqrt(446766); w2 <- sqrt(446865)
  Pmeta  <- 2*pnorm(-abs((Z[,1]*w1 + Z[,2]*w2)/sqrt(w1^2+w2^2)))
  idx$meta <- which(Pmeta < 1e-8)
  
  ## ----- summarise -------------------------------------------
  summarise <- function(v, truth) {
    n_rej <- length(v)
    n_alt <- sum(truth)
    TP <- sum(truth[v])
    FP <- n_rej - TP
    c(
      n_rej,
      if (n_alt == 0) NA_real_ else TP / n_alt,   # power wrt BOTH-traits only
      if (n_rej == 0) NA_real_ else FP / n_rej    # FDP
    )
  }
  stats <- list(
    t1   = summarise(idx$t1,   is_both),
    t2   = summarise(idx$t2,   is_both),
    t3   = summarise(idx$t3,   is_both),
    EM   = summarise(idx$EM,   is_both),
    meta = summarise(idx$meta, is_both)
  )

  powerRes[i, c("nRejThreshold1","powerThreshold1","fdpThreshold1")] <- stats$t1
  powerRes[i, c("nRejThreshold2","powerThreshold2","fdpThreshold2")] <- stats$t2
  powerRes[i, c("nRejThreshold3","powerThreshold3","fdpThreshold3")] <- stats$t3
  powerRes[i, c("nRejNew","powerNew","fdpNew")]                     <- stats$EM
  powerRes[i, c("nRejMeta","powerMeta","fdpMeta")]                  <- stats$meta
  
  powerRes$nCausal[i]    <- sum(is_both)     # count BOTH-traits
  powerRes$propCausal[i] <- mean(is_both)
  
  cat("scenario", sc_name, "done | causal prop:",
      round(powerRes$propCausal[i],4), "\n")
}

print(powerRes)

## ---------- write ---------------------------------------------
out_dir <- here::here("CRC_IBD", "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(out_dir, sprintf("pleio_approach_aID%d.txt", aID))
fwrite(powerRes, outfile, sep = '\t')

print(powerRes)
cat("\nSaved to", outfile, "\n")
