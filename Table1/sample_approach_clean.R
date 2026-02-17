#!/usr/bin/env Rscript
# sample_approach_clean.R

.libPaths("/home/swang25/R/ubuntu/4.4.1")  # HPC-specific library path

setwd("/home/swang25/other_summary_statistics/CRC_IBD/")
here::i_am("sample_approach_clean.R") 


library(data.table)
library(dplyr)
library(devtools)
library(ks)
library(csmGmm)
library(here)
library(locfdr)
library(purrr)

# args (HPC job array / command line)
# aID: analysis ID (1 = pleiotropy independent, 2 = correlation)
# Snum: unused below; keep if other scripts use it

args  <- commandArgs(trailingOnly = TRUE)
aID   <- as.integer(args[1])
Snum  <- as.integer(args[2]); if (is.na(Snum)) Snum <- 1  # SUGGESTION: remove if truly unused

# source helpers from SupportingCode/
codePath <- here::here("SupportingCode")
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

outputDir <- here::here("output")
fnameRoot <- file.path(outputDir, sprintf("Fig4_data_aID%d", aID))
summaryStatDir <- here::here("Data")  # SUGGESTION: name implies directory, good

# EM tolerance
oldEps <- 0.01
newEps <- 1e-5   # tighter EM stopping criterion

# flags
# cor = FALSE -> pleiotropy independent branch (needs P-values matrix P)
# cor = TRUE  -> correlated branch (X only)
cor <- FALSE

# aID == 1: pleiotropy (independent); uses both Z and P for HDMT/DACT/EB
# aID == 2: correlation; uses Z only + estimated rho
if (aID == 1) {
  # Pleiotropy (independent)
  cleanUKB <- fread(
    here::here(summaryStatDir, "crc_ibd_harmonized.cleaned.tsv.gz"),
    select = c("z_crc", "z_ibd", "p_crc", "p_ibd")
  )
  X <- as.matrix(cleanUKB[, .(z_crc, z_ibd)])  # test statistics (Z)
  P <- as.matrix(cleanUKB[, .(p_crc, p_ibd)])  # p-values (needed for null_estimation/DACT)
} else if (aID == 2) {
  # Correlation
  cor <- TRUE
  cleanUKB <- fread(
    here::here(summaryStatDir, "cor_crc_ibd.cleaned.tsv.gz"),
    select = c("Zcrc", "Zibd")
  )
  setnames(cleanUKB, c("Zcrc", "Zibd"), c("z_crc", "z_ibd"))
  X <- as.matrix(cleanUKB[, .(z_crc, z_ibd)])
} else {
  stop("aID only 1 and 2")
}

# drop non-finite values
keep <- is.finite(X[, 1]) & is.finite(X[, 2])
X <- X[keep, , drop = FALSE]
if (!cor) P <- P[keep, , drop = FALSE]
stopifnot(nrow(X) > 0L)

# ---------------- HDMT & DACT (pleiotropy only) ----------------
if (!cor) {
  # null proportion estimation from p-values
  nullprop <- tryCatch(null_estimation(P), error = function(e) e, warning = function(w) w)

  # NOTE: your logic treats warnings like errors because warning=function(w) w returns a "warning" object.
  # SUGGESTION (later): in production, don't capture warnings as return values; just let warnings print.

  if (is.list(nullprop)) {
    hdmtOut <- tryCatch(
      fdr_est(
        nullprop$alpha00, nullprop$alpha01, nullprop$alpha10,
        nullprop$alpha1, nullprop$alpha2, P
      ),
      error = function(e) e, warning = function(w) w
    )
  } else {
    hdmtOut <- rep(NA_real_, nrow(P))
  }

  # Ensure a dataframe is written even on failure
  if (!is.data.frame(hdmtOut)) hdmtOut <- data.frame(fixedFdr = rep(NA_real_, nrow(P)))

  write.table(
    hdmtOut,
    paste0(fnameRoot, "_hdmt.txt"),
    append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE
  )

  # DACT: handle errors only; don't return warnings as values
  DACTout <- tryCatch(
    DACT(
      p_a = P[, 1],
      p_b = P[, 2],
      nullEst = if (is.list(nullprop)) nullprop else NULL,
      correction = "JC"
    ),
    error = function(e) NULL
  )

  # Only accept real p_dact; otherwise write NAs so file isn't header-only
  DACTfreqp <- if (!is.null(DACTout) && !is.null(DACTout$p_dact)) {
    DACTout$p_dact
  } else {
    rep(NA_real_, nrow(P))
  }

  write.table(
    DACTfreqp,
    paste0(fnameRoot, "_DACTp.txt"),
    append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE
  )
}

# ---------------- EB kernel / df (pleiotropy only) ----------------
if (!cor) {
  # Empirical Bayes: kernel density version
  kfit <- emp_bayes_framework(
    summary_tab = X,
    sameDirAlt = FALSE,
    kernel = TRUE, joint = FALSE, ind = TRUE,
    dfFit = 7, Hdist_epsilon = 1e-2, checkpoint = TRUE
  )
  write.table(
    if (is.list(kfit)) kfit$lfdrVec else rep(NA_real_, nrow(X)),
    paste0(fnameRoot, "_kernel.txt"),
    append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE
  )

  # Empirical Bayes: t mixture df=7
  fit7 <- emp_bayes_framework(
    summary_tab = X,
    sameDirAlt = FALSE,
    kernel = FALSE, joint = FALSE, ind = TRUE,
    dfFit = 7, Hdist_epsilon = 1e-2, checkpoint = TRUE
  )
  write.table(
    if (is.list(fit7)) fit7$lfdrVec else rep(NA_real_, nrow(X)),
    paste0(fnameRoot, "_df7.txt"),
    append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE
  )

  # Empirical Bayes: t mixture df=50
  fit50 <- emp_bayes_framework(
    summary_tab = X,
    sameDirAlt = FALSE,
    kernel = FALSE, joint = FALSE, ind = TRUE,
    dfFit = 50, Hdist_epsilon = 1e-2, checkpoint = TRUE
  )
  write.table(
    if (is.list(fit50)) fit50$lfdrVec else rep(NA_real_, nrow(X)),
    paste0(fnameRoot, "_df50.txt"),
    append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE
  )
}

# ---------------- EM (2D): csmGmm symmetric mixture ----------------
# Initialization (pi: mixture weights; mu: component means)
initPiList <- list(
  c(0.97),
  c(0.01, 0.01),
  c(0.01, 0.01),
  c(0.01)
)

initMuList <- list(
  matrix(0, 2, 1),
  matrix(c(0, 1, 0, 2), 2),
  matrix(c(1, 0, 2, 0), 2),
  matrix(c(2, 2), 2, 1)
)

if (!cor) {
  # Independent version
  newRes <- symm_fit_ind_EM(
    testStats = X,
    sameDirAlt = FALSE,
    initMuList = initMuList,
    initPiList = initPiList,
    eps = newEps
  )
} else {
  # Correlated version: estimate rho from Zs
  rho <- stats::cor(X[, 1], X[, 2])
  rho <- max(min(rho, 0.99), -0.99)  # guard against singular covariance
  C <- matrix(c(1, rho, rho, 1), 2)

  newRes <- symm_fit_cor_EM(
    testStats = X,
    corMat = C,
    initMuList = initMuList,
    initPiList = initPiList,
    eps = newEps
  )
}

# save outputs
# make sure it exist: dir.create(outputDir, recursive=TRUE)
write.table(
  newRes$lfdrResults,
  paste0(fnameRoot, "_newlfdr.txt"),
  append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE
)
write.table(
  do.call(cbind, newRes$muInfo),
  paste0(fnameRoot, "_muInfo.txt"),
  append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE
)
write.table(
  do.call(cbind, newRes$piInfo),
  paste0(fnameRoot, "_piInfo.txt"),
  append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE
)


# You may also want to save rho (aID==2) for reproducibility:
# writeLines(as.character(rho), paste0(fnameRoot, "_rho.txt"))
