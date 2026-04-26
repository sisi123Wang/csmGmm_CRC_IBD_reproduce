#!/usr/bin/env Rscript
# Process two-way pleiotropy / correlation analysis output
#
# This script:
# 1. reads method-specific output files
# 2. applies rejection rules under preset FDR thresholds
# 3. saves rejected SNPs for each analysis
# 4. writes a summary table of rejection counts

.libPaths("/home/swang25/R/ubuntu/4.4.1")

# If here() fails, set the working directory manually to the project folder
# or set -cwd correctly in the .lsf file.
here::i_am("Fig1/two-way pleiotropy analysis.R")

suppressPackageStartupMessages({
  library(purrr)
  library(readr)
  library(data.table)
  library(dplyr)
  library(devtools)
  library(ks)
  library(csmGmm)
  library(here)
})

# ---------------------------
# Command-line arguments
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)
aID  <- as.numeric(args[1])   # currently not used below, but kept for compatibility
Snum <- as.numeric(args[2])

# ---------------------------
# Source helper scripts
# ---------------------------
codePath <- here::here("SupportingCode")
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(file.path(codePath, toBeSourced), source)

# ---------------------------
# Directories and output names
# ---------------------------
outputDir <- here::here("Fig1", "output")
dataDir   <- here::here("Data")

fnameOut <- file.path(outputDir, paste0("processed_ukb_data_S", Snum, ".txt"))
rejectFnameRoot <- file.path(outputDir, paste0("reject_bmi_with_overall_neg5_reject_S_", Snum))

# ---------------------------
# Nominal FDR thresholds
# Snum = 1 uses a stricter setup for several methods
# ---------------------------
if (Snum == 1) {
  fdrLimitHDMTi   <- 0.1
  fdrLimitDACTi   <- 0.01
  fdrLimitKerneli <- 0.01
  fdrLimit7i      <- 0.01
  fdrLimit50i     <- 0.01
  fdrLimitNew     <- 0.1

  fdrLimitHDMTr   <- 0.01
  fdrLimitDACTr   <- 0.01
  fdrLimitKernelr <- 0.01
  fdrLimit7r      <- 0.01
  fdrLimit50r     <- 0.01
} else {
  fdrLimitHDMTi   <- 0.1
  fdrLimitDACTi   <- 0.1
  fdrLimitKerneli <- 0.1
  fdrLimit7i      <- 0.1
  fdrLimit50i     <- 0.1
  fdrLimitNew     <- 0.1

  fdrLimitHDMTr   <- 0.1
  fdrLimitDACTr   <- 0.1
  fdrLimitKernelr <- 0.1
  fdrLimit7r      <- 0.1
  fdrLimit50r     <- 0.1
}

# ---------------------------
# Output file roots
# ---------------------------
fnameRoot   <- paste0(outputDir, "/Fig4_data_aID", 1:9)
fnameDACT   <- paste0(fnameRoot, "_DACTp.txt")
fnameHDMT   <- paste0(fnameRoot, "_hdmt.txt")
fnameKernel <- paste0(fnameRoot, "_kernel.txt")
fname7      <- paste0(fnameRoot, "_df7.txt")
fname50     <- paste0(fnameRoot, "_df50.txt")   # this must exist if used below
fnameNew    <- paste0(fnameRoot, "_newlfdr.txt")

# ---------------------------
# Which Z columns correspond to each aID/file
# ---------------------------
selections <- list(
  c("Zcad", "Zbmi"),
  c("Zoverall", "Zcad"),
  c("Zoverall", "Zbmi"),
  c("Zoverall", "Zlcukb"),
  c("Zcad_cardio", "Zcadukb"),
  c("Zoverall", "Zcad", "Zbmi"),
  c("Zlcukb", "Zcadukb"),
  c("Zlcukb", "Zbmi"),
  c("Zcadukb", "Zbmi")
)

# ---------------------------
# Main loop
# ---------------------------
allResults <- data.frame()

for (file_it in 1:9) {

  # -------------------------
  # Load the matching cleaned Z data
  # -------------------------
  if (file_it == 4) {
    cleanZ <- fread(here::here(dataDir, "replication_with_lcoverall.txt"))
    fdrLimitHDMT   <- fdrLimitHDMTr
    fdrLimitDACT   <- fdrLimitDACTr
    fdrLimitKernel <- fdrLimitKernelr
    fdrLimit50     <- fdrLimit50r
    fdrLimit7      <- fdrLimit7r
  } else if (file_it == 5) {
    cleanZ <- fread(here::here(dataDir, "cad_for_replication.txt"))
    fdrLimitHDMT   <- fdrLimitHDMTr
    fdrLimitDACT   <- fdrLimitDACTr
    fdrLimitKernel <- fdrLimitKernelr
    fdrLimit50     <- fdrLimit50r
    fdrLimit7      <- fdrLimit7r
  } else {
    cleanZ <- fread(here::here(dataDir, "bmi_with_overall.txt"))
    fdrLimitHDMT   <- fdrLimitHDMTi
    fdrLimitDACT   <- fdrLimitDACTi
    fdrLimitKernel <- fdrLimitKerneli
    fdrLimit50     <- fdrLimit50i
    fdrLimit7      <- fdrLimit7i
  }

  # Store rejection counts for this analysis
  tempRes <- data.frame(
    Method = c("New", "Kernel", "df50", "df7", "HDMT", "DACT"),
    numReject = NA
  )

  # -------------------------
  # New method
  # Rejection rule: rank lfdr, take cumulative mean, reject while below threshold
  # -------------------------
  tempNew <- fread(fnameNew[file_it], header = TRUE, data.table = FALSE)

  tempDat <- cleanZ %>%
    select(all_of(selections[[file_it]]), chrpos) %>%
    mutate(origIdx = 1:nrow(.)) %>%
    mutate(newLfdr = tempNew$x) %>%
    arrange(newLfdr) %>%
    mutate(cumNew = cummean(newLfdr),
           rejNew = ifelse(cumNew < fdrLimitNew, 1, 0)) %>%
    arrange(origIdx)

  tempRes$numReject[1] <- sum(tempDat$rejNew)

  # For file_it 7:9, only the new method is processed
  if (file_it >= 7) {
    rejectDat <- tempDat %>% filter(rejNew == 1)

    allResults <- rbind(allResults, tempRes %>% mutate(aID = file_it))

    write.table(
      rejectDat,
      paste0(rejectFnameRoot, "_aID", file_it, ".txt"),
      append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
    )
    next
  }

  # -------------------------
  # Kernel method
  # -------------------------
  tempKernel <- fread(fnameKernel[file_it], header = TRUE, data.table = FALSE)

  tempDat <- tempDat %>%
    mutate(kernelLfdr = tempKernel$x) %>%
    arrange(kernelLfdr) %>%
    mutate(cumKernel = cummean(kernelLfdr),
           rejKernel = ifelse(cumKernel < fdrLimitKernel, 1, 0)) %>%
    arrange(origIdx)

  tempRes$numReject[2] <- sum(tempDat$rejKernel)

  # -------------------------
  # df50 method
  # -------------------------
  tempdf50 <- fread(fname50[file_it], header = TRUE, data.table = FALSE)

  tempDat <- tempDat %>%
    mutate(df50Lfdr = tempdf50$x) %>%
    arrange(df50Lfdr) %>%
    mutate(cumdf50 = cummean(df50Lfdr),
           rejdf50 = ifelse(cumdf50 < fdrLimit50, 1, 0)) %>%
    arrange(origIdx)

  tempRes$numReject[3] <- sum(tempDat$rejdf50)

  # -------------------------
  # df7 method
  # -------------------------
  tempdf7 <- fread(fname7[file_it], header = TRUE, data.table = FALSE)

  tempDat <- tempDat %>%
    mutate(df7Lfdr = tempdf7$x) %>%
    arrange(df7Lfdr) %>%
    mutate(cumdf7 = cummean(df7Lfdr),
           rejdf7 = ifelse(cumdf7 < fdrLimit7, 1, 0)) %>%
    arrange(origIdx)

  tempRes$numReject[4] <- sum(tempDat$rejdf7)

  # For file_it 4:6, only New + Kernel + df50 + df7 are processed
  if (file_it >= 4 & file_it <= 6) {
    tempDat <- tempDat %>%
      mutate(rejAny = ifelse(
        rejdf7 == 1 | rejdf50 == 1 | rejKernel == 1 | rejNew == 1, 1, 0
      ))

    rejectDat <- tempDat %>% filter(rejAny == 1)

    allResults <- rbind(allResults, tempRes %>% mutate(aID = file_it))

    write.table(
      rejectDat,
      paste0(rejectFnameRoot, "_aID", file_it, ".txt"),
      append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
    )
    next
  }

  # -------------------------
  # HDMT method
  # -------------------------
  tempHDMT <- fread(fnameHDMT[file_it], header = TRUE, data.table = FALSE)

  tempDat <- tempDat %>%
    mutate(hdmtFdr = tempHDMT$fixedFdr,
           rejHDMT = ifelse(hdmtFdr < fdrLimitHDMT, 1, 0)) %>%
    arrange(origIdx)

  tempRes$numReject[5] <- sum(tempDat$rejHDMT)

  # -------------------------
  # DACT method
  # BH-style thresholding on ordered p-values
  # -------------------------
  tempDACT <- fread(fnameDACT[file_it], header = TRUE, data.table = FALSE)

  tempDat <- tempDat %>%
    mutate(DACTp = tempDACT$x) %>%
    arrange(DACTp) %>%
    mutate(
      rankedIdxP = 1:nrow(.),
      km = (1:nrow(.)) / nrow(.),
      RHS = km * fdrLimitDACT
    )

  rejected <- which(tempDat$DACTp <= tempDat$RHS)
  maxIdx <- if (length(rejected) == 0) 0 else max(rejected)

  tempDat <- tempDat %>%
    mutate(rejDACT = ifelse(rankedIdxP <= maxIdx, 1, 0)) %>%
    arrange(origIdx)

  tempRes$numReject[6] <- sum(tempDat$rejDACT)

  # -------------------------
  # Save SNPs rejected by any method
  # -------------------------
  tempDat <- tempDat %>%
    mutate(rejAny = ifelse(
      rejDACT == 1 | rejHDMT == 1 | rejdf7 == 1 |
      rejdf50 == 1 | rejKernel == 1 | rejNew == 1, 1, 0
    ))

  rejectDat <- tempDat %>% filter(rejAny == 1)

  allResults <- rbind(allResults, tempRes %>% mutate(aID = file_it))

  write.table(
    rejectDat,
    paste0(rejectFnameRoot, "_aID", file_it, ".txt"),
    append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
  )

  cat(file_it)
}

# ---------------------------
# Save summary table
# ---------------------------
write.table(
  allResults,
  fnameOut,
  append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)
