# summarize_ukb_analysis.R 
# Process raw analysis of UKB data
.libPaths("/home/swang25/R/ubuntu/4.4.1")  # HPC R library path

# setwd("/path/to/csmGmm_reproduce/Fig4/summarize_ukb_analysis.R") or set the path after the -cwd flag if error to locate file
here::i_am("CRC_IBD/summarize_ukb_analysis.R")
setwd("/home/swang25/other_summary_statistics/CRC_IBD/")

library(data.table)
library(dplyr)
library(devtools)
library(ks)
library(csmGmm)
library(here)

# record
args <- commandArgs(trailingOnly = TRUE)
aID  <- as.numeric(args[1])
Snum <- as.numeric(args[2])  # scenario index controlling FDR thresholds

# source the .R scripts from the SupportingCode
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)


outputDir <- here::here("CRC_IBD", "output")
dataDir   <- here::here("CRC_IBD", "Data")
fnameOut        <- paste0(outputDir, "/processed_ukb_data_S", Snum, ".txt")
rejectFnameRoot <- paste0(outputDir, "/reject_bmi_with_overall_neg5_reject_S_", Snum)

# dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)

# nominal FDR thresholds
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


fnameRoot   <- paste0(outputDir, "/Fig4_data_aID", 1:2)
fnameDACT   <- paste0(fnameRoot, "_DACTp.txt")
fnameHDMT   <- paste0(fnameRoot, "_hdmt.txt")
fnameKernel <- paste0(fnameRoot, "_kernel.txt")
fname7      <- paste0(fnameRoot, "_df7.txt")
fname50     <- paste0(fnameRoot, "_df50.txt")
fnameNew    <- paste0(fnameRoot, "_newlfdr.txt")

selections <- list()
selections[[1]] <- c("z_crc", "z_ibd")  # aID 1: harmonized independent
selections[[2]] <- c("z_crc", "z_ibd")  # aID 2: correlated (renamed below)


cummean_reject <- function(x, q) {
  o  <- order(x)
  cm <- cumsum(x[o]) / seq_along(o)
  r  <- max(which(cm <= q), 0)
  rej <- integer(length(x)); if (r > 0) rej[o[seq_len(r)]] <- 1L
  rej
}

allResults <- c()  # SUGGESTION: use list() then bind_rows for speed, but ok as-is

for (file_it in 1:2) {
  # file_it==1: harmonized independent (z_crc,z_ibd,p_crc,p_ibd)
  # file_it==2: correlated same-cohort (Zcrc/Zibd/pCRC/pIBD) -> renamed
  if (file_it == 1) {
    # aID1: harmonized, independent
    cleanZ <- fread(
      here::here(dataDir, "crc_ibd_harmonized.cleaned.tsv.gz"),
      select = c("chromosome", "base_pair_location",
                 "z_crc", "z_ibd", "p_crc", "p_ibd")
    )
    cleanZ <- cleanZ %>%
      mutate(chrpos = paste0(chromosome, ":", base_pair_location)) %>%
      select(z_crc, z_ibd, p_crc, p_ibd, chrpos)

    fdrLimitHDMT   <- fdrLimitHDMTi
    fdrLimitDACT   <- fdrLimitDACTi
    fdrLimitKernel <- fdrLimitKerneli
    fdrLimit50     <- fdrLimit50i
    fdrLimit7      <- fdrLimit7i

  } else if (file_it == 2) {
    # aID2: correlated same-cohort
    cleanZ <- fread(
      here::here(dataDir, "cor_crc_ibd.cleaned.tsv.gz"),
      select = c("chromosome", "base_pair_location",
                 "Zcrc", "Zibd", "pCRC", "pIBD")
    )
    cleanZ <- cleanZ %>%
      rename(z_crc = Zcrc, z_ibd = Zibd,
             p_crc = pCRC, p_ibd = pIBD) %>%
      mutate(chrpos = paste0(chromosome, ":", base_pair_location)) %>%
      select(z_crc, z_ibd, p_crc, p_ibd, chrpos)

    fdrLimitHDMT   <- fdrLimitHDMTr
    fdrLimitDACT   <- fdrLimitDACTr
    fdrLimitKernel <- fdrLimitKernelr
    fdrLimit50     <- fdrLimit50r
    fdrLimit7      <- fdrLimit7r
  }

  # Methods: New (csmGmm), Kernel, df50, df7, HDMT, DACT
  tempRes <- data.frame(
    Method    = c("New", "Kernel", "df50", "df7", "HDMT", "DACT"),
    numReject = NA_integer_
  )

  # --- New (csmGmm lfdr) ---
 
  tempNew <- fread(fnameNew[file_it], header = TRUE, data.table = FALSE)
  x_new   <- tempNew[[1]]  # first column = lfdr vector

  tempDat <- cleanZ %>%
    select(all_of(selections[[file_it]]), chrpos) %>%
    mutate(
      origIdx = row_number(),
      newLfdr = x_new
    ) %>%
    arrange(newLfdr) %>%
    mutate(
      cumNew = cummean(newLfdr),
      rejNew = ifelse(cumNew < fdrLimitNew, 1L, 0L)
    ) %>%
    arrange(origIdx)

  tempRes$numReject[1] <- sum(tempDat$rejNew, na.rm = TRUE)

  # --- Kernel lfdr ---
  if (file.exists(fnameKernel[file_it])) {
    tempKernel <- fread(fnameKernel[file_it], header = TRUE, data.table = FALSE)
    tempDat <- tempDat %>%
      mutate(kernelLfdr = tempKernel[[1]]) %>%
      arrange(kernelLfdr) %>%
      mutate(
        cumKernel = cummean(kernelLfdr),
        rejKernel = ifelse(cumKernel < fdrLimitKernel, 1L, 0L)
      ) %>%
      arrange(origIdx)
    tempRes$numReject[2] <- sum(tempDat$rejKernel, na.rm = TRUE)
  } else {
    # SUGGESTION: optionally record 0 instead of NA when file missing
    # tempRes$numReject[2] <- 0L
  }

  # --- df50 lfdr ---
  if (file.exists(fname50[file_it])) {
    tempdf50 <- fread(fname50[file_it], header = TRUE, data.table = FALSE)
    tempDat <- tempDat %>%
      mutate(df50Lfdr = tempdf50[[1]]) %>%
      arrange(df50Lfdr) %>%
      mutate(
        cumdf50 = cummean(df50Lfdr),
        rejdf50 = ifelse(cumdf50 < fdrLimit50, 1L, 0L)
      ) %>%
      arrange(origIdx)
    tempRes$numReject[3] <- sum(tempDat$rejdf50, na.rm = TRUE)
  }

  # --- df7 lfdr ---
  if (file.exists(fname7[file_it])) {
    tempdf7 <- fread(fname7[file_it], header = TRUE, data.table = FALSE)
    tempDat <- tempDat %>%
      mutate(df7Lfdr = tempdf7[[1]]) %>%
      arrange(df7Lfdr) %>%
      mutate(
        cumdf7 = cummean(df7Lfdr),
        rejdf7 = ifelse(cumdf7 < fdrLimit7, 1L, 0L)
      ) %>%
      arrange(origIdx)
    tempRes$numReject[4] <- sum(tempDat$rejdf7, na.rm = TRUE)
  }

  # --- HDMT ---

  if (file.exists(fnameHDMT[file_it])) {
    tempHDMT <- fread(fnameHDMT[file_it], header = TRUE, data.table = FALSE)
    if ("fixedFdr" %in% names(tempHDMT)) {
      tempDat <- tempDat %>%
        mutate(
          hdmtFdr = tempHDMT$fixedFdr,
          rejHDMT = ifelse(hdmtFdr < fdrLimitHDMT, 1L, 0L)
        ) %>%
        arrange(origIdx)
      tempRes$numReject[5] <- sum(tempDat$rejHDMT, na.rm = TRUE)
    }
  }

  # --- DACT ---
  if (file.exists(fnameDACT[file_it])) {
    tempDACT <- fread(fnameDACT[file_it], header = TRUE, data.table = FALSE)
    tempDat <- tempDat %>%
      mutate(DACTp = tempDACT[[1]]) %>%
      arrange(DACTp) %>%
      mutate(
        rankedIdxP = row_number(),
        km = rankedIdxP / n(),
        RHS = km * fdrLimitDACT
      )

    rejected <- which(tempDat$DACTp <= tempDat$RHS)
    maxIdx <- if (length(rejected)) max(rejected) else 0

    tempDat <- tempDat %>%
      mutate(rejDACT = ifelse(rankedIdxP <= maxIdx, 1L, 0L)) %>%
      arrange(origIdx)

    tempRes$numReject[6] <- sum(tempDat$rejDACT, na.rm = TRUE)
  }

  
  tempDat <- tempDat %>%
    mutate(rejAny = ifelse(
      coalesce(rejDACT,   0L) == 1L |
        coalesce(rejHDMT, 0L) == 1L |
        coalesce(rejdf7,  0L) == 1L |
        coalesce(rejdf50, 0L) == 1L |
        coalesce(rejKernel, 0L) == 1L |
        coalesce(rejNew,  0L) == 1L,
      1L, 0L
    ))

  allResults <- rbind(allResults, tempRes %>% mutate(aID = file_it))

  
  write.table(
    tempDat %>% filter(rejAny == 1),
    paste0(rejectFnameRoot, "_aID", file_it, ".txt"),
    append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
  )

  cat(file_it)  # prints progress (1 then 2)
}

# save summary table (method counts by aID)
write.table(
  allResults,
  fnameOut,
  append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"
)
