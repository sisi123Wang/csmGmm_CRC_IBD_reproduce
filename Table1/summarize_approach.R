# Process raw analysis of UKB data

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig4/summarize_approach.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("Fig4/summarize_approach.R")

# load libraries
library(data.table)
library(dplyr)
library(devtools)
library(ks)
library(csmGmm)
library(here)

# record input - controls seed, parameters, etc.
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("output")
dataDir   <- here::here("Data")
rejectFnameRoot <- file.path(outputDir, paste0("reject_crc_ibd_S_", Snum))
fnameOut <- paste0(outputDir, "/processed_ukb_data_S", Snum, ".txt")

# nominal fdr
if (Snum == 1) {
  fdrLimitHDMTi <- 0.1
  fdrLimitDACTi <- 0.01
  fdrLimitKerneli <- 0.01
  fdrLimit7i <- 0.01
  fdrLimit50i <- 0.01
  fdrLimitNew <- 0.1
  fdrLimitHDMTr <- 0.01
  fdrLimitDACTr <- 0.01
  fdrLimitKernelr <- 0.01
  fdrLimit7r <- 0.01
  fdrLimit50r <- 0.01
} else {
  fdrLimitHDMTi <- 0.1
  fdrLimitDACTi <- 0.1
  fdrLimitKerneli <- 0.1
  fdrLimit7i <- 0.1
  fdrLimit50i <- 0.1
  fdrLimitNew <- 0.1
  fdrLimitHDMTr <- 0.1
  fdrLimitDACTr <- 0.1
  fdrLimitKernelr <- 0.1
  fdrLimit7r <- 0.1
  fdrLimit50r <- 0.1
}

# raw output files
ids <- c(1, 2)
fnameRoot  <- file.path(outputDir, sprintf("Fig4_data_aID%d", ids))
fnameDACT  <- paste0(fnameRoot, "_DACTp.txt")
fnameHDMT  <- paste0(fnameRoot, "_hdmt.txt")
fnameKernel<- paste0(fnameRoot, "_kernel.txt")
fname7     <- paste0(fnameRoot, "_df7.txt")
fname50    <- paste0(fnameRoot, "_df50.txt")
fnameNew   <- paste0(fnameRoot, "_newlfdr.txt")

selections <- list(
  c("z_crc","z_ibd"),  # aID=1
  c("z_crc","z_ibd")   # aID=2
)


# results
allResults <- data.frame()
for (i in seq_along(ids)) {
  file_it <- ids[i]
  
  if (file_it == 1) {
    cleanZ <- fread(here::here(dataDir, "crc_ibd_harmonized.cleaned.tsv.gz"))
    fdrLimitHDMT <- fdrLimitHDMTi
    fdrLimitDACT <- fdrLimitDACTi
    fdrLimitKernel <- fdrLimitKerneli
    fdrLimit50 <- fdrLimit50i
    fdrLimit7  <- fdrLimit7i
  } else { # file_it == 2 (correlated case)
    cleanZ <- fread(here::here(dataDir, "cor_crc_ibd.cleaned.tsv.gz"))
    fdrLimitHDMT <- fdrLimitHDMTi   # won’t be used;
    fdrLimitDACT <- fdrLimitDACTi   # won’t be used;
    fdrLimitKernel <- fdrLimitKerneli
    fdrLimit50 <- fdrLimit50i
    fdrLimit7  <- fdrLimit7i
  }
  
  # ensure there is a join key; fabricate if missing
  if (!"chrpos" %in% names(cleanZ)) {
    if (all(c("CHR","POS") %in% names(cleanZ))) {
      cleanZ[, chrpos := paste0(CHR, ":", POS)]
    } else {
      cleanZ[, chrpos := as.character(seq_len(.N))]
    }
  }
  
  # ---- NEW ----
  tempNew <- fread(fnameNew[i], header=TRUE, data.table=FALSE)
  tempRes <- data.frame(Method=c("New","Kernel","df50","df7","HDMT","DACT"), numReject=NA_real_)
  
  tempDat <- cleanZ %>%
    dplyr::select(dplyr::all_of(selections[[i]]), chrpos) %>%
    dplyr::mutate(origIdx = dplyr::row_number(),
                  newLfdr = tempNew$x) %>%
    dplyr::arrange(newLfdr) %>%
    dplyr::mutate(cumNew = cummean(newLfdr),
                  rejNew = as.integer(cumNew < fdrLimitNew)) %>%
    dplyr::arrange(origIdx)
  tempRes$numReject[1] <- sum(tempDat$rejNew, na.rm=TRUE)
  
  # For aID=2 (correlated), stop after NEW (skip HDMT/DACT)
  if (file_it == 2) {
    rejectDat <- dplyr::filter(tempDat, rejNew == 1)
    allResults <- rbind(allResults, transform(tempRes, aID=file_it))
    write.table(rejectDat, sprintf("%s_aID%d.txt", rejectFnameRoot, file_it),
                quote=FALSE, sep='\t', row.names=FALSE)
    next
  }
  
  # ---- kernel / df50 / df7 (present for aID=1) ----
  tempKernel <- fread(fnameKernel[i], header=TRUE, data.table=FALSE)
  tempDat <- tempDat %>%
    dplyr::mutate(kernelLfdr = tempKernel$x) %>%
    dplyr::arrange(kernelLfdr) %>%
    dplyr::mutate(cumKernel = cummean(kernelLfdr),
                  rejKernel = as.integer(cumKernel < fdrLimitKernel)) %>%
    dplyr::arrange(origIdx)
  tempRes$numReject[2] <- sum(tempDat$rejKernel, na.rm=TRUE)
  
  tempdf50 <- fread(fname50[i], header=TRUE, data.table=FALSE)
  tempDat <- tempDat %>%
    dplyr::mutate(df50Lfdr = tempdf50$x) %>%
    dplyr::arrange(df50Lfdr) %>%
    dplyr::mutate(cumdf50 = cummean(df50Lfdr),
                  rejdf50 = as.integer(cumdf50 < fdrLimit50)) %>%
    dplyr::arrange(origIdx)
  tempRes$numReject[3] <- sum(tempDat$rejdf50, na.rm=TRUE)
  
  tempdf7 <- fread(fname7[i], header=TRUE, data.table=FALSE)
  tempDat <- tempDat %>%
    dplyr::mutate(df7Lfdr = tempdf7$x) %>%
    dplyr::arrange(df7Lfdr) %>%
    dplyr::mutate(cumdf7 = cummean(df7Lfdr),
                  rejdf7 = as.integer(cumdf7 < fdrLimit7)) %>%
    dplyr::arrange(origIdx)
  tempRes$numReject[4] <- sum(tempDat$rejdf7, na.rm=TRUE)
  
  # ---- HDMT (aID=1 only) ----
  tempHDMT <- fread(fnameHDMT[i], header=TRUE, data.table=FALSE)
  tempDat <- tempDat %>%
    dplyr::mutate(hdmtFdr = tempHDMT$fixedFdr,
                  rejHDMT = as.integer(hdmtFdr < fdrLimitHDMT)) %>%
    dplyr::arrange(origIdx)
  tempRes$numReject[5] <- sum(tempDat$rejHDMT, na.rm=TRUE)
  
  # ---- DACT (aID=1 only) ----
  tempDACT <- fread(fnameDACT[i], header=TRUE, data.table=FALSE)
  tempDat <- tempDat %>%
    dplyr::mutate(DACTp = tempDACT$x) %>%
    dplyr::arrange(DACTp) %>%
    dplyr::mutate(rankedIdxP = dplyr::row_number(),
                  km = rankedIdxP / n(),
                  RHS = km * fdrLimitDACT)
  rejected <- which(tempDat$DACTp <= tempDat$RHS)
  maxIdx <- if (length(rejected)) max(rejected) else 0
  tempDat <- tempDat %>%
    dplyr::mutate(rejDACT = as.integer(rankedIdxP <= maxIdx)) %>%
    dplyr::arrange(origIdx)
  tempRes$numReject[6] <- sum(tempDat$rejDACT, na.rm=TRUE)
  
  # ---- any rejection & save ----
  tempDat <- tempDat %>%
    dplyr::mutate(rejAny = as.integer(rowSums(cbind(
      rejDACT, rejHDMT, rejdf7, rejdf50, rejKernel, rejNew
    )) > 0))
  rejectDat <- dplyr::filter(tempDat, rejAny == 1)
  
  allResults <- rbind(allResults, transform(tempRes, aID=file_it))
  write.table(rejectDat, sprintf("%s_aID%d.txt", rejectFnameRoot, file_it),
              quote=FALSE, sep='\t', row.names=FALSE)
}

# save
write.table(allResults, fnameOut, append=F, quote=F, row.names=F, col.names=T, sep='\t')




