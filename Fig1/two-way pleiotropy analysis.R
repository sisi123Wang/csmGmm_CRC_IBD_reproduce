# Process raw analysis of UKB data

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig4/summarize_ukb_analysis.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
.libPaths("/home/swang25/R/ubuntu/4.4.1")
#setwd("/home/swang25/other_summary_statistics") 

here::i_am("Fig1/two-way pleiotropy analysis.R")

# analysis catalogue keyed by aID 
cfg <- list(
  `1` = list(pair = "crc_ob",  type = "pleio",
             zcols = c("Zcrc","Zob")),
  `2` = list(pair = "crc_ibd", type = "pleio",
             zcols = c("Zcrc","Zibd")),
  `3` = list(pair = "crc_ob",  type = "corr",
             zcols = c("Zcrc","Zob")),
  `4` = list(pair = "crc_ibd", type = "corr",
             zcols = c("Zcrc","Zibd"))
)




# load
crc_ob_harmonized <- readr::read_tsv("Data/crc_obesity.txt",
                                     show_col_types = FALSE)

crc_ibd_harmonized <- readr::read_tsv("Data/crc_ibd.txt",
                                     show_col_types = FALSE)
cor_crc_ibd <- readr::read_tsv("Data/cor_crc_ibd.txt",
                                     show_col_types = FALSE)

cor_crc_ob <- readr::read_tsv("Data/cor_crc_obesity.txt",
                                      show_col_types = FALSE)
library(purrr)
library(readr)
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
outputDir <- here::here("Fig1", "output")
dataDir <- here::here("Data")
fnameOut <- paste0(outputDir, "/processed_ukb_data_S", Snum, ".txt")
rejectFnameRoot <- paste0(outputDir, "/reject_bmi_with_overall_neg5_reject_S_", Snum)

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

# how the raw output files are named
fnameRoot <- paste0(outputDir, "/Fig4_data_aID", 1:9)
fnameDACT <- paste0(fnameRoot, "_DACTp.txt")
fnameHDMT <- paste0(fnameRoot, "_hdmt.txt")
fnameKernel <- paste0(fnameRoot, "_kernel.txt")
fname7 <- paste0(fnameRoot, "_df7.txt")
#fname50 <- paste0(fnameRoot, "_df50.txt")
fnameNew <- c(paste0(fnameRoot, "_newlfdr.txt"))

selections <- list()
selections[[1]] <- c("Zcad", "Zbmi")
selections[[2]] <- c("Zoverall", "Zcad")
selections[[3]] <- c("Zoverall", "Zbmi")
selections[[4]] <- c("Zoverall", "Zlcukb")
selections[[5]] <- c("Zcad_cardio", "Zcadukb")
selections[[6]] <- c("Zoverall", "Zcad", "Zbmi")
selections[[7]] <- c("Zlcukb", "Zcadukb")
selections[[8]] <- c("Zlcukb", "Zbmi")
selections[[9]] <- c("Zcadukb", "Zbmi")


# results
allResults <- c()
for (file_it in 1:9) {
  
  if (file_it == 4) {
    cleanZ <- fread(here::here(dataDir, "replication_with_lcoverall.txt"))
    fdrLimitHDMT <- fdrLimitHDMTr
    fdrLimitDACT <- fdrLimitDACTr
    fdrLimitKernel <- fdrLimitKernelr
    fdrLimit50 <- fdrLimit50r
    fdrLimit7 <- fdrLimit7r
  } else if (file_it == 5) {
    cleanZ <- fread(here::here(dataDir, "cad_for_replication.txt"))
    fdrLimitHDMT <- fdrLimitHDMTr
    fdrLimitDACT <- fdrLimitDACTr
    fdrLimitKernel <- fdrLimitKernelr
    fdrLimit50 <- fdrLimit50r
    fdrLimit7 <- fdrLimit7r
  } else {
    cleanZ <- fread(here::here(dataDir, "bmi_with_overall.txt"))
    fdrLimitHDMT <- fdrLimitHDMTi
    fdrLimitDACT <- fdrLimitDACTi
    fdrLimitKernel <- fdrLimitKerneli
    fdrLimit50 <- fdrLimit50i
    fdrLimit7 <- fdrLimit7i  
  }
  
  # hold temporary results
  tempRes <- data.frame(Method=c("New", "Kernel", "df50", "df7", "HDMT", "DACT"), numReject=NA)
  
  # new 
  tempNew <- fread(fnameNew[file_it], header=T, data.table=F)
  tempDat <- cleanZ %>% select(all_of(selections[[file_it]]), chrpos) %>%
    mutate(origIdx = 1:nrow(.)) %>%   
    mutate(newLfdr = tempNew$x) %>%
    arrange(newLfdr) %>%
    mutate(cumNew = cummean(newLfdr)) %>%
    mutate(rejNew = ifelse(cumNew < fdrLimitNew, 1, 0)) %>%
    arrange(origIdx)
  tempRes$numReject[1] <- sum(tempDat$rejNew)
  
  if (file_it >= 7) {
    # any rejection
    rejectDat <- tempDat %>% filter(rejNew == 1)
    
    # allResults
    allResults <- rbind(allResults, tempRes %>% mutate(aID = file_it))
    
    # save
    write.table(rejectDat, paste0(rejectFnameRoot, "_aID", file_it, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
    next
  }
  
  # kernel
  tempKernel <- fread(fnameKernel[file_it], header=T, data.table=F)
  tempDat <- tempDat %>% mutate(kernelLfdr = tempKernel$x) %>%
    arrange(kernelLfdr) %>%
    mutate(cumKernel = cummean(kernelLfdr)) %>%
    mutate(rejKernel = ifelse(cumKernel < fdrLimitKernel, 1, 0)) %>%
    arrange(origIdx)
  tempRes$numReject[2] <- sum(tempDat$rejKernel)
  
  # df50
  tempdf50 <- fread(fname50[file_it], header=T, data.table=F)
  tempDat <- tempDat %>% mutate(df50Lfdr = tempdf50$x) %>%
    arrange(df50Lfdr) %>%
    mutate(cumdf50 = cummean(df50Lfdr)) %>%
    mutate(rejdf50 = ifelse(cumdf50 < fdrLimit50, 1, 0)) %>%
    arrange(origIdx)
  tempRes$numReject[3] <- sum(tempDat$rejdf50)
  
  # df7
  tempdf7 <- fread(fname7[file_it], header=T, data.table=F)
  tempDat <- tempDat %>% mutate(df7Lfdr = tempdf7$x) %>%
    arrange(df7Lfdr) %>%
    mutate(cumdf7 = cummean(df7Lfdr)) %>%
    mutate(rejdf7 = ifelse(cumdf7 < fdrLimit7, 1, 0)) %>%
    arrange(origIdx)
  tempRes$numReject[4] <- sum(tempDat$rejdf7)
  
  if (file_it >= 4 & file_it <= 6) {
    # any rejection
    tempDat <- tempDat %>% mutate(rejAny = ifelse(rejdf7 == 1 |  rejdf50 == 1 | rejKernel == 1 | rejNew == 1, 1, 0))
    rejectDat <- tempDat %>% filter(rejAny == 1)
    
    # allResults
    allResults <- rbind(allResults, tempRes %>% mutate(aID = file_it))
    
    # save
    write.table(rejectDat, paste0(rejectFnameRoot, "_aID", file_it, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
    next
  }
  
  # HDMT
  tempHDMT <- fread(fnameHDMT[file_it], header=T, data.table=F)
  tempDat <- tempDat %>% mutate(hdmtFdr = tempHDMT$fixedFdr) %>%
    mutate(rejHDMT = ifelse(hdmtFdr < fdrLimitHDMT, 1, 0)) %>%
    arrange(origIdx)
  tempRes$numReject[5] <- sum(tempDat$rejHDMT)
  
  # DACT
  tempDACT <- fread(fnameDACT[file_it], header=T, data.table=F)
  tempDat <- tempDat %>% mutate(DACTp = tempDACT$x) %>%
    arrange(DACTp) %>%
    mutate(rankedIdxP = 1:nrow(.)) %>%
    mutate(km = 1:nrow(.) / nrow(.)) %>%
    mutate(RHS = km * fdrLimitDACT)
  rejected <- which(tempDat$DACTp <= tempDat$RHS)
  if (length(rejected) == 0) {
    maxIdx <- 0
  } else {maxIdx <- max(rejected)}
  tempDat <- tempDat %>% mutate(rejDACT = ifelse(rankedIdxP <= maxIdx, 1, 0)) %>%
    arrange(origIdx)
  tempRes$numReject[6] <- sum(tempDat$rejDACT)
  
  # any rejection
  tempDat <- tempDat %>% mutate(rejAny = ifelse(rejDACT == 1 | rejHDMT == 1 | rejdf7 == 1 | 
                                                  rejdf50 == 1 | rejKernel == 1 | rejNew == 1, 1, 0))
  
  rejectDat <- tempDat %>% filter(rejAny == 1) 
  
  # allResults
  allResults <- rbind(allResults, tempRes %>% mutate(aID = file_it))
  
  # save rejections
  write.table(rejectDat, paste0(rejectFnameRoot, "_aID", file_it, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
  
  cat(file_it)
}

# save
write.table(allResults, fnameOut, append=F, quote=F, row.names=F, col.names=T, sep='\t')






