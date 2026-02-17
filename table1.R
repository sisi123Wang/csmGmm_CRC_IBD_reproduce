.libPaths("/home/swang25/R/ubuntu/4.4.1")

# Make Table 1
# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig4/plot_data_analysis.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
#setwd("/home/swang25/other_summary_statistics/CRC_IBD/")
here::i_am("CRC_IBD//table1.R")

library(utils)
library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(cowplot)
library(data.table)
library(xtable)
library(devtools)
library(wesanderson)
library(here)
library(manhattan)
library(ggpubr)

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("CRC_IBD", "output")
outputDir_Fig <- here::here("CRC_IBD", "figures")
dataDir <- here::here("Data")

#----------------------------------------------------------------------------------------#
#Table for simulation
final_result <- data.frame(nCausal =NA, propCausal = NA , nRejNew = NA, nRejMeta = NA, nRejThreshold1 = NA, nRejThreshold2 = NA, nRejThreshold3 = NA,
                           powerNew = NA, powerMeta = NA, powerThreshold1 = NA, powerThreshold2 = NA, powerThreshold3 = NA,
                           fdpNew = NA, fdpMeta = NA,fdpThreshold1 = NA,fdpThreshold2 = NA,fdpThreshold3 = NA)

sim_results <- fread(paste0(outputDir, "/pleio_approach_aID", 1, ".txt")) %>%
  select("nCausal","propCausal", "nRejNew", "nRejMeta", "nRejThreshold1", "nRejThreshold2", "nRejThreshold3", 
         "powerNew","powerMeta","powerThreshold1", "powerThreshold2","powerThreshold3",
         "fdpNew", "fdpMeta" ,"fdpThreshold1" ,"fdpThreshold2","fdpThreshold3")

for (aID in 2:100) {
  tempResults <- fread(paste0(outputDir, "/pleio_approach_aID", aID, ".txt")) %>%
    select("nCausal","propCausal", "nRejNew", "nRejMeta", "nRejThreshold1", "nRejThreshold2", "nRejThreshold3", 
           "powerNew","powerMeta","powerThreshold1", "powerThreshold2","powerThreshold3",
           "fdpNew", "fdpMeta" ,"fdpThreshold1" ,"fdpThreshold2","fdpThreshold3")
  sim_results <- rbind(sim_results,tempResults)
}

targets <- c(4e-04, 2e-03, 4e-03)   # s1, s2, s3
tol <- 5e-05                         # tolerance range (±0.00005)

pCausal1 <- subset(sim_results, abs(propCausal - targets[1]) < tol)
pCausal2 <- subset(sim_results, abs(propCausal - targets[2]) < tol)
pCausal3 <- subset(sim_results, abs(propCausal - targets[3]) < tol)

pCausal1 <- sapply(pCausal1, mean, na.rm = TRUE)
pCausal2 <- sapply(pCausal2, mean, na.rm = TRUE)
pCausal3 <- sapply(pCausal3, mean, na.rm = TRUE)

final_result[1,] = pCausal1
final_result[2,] = pCausal2
final_result[3,] = pCausal3
write.table(final_result, paste0(outputDir_Fig, "/Tab1.txt"), append=F, quote=F, row.names=F, col.names=T)

setwd("/home/swang25/other_summary_statistics/CRC_IBD/figures/")
Tab1 <- read.table("Tab1.txt", header=TRUE)
View(Tab1)
