# CRC Pleiotropy Analysis - Simple Version
.libPaths("/home/swang25/R/ubuntu/4.4.1")

library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(cowplot)
library(xtable)
library(devtools)
library(here)

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("Fig1", "output")
dataDir <- here::here("Data")

# for colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# plot manhattan function
plotManhattan <- function(plotRes, chrCounts, colValues, shapeValues, ylimits, legName) {
  # arrange data by chromosome
  plotRes <- plotRes %>% arrange(Chr)
  uniqueChrs <- sort(unique(plotRes$Chr))
  
  # add true positions
  truePos <- rep(NA, nrow(plotRes))
  counter <- 1
  for (chr_it in 1:length(uniqueChrs)) {
    tempChr <- uniqueChrs[chr_it]
    tempDat <- plotRes %>% filter(Chr == tempChr)
    truePos[counter:(counter + nrow(tempDat) - 1)] <- rep(sum(chrCounts[1:tempChr]), nrow(tempDat)) + tempDat$BP
    counter <- counter + nrow(tempDat)
  }
  
  # plot
  xBreaks <- cumsum(chrCounts[-1])
  xBreaksLabs <- 1:22
  xBreaksLabs[c(9, 11, 13, 15, 16, 17, 19, 20, 21)] <- ""
  
  plotDat <- plotRes %>% mutate(truePos = truePos)
  
  returnPlot <- ggplot(plotDat, aes(x=truePos, y=-log10(newLfdr), color=as.factor(cat), shape=as.factor(cat))) +
    geom_point() +
    xlab("Chromosome") + ylab("-log10(lfdr)") +
    scale_color_manual(name=legName, values=colValues) +
    scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(name="Chr", breaks=xBreaks, labels=xBreaksLabs) +
    ylim(ylimits) +
    theme_cowplot() +
    theme(axis.text=element_text(size=18), axis.title=element_text(size=18),
          legend.title=element_text(size=18), legend.text=element_text(size=18)) +
    guides(colour = guide_legend(override.aes = list(size=4)))
  
  return(returnPlot)
}

# load data
aID1 <- fread(here::here(outputDir, "Fig4_data_aID_ob21_newlfdr.txt"))
aID2 <- fread(here::here(outputDir, "Fig4_data_aID2_newlfdr.txt"))

# add position information (no filtering needed since rejNew doesn't exist)
aID1_new <- aID1 %>%
  mutate(Chr = as.numeric(chr), BP = as.numeric(bp)) %>%
  mutate(chrpos = paste(Chr, BP, sep = ":")) %>%
  mutate(newLfdr = lfdr_new)

aID2_new <- aID2 %>%
  mutate(Chr = as.numeric(chr), BP = as.numeric(bp)) %>%
  mutate(chrpos = paste(Chr, BP, sep = ":")) %>%
  mutate(newLfdr = lfdr_new)

# for plotting axes - calculate from existing data instead of bmi_with_overall.txt (did not create it yet)
allZ <- rbind(
  aID1_new %>% select(Chr, BP),
  aID2_new %>% select(Chr, BP)
) %>% distinct()

chrCounts <- rep(0, 23)
for (chr_it in 1:22) {
  tempDat <- allZ %>% filter(Chr == chr_it)
  if (nrow(tempDat) > 0) {
    maxPos <- max(tempDat$BP)
    chrCounts[chr_it + 1] <- maxPos
  }
}

# data for manhattan plot - CRC pleiotropy
manDataTwo <- rbind(aID1_new %>% select(chrpos, Chr, BP, newLfdr) %>% mutate(cat = "CRC,Obesity"),
                    aID2_new %>% select(chrpos, Chr, BP, newLfdr) %>% mutate(cat = "CRC,IBD")) %>%
  as.data.frame(.) %>%
  mutate(Chr = as.numeric(Chr), BP = as.numeric(BP)) %>%
  arrange(newLfdr) %>% # lowest lfdr first
  #distinct(., chrpos, .keep_all = TRUE)
  distinct(chrpos, cat, .keep_all = TRUE)   # keep one row per pair--show all significants  at same locus,



manPlotTwo <- plotManhattan(plotRes = manDataTwo, chrCounts,
                            colValues = gg_color_hue(2), shapeValues=c(16,17), 
                            ylimits=c(0, 12.5), legName="CRC Pleiotropy")
manPlotTwo

ggsave(filename = paste0(outputDir, "/CRC_Pleiotropy_Manhattan.pdf"),
       plot = manPlotTwo, width = 18, height = 12)



