.libPaths(c("/home/swang25/R/ubuntu/4.4.1", Sys.getenv("R_LIBS_USER"), .libPaths()))
library(ggplot2)
library(dplyr)
library(data.table)
library(readr)
library(cowplot)

prepare_manhattan_data <- function(gwas_data, trait_name) {
  gwas_data %>%
    filter(!is.na(chromosome), !is.na(base_pair_location), !is.na(p_value)) %>%
    mutate(
      Chr = as.numeric(chromosome),
      BP = as.numeric(base_pair_location),
      P = as.numeric(p_value),
      chrpos = paste(Chr, BP, sep = ":"),
      logP = -log10(P),
      trait = trait_name
    ) %>%
    filter(P > 0, P <= 1, is.finite(logP)) %>%
    # Sort by significance
    arrange(P) %>%
    # Remove duplicates, keep most significant per position
    distinct(chrpos, .keep_all = TRUE) %>%
    # Select relevant columns
    select(chrpos, Chr, BP, P, logP, trait)
}

# OBESITY
crc_ob_harmonized <- fread("Data/crc_obesity_reharmonized.txt", na.strings = c("NA", "#NA"))
cat("Original Obesity data rows:", nrow(crc_ob_harmonized), "\n")

ob_manhattan_data <- crc_ob_harmonized %>%
  select(chromosome, base_pair_location, p_ob) %>%
  rename(p_value = p_ob) %>%
  prepare_manhattan_data(., "Obesity")

cat("Obesity Manhattan data:", nrow(ob_manhattan_data), "SNPs\n")

#chrCounts
chrCounts <- rep(0, 23)
for (chr_it in 1:22) {
  tempDat <- ob_manhattan_data %>% filter(Chr == chr_it)
  if (nrow(tempDat) > 0) {
    maxPos <- max(tempDat$BP)
    chrCounts[chr_it + 1] <- maxPos
  }
}

plotManhattan_GWAS <- function(plotRes, chrCounts, colValues, shapeValues, ylimits, legName) {
  plotRes <- plotRes %>% arrange(Chr)
  uniqueChrs <- sort(unique(plotRes$Chr))
  
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
  
  returnPlot <- ggplot(plotDat, aes(x=truePos, y=logP, color=trait)) +
    geom_point(size=0.5) +
    xlab("Chromosome") + ylab("-log10(P)") +
    scale_color_manual(name=legName, values=colValues) +
    scale_x_continuous(name="Chr", breaks=xBreaks, labels=xBreaksLabs) +
    ylim(ylimits) +
    theme_cowplot() +
    theme(axis.text=element_text(size=18), axis.title=element_text(size=18),
          legend.title=element_text(size=18), legend.text=element_text(size=18)) +
    guides(colour = guide_legend(override.aes = list(size=4)))
  
  return(returnPlot)
}

# Obesity Manhattan Plot
manPlot_Obesity <- plotManhattan_GWAS(
  plotRes = ob_manhattan_data,
  chrCounts = chrCounts, 
  colValues = "#FF7F00",      # Orange color for Obesity
  shapeValues = 17,           # Triangle shape
  ylimits = c(0, max(ob_manhattan_data$logP, na.rm = TRUE) + 1),
  legName = "Obesity GWAS"
)

outputDir <- "output"
dir.create(outputDir, showWarnings = FALSE)

ggsave(
  filename = paste0(outputDir, "/Obesity_GWAS_Manhattan.pdf"),
  plot = manPlot_Obesity,
  width = 16, height = 10
)


# Obesity Statistics
obesity_genome_wide_sig <- sum(ob_manhattan_data$P < 5e-8, na.rm = TRUE)
obesity_suggestive_sig <- sum(ob_manhattan_data$P < 1e-5, na.rm = TRUE)

cat("Obesity GWAS:\n")
cat("  Total SNPs:", nrow(ob_manhattan_data), "\n")
cat("  Genome-wide significant (p < 5e-8):", obesity_genome_wide_sig, "\n")
cat("  Suggestive significant (p < 1e-5):", obesity_suggestive_sig, "\n")
cat("  Most significant p-value:", min(ob_manhattan_data$P, na.rm = TRUE), "\n\n")

# Top 10 Obesity hits
top_obesity_hits <- ob_manhattan_data %>%
  arrange(P) %>%
  head(10) %>%
  select(Chr, BP, P, logP)

cat("TOP 10 OBESITY HITS")
print(top_obesity_hits)



