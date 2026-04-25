library(ggplot2)
library(dplyr)
library(data.table)
library(readr)
library(cowplot)
library(gridExtra)
.libPaths("/home/swang25/R/ubuntu/4.4.1")

prepare_manhattan_data <- function(gwas_data, trait_name) {
  gwas_data %>%
    # Remove any rows with missing essential columns
    filter(!is.na(chromosome), !is.na(base_pair_location), !is.na(p_value)) %>%
    # Ensure proper data types
    mutate(
      Chr = as.numeric(chromosome),
      BP = as.numeric(base_pair_location),
      P = as.numeric(p_value),
      # Create position identifier - remove duplicates properly
      chrpos = paste(Chr, BP, sep = ":"),
      # Calculate -log10(p)
      logP = -log10(P),
      # Add trait label
      trait = trait_name
    ) %>%
    # Remove invalid p-values
    filter(P > 0, P <= 1, is.finite(logP)) %>%
    # Sort by significance
    arrange(P) %>%
    # Remove duplicates, keep most significant per position
    distinct(chrpos, .keep_all = TRUE) %>%
    # Select relevant columns
    select(chrpos, Chr, BP, P, logP, trait)
}

# IBD
crc_ibd_harmonized <- fread("Data/crc_ibd.txt", na.strings = c("NA", "#NA"))
cat("Original IBD data rows:", nrow(crc_ibd_harmonized), "\n")

ibd_manhattan_data <- crc_ibd_harmonized %>%
  select(chromosome, base_pair_location, p_ibd) %>%
  rename(p_value = p_ibd) %>%
  prepare_manhattan_data(., "IBD")

cat("IBD Manhattan data:", nrow(ibd_manhattan_data), "SNPs\n")

# OBESITY  
crc_ob_harmonized <- fread("Data/crc_obesity_reharmonized.txt", na.strings = c("NA", "#NA"))
cat("Original Obesity data rows:", nrow(crc_ob_harmonized), "\n")

ob_manhattan_data <- crc_ob_harmonized %>%
  select(chromosome, base_pair_location, p_ob) %>%
  rename(p_value = p_ob) %>%
  prepare_manhattan_data(., "Obesity")

cat("Obesity Manhattan data:", nrow(ob_manhattan_data), "SNPs\n")



# chromosome coverage
allZ <- rbind(
  ibd_manhattan_data %>% select(Chr, BP),
  ob_manhattan_data %>% select(Chr, BP)
) %>% distinct()

# chrCounts
chrCounts <- rep(0, 23)
for (chr_it in 1:22) {
  tempDat <- allZ %>% filter(Chr == chr_it)
  if (nrow(tempDat) > 0) {
    maxPos <- max(tempDat$BP)
    chrCounts[chr_it + 1] <- maxPos
  }
}


plotManhattan_GWAS <- function(plotRes, chrCounts, colValues, shapeValues, ylimits, legName) {
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

# IBD Manhattan Plot
manPlot_IBD <- plotManhattan_GWAS(
  plotRes = ibd_manhattan_data,  # Use consistent data object
  chrCounts = chrCounts,
  colValues = "#E31A1C",      # Red color for IBD
  shapeValues = 16,           # Circle shape
  ylimits = c(0, max(ibd_manhattan_data$logP, na.rm = TRUE) + 1),  # Consistent object
  legName = "IBD GWAS"
)

# Obesity Manhattan Plot
manPlot_Obesity <- plotManhattan_GWAS(
  plotRes = ob_manhattan_data,  # Use consistent data object
  chrCounts = chrCounts, 
  colValues = "#FF7F00",      # Orange color for Obesity
  shapeValues = 17,           # Triangle shape
  ylimits = c(0, max(ob_manhattan_data$logP, na.rm = TRUE) + 1),  # Consistent object
  legName = "Obesity GWAS"
)

outputDir <- "output"
dir.create(outputDir, showWarnings = FALSE)

ggsave(
  filename = paste0(outputDir, "/IBD_GWAS_Manhattan.pdf"),
  plot = manPlot_IBD,
  width = 18, height = 8
)

ggsave(
  filename = paste0(outputDir, "/Obesity_GWAS_Manhattan.pdf"),
  plot = manPlot_Obesity,
  width = 18, height = 8
)

# GWAS SUMMARY STATISTICS

# IBD Statistics
ibd_genome_wide_sig <- sum(ibd_manhattan_data$P < 5e-8, na.rm = TRUE)
ibd_suggestive_sig <- sum(ibd_manhattan_data$P < 1e-5, na.rm = TRUE)

cat("IBD GWAS:\n")
cat("  Total SNPs:", nrow(ibd_manhattan_data), "\n")
cat("  Genome-wide significant (p < 5e-8):", ibd_genome_wide_sig, "\n")
cat("  Suggestive significant (p < 1e-5):", ibd_suggestive_sig, "\n")
cat("  Most significant p-value:", min(ibd_manhattan_data$P, na.rm = TRUE), "\n\n")

# Obesity Statistics
obesity_genome_wide_sig <- sum(ob_manhattan_data$P < 5e-8, na.rm = TRUE)
obesity_suggestive_sig <- sum(ob_manhattan_data$P < 1e-5, na.rm = TRUE)

cat("Obesity GWAS:\n")
cat("  Total SNPs:", nrow(ob_manhattan_data), "\n")
cat("  Genome-wide significant (p < 5e-8):", obesity_genome_wide_sig, "\n")
cat("  Suggestive significant (p < 1e-5):", obesity_suggestive_sig, "\n")
cat("  Most significant p-value:", min(ob_manhattan_data$P, na.rm = TRUE), "\n\n")

# ============================================================================
# TOP HITS FOR EACH TRAIT
# ============================================================================

# Top 10 IBD hits
top_ibd_hits <- ibd_manhattan_data %>%
  arrange(P) %>%
  head(10) %>%
  select(Chr, BP, P, logP)

# Top 10 Obesity hits
top_obesity_hits <- ob_manhattan_data %>%
  arrange(P) %>%
  head(10) %>%
  select(Chr, BP, P, logP)

cat("=== TOP 10 HITS ===\n")
cat("IBD Top Hits:\n")
print(top_ibd_hits)

cat("\nObesity Top Hits:\n")
print(top_obesity_hits)

# ENHANCED PLOTS WITH GWAS SIGNIFICANCE LINES

# Function to add GWAS significance lines
add_gwas_significance <- function(plot_obj) {
  plot_obj + 
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red", alpha = 0.8) +
    geom_hline(yintercept = -log10(1e-5), linetype = "dashed", color = "blue", alpha = 0.8) +
    annotate("text", x = Inf, y = -log10(5e-8), label = "Genome-wide sig.", 
             hjust = 1.1, vjust = -0.5, color = "red", size = 3) +
    annotate("text", x = Inf, y = -log10(1e-5), label = "Suggestive sig.", 
             hjust = 1.1, vjust = -0.5, color = "blue", size = 3)
}

# Enhanced plots
manPlot_IBD_Enhanced <- add_gwas_significance(manPlot_IBD)
manPlot_Obesity_Enhanced <- add_gwas_significance(manPlot_Obesity)

# Save enhanced plots
ggsave(
  filename = paste0(outputDir, "/IBD_GWAS_Manhattan_Enhanced.pdf"),
  plot = manPlot_IBD_Enhanced,
  width = 18, height = 8
)

ggsave(
  filename = paste0(outputDir, "/Obesity_GWAS_Manhattan_Enhanced.pdf"),
  plot = manPlot_Obesity_Enhanced,
  width = 18, height = 8
)

cat("\n=== PLOTS SAVED ===\n")
cat("Files created:\n")
cat("- IBD_GWAS_Manhattan.pdf\n")
cat("- Obesity_GWAS_Manhattan.pdf\n") 
cat("- IBD_Obesity_GWAS_Comparison.pdf\n")
cat("- IBD_GWAS_Manhattan_Enhanced.pdf\n")
cat("- Obesity_GWAS_Manhattan_Enhanced.pdf\n")