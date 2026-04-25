#!/usr/bin/env Rscript
# Z-score scatter plots for pleiotropy analysis
# Creates scatter plots of Z1 vs Z2 highlighting pleiotropic hits

.libPaths("/home/swang25/R/ubuntu/4.4.1")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(here)
  library(viridis)
})

create_zscore_scatter <- function(dt, trait1_name, trait2_name, 
                                  lfdr_threshold = 0.1, alpha_all = 0.3) {
  
  message("Processing ", nrow(dt), " SNPs for ", trait1_name, " vs ", trait2_name)
  
  # Classify SNPs
  dt <- dt %>%
    mutate(
      pleiotropic = lfdr_new < lfdr_threshold,
      quadrant = case_when(
        Z1 > 0 & Z2 > 0 ~ "Q1: +/+",
        Z1 < 0 & Z2 > 0 ~ "Q2: -/+", 
        Z1 < 0 & Z2 < 0 ~ "Q3: -/-",
        Z1 > 0 & Z2 < 0 ~ "Q4: +/-",
        TRUE ~ "Origin"
      ),
      point_type = ifelse(pleiotropic, "Pleiotropic", "Non-pleiotropic")
    )
  
  # Summary statistics
  n_pleio <- sum(dt$pleiotropic)
  correlation <- cor(dt$Z1, dt$Z2, use = "complete.obs")
  
  message("  Pleiotropic SNPs: ", n_pleio, " (", round(100*n_pleio/nrow(dt), 2), "%)")
  message("  Z-score correlation: ", round(correlation, 3))
  
  # Quadrant counts for pleiotropic SNPs
  quadrant_counts <- dt %>%
    filter(pleiotropic) %>%
    count(quadrant) %>%
    arrange(desc(n))
  
  message("  Pleiotropic SNP distribution by quadrant:")
  print(quadrant_counts)
  
  # Create the plot
  p <- ggplot(dt, aes(x = Z1, y = Z2)) +
    # Background points (non-pleiotropic)
    geom_point(data = filter(dt, !pleiotropic), 
               alpha = alpha_all, size = 0.5, color = "grey70") +
    # Pleiotropic points
    geom_point(data = filter(dt, pleiotropic),
               aes(color = -log10(lfdr_new)), 
               size = 1.2, alpha = 0.8) +
    # Color scale for pleiotropic points
    scale_color_viridis_c(name = "-log10(lfdr)", 
                          option = "plasma", 
                          trans = "sqrt") +
    # Add quadrant lines
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    # Labels and theme
    labs(
      x = paste0("Z-score ", trait1_name),
      y = paste0("Z-score ", trait2_name),
      title = paste0("Z-score Correlation: ", trait1_name, " vs ", trait2_name),
      subtitle = paste0("n = ", format(nrow(dt), big.mark = ","), " SNPs | ", 
                        n_pleio, " pleiotropic (lfdr < ", lfdr_threshold, ") | ",
                        "r = ", round(correlation, 3)),
      caption = "Grey points: non-pleiotropic SNPs\nColored points: pleiotropic SNPs (colored by -log10(lfdr))"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      plot.caption = element_text(size = 10, hjust = 0),
      axis.title = element_text(size = 12),
      legend.position = "right"
    ) +
    # Equal scaling for both axes to show true correlation
    coord_fixed(ratio = 1, xlim = c(-8.5, 8.5), ylim = c(-8.5, 8.5))
  
  # Add quadrant labels
  p <- p + 
    annotate("text", x = 6, y = 6, label = "Q1: +/+", 
             size = 3, alpha = 0.7, fontface = "italic") +
    annotate("text", x = -6, y = 6, label = "Q2: -/+", 
             size = 3, alpha = 0.7, fontface = "italic") +
    annotate("text", x = -6, y = -6, label = "Q3: -/-", 
             size = 3, alpha = 0.7, fontface = "italic") +
    annotate("text", x = 6, y = -6, label = "Q4: +/-", 
             size = 3, alpha = 0.7, fontface = "italic")
  
  return(list(plot = p, data_summary = list(
    n_total = nrow(dt),
    n_pleiotropic = n_pleio,
    correlation = correlation,
    quadrant_counts = quadrant_counts
  )))
}

# Function to create density overlay
add_density_overlay <- function(base_plot, dt) {
  dt_pleio <- filter(dt, lfdr_new < 0.01)
  if(nrow(dt_pleio) > 10) {  # Need enough points for density
    base_plot + 
      stat_density_2d(data = dt_pleio, 
                      aes(x = Z1, y = Z2), 
                      color = "red", alpha = 0.6, linewidth = 0.5)
  } else {
    message("Not enough pleiotropic SNPs for density overlay")
    return(base_plot)
  }
}

# Set up file paths
outdir <- here("Fig1", "output")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# File paths to the actual results files (not raw data)
crc_obesity_file <- file.path(outdir, "Fig4_data_aID1_newlfdr.txt")
crc_ibd_file <- file.path(outdir, "Fig4_data_aID2_newlfdr.txt")

# Check if files exist and process
message("CRC-Obesity file: ", crc_obesity_file)
message("CRC-IBD file: ", crc_ibd_file)

# Create plots for CRC-Obesity
if(file.exists(crc_obesity_file)) {
  message("Creating CRC-Obesity Z-score scatter plot...")
  crc_obesity_data <- fread(crc_obesity_file, na.strings = c("NA", "#NA", ""))
  
  # Check if required columns exist
  required_cols <- c("Z1", "Z2", "lfdr_new")
  if(all(required_cols %in% names(crc_obesity_data))) {
    obesity_result <- create_zscore_scatter(
      dt = crc_obesity_data,
      trait1_name = "CRC",
      trait2_name = "Obesity",
      lfdr_threshold = 0.1
    )
    
    # Save plot
    ggsave(file.path(outdir, "CRC_Obesity_Zscore_Scatter.png"), 
           obesity_result$plot, 
           width = 10, height = 8, dpi = 300)
    
    # Save with density overlay
    obesity_density <- add_density_overlay(obesity_result$plot, crc_obesity_data)
    ggsave(file.path(outdir, "CRC_Obesity_Zscore_Scatter_Density.png"), 
           obesity_density, 
           width = 10, height = 8, dpi = 300)
  } else {
    message("Missing required columns in CRC-Obesity file. Found: ", paste(names(crc_obesity_data), collapse = ", "))
  }
} else {
  message("CRC-Obesity results file not found: ", crc_obesity_file)
}

# Create plots for CRC-IBD
if(file.exists(crc_ibd_file)) {
  message("Creating CRC-IBD Z-score scatter plot...")
  crc_ibd_data <- fread(crc_ibd_file, na.strings = c("NA", "#NA", ""))
  
  # Check if required columns exist
  required_cols <- c("Z1", "Z2", "lfdr_new")
  if(all(required_cols %in% names(crc_ibd_data))) {
    ibd_result <- create_zscore_scatter(
      dt = crc_ibd_data,
      trait1_name = "CRC",
      trait2_name = "IBD", 
      lfdr_threshold = 0.01
    )
    
    # Save plot
    ggsave(file.path(outdir, "CRC_IBD_Zscore_Scatter.png"), 
           ibd_result$plot, 
           width = 10, height = 8, dpi = 300)
    
    # Save with density overlay  
    ibd_density <- add_density_overlay(ibd_result$plot, crc_ibd_data)
    ggsave(file.path(outdir, "CRC_IBD_Zscore_Scatter_Density.png"), 
           ibd_density, 
           width = 10, height = 8, dpi = 300)
  } else {
    message("Missing required columns in CRC-IBD file. Found: ", paste(names(crc_ibd_data), collapse = ", "))
  }
} else {
  message("CRC-IBD results file not found: ", crc_ibd_file)
}

# Create summary table
create_summary_table <- function() {
  summary_data <- data.frame(
    Trait_Pair = character(),
    Total_SNPs = integer(),
    Pleiotropic_SNPs = integer(),
    Pleiotropy_Rate = numeric(),
    Z_Correlation = numeric(),
    Q1_Plus_Plus = integer(),
    Q2_Minus_Plus = integer(), 
    Q3_Minus_Minus = integer(),
    Q4_Plus_Minus = integer()
  )
  
  if(exists("obesity_result")) {
    quad_counts <- obesity_result$data_summary$quadrant_counts
    summary_data <- rbind(summary_data, data.frame(
      Trait_Pair = "CRC-Obesity",
      Total_SNPs = obesity_result$data_summary$n_total,
      Pleiotropic_SNPs = obesity_result$data_summary$n_pleiotropic,
      Pleiotropy_Rate = round(100 * obesity_result$data_summary$n_pleiotropic / 
                                obesity_result$data_summary$n_total, 2),
      Z_Correlation = round(obesity_result$data_summary$correlation, 3),
      Q1_Plus_Plus = ifelse("Q1: +/+" %in% quad_counts$quadrant, 
                            quad_counts$n[quad_counts$quadrant == "Q1: +/+"], 0),
      Q2_Minus_Plus = ifelse("Q2: -/+" %in% quad_counts$quadrant,
                             quad_counts$n[quad_counts$quadrant == "Q2: -/+"], 0),
      Q3_Minus_Minus = ifelse("Q3: -/-" %in% quad_counts$quadrant,
                              quad_counts$n[quad_counts$quadrant == "Q3: -/-"], 0),
      Q4_Plus_Minus = ifelse("Q4: +/-" %in% quad_counts$quadrant,
                             quad_counts$n[quad_counts$quadrant == "Q4: +/-"], 0)
    ))
  }
  
  if(exists("ibd_result")) {
    quad_counts <- ibd_result$data_summary$quadrant_counts
    summary_data <- rbind(summary_data, data.frame(
      Trait_Pair = "CRC-IBD",
      Total_SNPs = ibd_result$data_summary$n_total,
      Pleiotropic_SNPs = ibd_result$data_summary$n_pleiotropic,
      Pleiotropy_Rate = round(100 * ibd_result$data_summary$n_pleiotropic / 
                                ibd_result$data_summary$n_total, 2),
      Z_Correlation = round(ibd_result$data_summary$correlation, 3),
      Q1_Plus_Plus = ifelse("Q1: +/+" %in% quad_counts$quadrant,
                            quad_counts$n[quad_counts$quadrant == "Q1: +/+"], 0),
      Q2_Minus_Plus = ifelse("Q2: -/+" %in% quad_counts$quadrant,
                             quad_counts$n[quad_counts$quadrant == "Q2: -/+"], 0),
      Q3_Minus_Minus = ifelse("Q3: -/-" %in% quad_counts$quadrant,
                              quad_counts$n[quad_counts$quadrant == "Q3: -/-"], 0),
      Q4_Plus_Minus = ifelse("Q4: +/-" %in% quad_counts$quadrant,
                             quad_counts$n[quad_counts$quadrant == "Q4: +/-"], 0)
    ))
  }
  
  if(nrow(summary_data) > 0) {
    write.table(summary_data, 
                file.path(outdir, "Zscore_Scatter_Summary.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  return(summary_data)
}

# Generate summary
summary_table <- create_summary_table()
if(nrow(summary_table) > 0) {
  message("\n=== PLEIOTROPY SUMMARY ===")
  print(summary_table)
} else {
  message("No data processed - check file paths and column names")
}

message("\n=== Z-SCORE SCATTER ANALYSIS COMPLETE ===")
message("Output directory: ", outdir)
message("Check for the following files:")
message("- CRC_Obesity_Zscore_Scatter.png")  
message("- CRC_IBD_Zscore_Scatter.png")
message("- CRC_Obesity_Zscore_Scatter_Density.png")
message("- CRC_IBD_Zscore_Scatter_Density.png")
message("- Zscore_Scatter_Summary.txt")