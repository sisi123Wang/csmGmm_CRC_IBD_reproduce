#!/usr/bin/env Rscript
# Regional analysis for pleiotropy signals

# - Panel 1: Pleiotropy signal (-log10 lfdr)
# - Panel 2: CRC Z-scores  
# - Panel 3: Obesity/IBD Z-scores
.libPaths("/home/swang25/R/ubuntu/4.4.1")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(here)
  library(viridis)
  library(gridExtra)
})

# Function to identify top pleiotropic regions
identify_top_regions <- function(dt, trait_name, lfdr_threshold = 0.001, min_distance = 1e6) {
  # significant SNPs by lfdr below 0.001 here
  sig_snps <- dt[lfdr_new < lfdr_threshold][order(lfdr_new)]
  
  if(nrow(sig_snps) == 0) {
    message("No significant SNPs found for ", trait_name, " at threshold ", lfdr_threshold)
    return(NULL)
  }
  
  message("Found ", nrow(sig_snps), " significant SNPs for ", trait_name)
  
  # Group nearby SNPs into regions
  regions <- sig_snps %>%
    arrange(chr, bp) %>% #Sort significant snps by chromosome & position
    group_by(chr) %>%
    mutate(
      # Create region groups based on distance
      region_group = cumsum(c(1, diff(bp) > min_distance)) # Cluster nearby hits into one locus, count as next for gap exceeds min_distance
    ) %>%
    group_by(chr, region_group) %>%
    summarise(
      start_bp = min(bp),
      end_bp = max(bp),
      center_bp = round(mean(bp)),
      n_snps = n(),# number of significant SNPs that each locus contain
      min_lfdr = min(lfdr_new),# use it for the best SNP in the region ranking
      max_neg_log_lfdr = max(-log10(lfdr_new)),
      lead_snp_bp = bp[which.min(lfdr_new)][1],
      lead_snp_lfdr = min(lfdr_new),
      .groups = 'drop'
    ) %>%
    arrange(min_lfdr) %>% # smallest lfdr is placed first!
    mutate(
      region_id = paste0("Chr", chr, "_", round(center_bp/1e6, 1), "Mb"),
      flank_start = pmax(1, start_bp - 500000),  # 500kb flanking
      flank_end = end_bp + 500000
    )
  
  return(regions)
}

# Function to create regional plot
create_regional_plot <- function(dt, chr, start_bp, end_bp, trait1_name, trait2_name, 
                                 region_id, lfdr_threshold = 0.01) {
  
  # Extract regional data
  region_data <- dt %>%
    filter(chr == !!chr, bp >= start_bp, bp <= end_bp) %>%
    mutate(
      position_mb = bp / 1e6,
      significant = lfdr_new < lfdr_threshold,
      neg_log_lfdr = -log10(lfdr_new)
    )
  
  if(nrow(region_data) == 0) {
    message("No data in region Chr", chr, ":", start_bp, "-", end_bp)
    return(NULL)
  }
  
  message("Creating regional plot for ", region_id, " with ", nrow(region_data), " SNPs")
  
  # Calculate plot limits
  x_min <- min(region_data$position_mb)
  x_max <- max(region_data$position_mb)
  
  # Plot 1: Pleiotropy signal (-log10 lfdr)
  p1 <- ggplot(region_data, aes(x = position_mb, y = neg_log_lfdr)) +
    geom_point(aes(color = significant), alpha = 0.7, size = 1) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red"),
                       name = paste0("lfdr < ", lfdr_threshold)) +
    geom_hline(yintercept = -log10(lfdr_threshold), linetype = "dashed", color = "red", alpha = 0.7) +
    labs(
      title = paste0("Pleiotropy Signal: ", trait1_name, " vs ", trait2_name),
      subtitle = region_id,
      x = "",  # No x-axis label for top plot
      y = "-log10(lfdr)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_blank(),  # Remove x-axis labels
      legend.position = "right"
    ) +
    xlim(x_min, x_max)
  
  # Plot 2: Individual trait 1 Z-scores
  p2 <- ggplot(region_data, aes(x = position_mb, y = abs(Z1))) +
    geom_point(aes(color = Z1 > 0), alpha = 0.7, size = 1) +
    scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red"),
                       name = paste0(trait1_name, " effect"),
                       labels = c("Negative", "Positive")) +
    geom_hline(yintercept = 3, linetype = "dashed", alpha = 0.5) +  # |Z| > 3 threshold
    labs(
      x = "",
      y = paste0("|Z-score| ", trait1_name)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      legend.position = "right"
    ) +
    xlim(x_min, x_max)
  
  # Plot 3: Individual trait 2 Z-scores  
  p3 <- ggplot(region_data, aes(x = position_mb, y = abs(Z2))) +
    geom_point(aes(color = Z2 > 0), alpha = 0.7, size = 1) +
    scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red"),
                       name = paste0(trait2_name, " effect"),
                       labels = c("Negative", "Positive")) +
    geom_hline(yintercept = 3, linetype = "dashed", alpha = 0.5) +
    labs(
      x = paste0("Position (Mb) - Chromosome ", chr),
      y = paste0("|Z-score| ", trait2_name)
    ) +
    theme_minimal() +
    theme(legend.position = "right") +
    xlim(x_min, x_max)
  
  # Combine plots
  combined_plot <- grid.arrange(p1, p2, p3, ncol = 1, heights = c(1, 0.8, 0.8))
  
  return(list(
    plot = combined_plot,
    data = region_data,
    summary = list(
      region = region_id,
      chr = chr,
      start_mb = x_min,
      end_mb = x_max,
      n_snps = nrow(region_data),
      n_significant = sum(region_data$significant),
      max_neg_log_lfdr = max(region_data$neg_log_lfdr),
      trait1_cor = cor(region_data$Z1, region_data$Z2, use = "complete.obs")
    )
  ))
}

# Function to analyze top SNPs in region
analyze_region_snps <- function(region_data, top_n = 10) {
  top_snps <- region_data %>%
    arrange(lfdr_new) %>%
    head(top_n) %>%
    select(chr, bp, Z1, Z2, lfdr_new) %>%
    mutate(
      position_mb = round(bp / 1e6, 2),
      neg_log_lfdr = round(-log10(lfdr_new), 2),
      Z1 = round(Z1, 3),
      Z2 = round(Z2, 3),
      effect_direction = case_when(
        Z1 > 0 & Z2 > 0 ~ "Same (+/+)",
        Z1 < 0 & Z2 < 0 ~ "Same (-/-)", 
        Z1 > 0 & Z2 < 0 ~ "Opposite (+/-)",
        Z1 < 0 & Z2 > 0 ~ "Opposite (-/+)",
        TRUE ~ "Mixed"
      )
    ) %>%
    select(chr, position_mb, Z1, Z2, neg_log_lfdr, effect_direction)
  
  return(top_snps)
}

# Main analysis
main_regional_analysis <- function() {
  
  outdir <- here("Fig1", "output")
  
  # Load pleiotropy results
  message("Loading pleiotropy results...")
  
  # CRC-Obesity
  obesity_file <- file.path(outdir, "Fig4_data_aID_less_strict1_newlfdr.txt")
  if(file.exists(obesity_file)) {
    dt_obesity <- fread(obesity_file)
    message("Loaded CRC-Obesity data: ", nrow(dt_obesity), " SNPs")
  } else {
    stop("CRC-Obesity file not found: ", obesity_file)
  }
  
  # CRC-IBD  
  ibd_file <- file.path(outdir, "Fig4_data_aID_less_strict2_newlfdr.txt")
  if(file.exists(ibd_file)) {
    dt_ibd <- fread(ibd_file)
    message("Loaded CRC-IBD data: ", nrow(dt_ibd), " SNPs")
  } else {
    stop("CRC-IBD file not found: ", ibd_file)
  }
  
  # Identify top regions
  message("\n=== IDENTIFYING TOP PLEIOTROPIC REGIONS ===")
  
  obesity_regions <- identify_top_regions(dt_obesity, "CRC-Obesity", lfdr_threshold = 0.001)
  ibd_regions <- identify_top_regions(dt_ibd, "CRC-IBD", lfdr_threshold = 0.001)
  #	Significant SNPs here are identified as those whose lfdr is below 0.001.
  if(!is.null(obesity_regions)) {
    message("\nTop CRC-Obesity regions:")
    print(obesity_regions[1:min(5, nrow(obesity_regions)), ])
  }
  
  if(!is.null(ibd_regions)) {
    message("\nTop CRC-IBD regions:")
    print(ibd_regions[1:min(5, nrow(ibd_regions)), ])
  }
  
  #Found 10093 significant SNPs for CRC-Obesity
  #Found 2462 significant SNPs for CRC-IBD
  
  # Create regional plots for top regions
  message("\n=== CREATING REGIONAL PLOTS ===")
  
  regional_results <- list()
  
  # Plot top CRC-Obesity regions
  if(!is.null(obesity_regions) && nrow(obesity_regions) > 0) {
    for(i in 1:min(3, nrow(obesity_regions))) {  # Top 3 regions
      region <- obesity_regions[i, ]
      
      result <- create_regional_plot(
        dt = dt_obesity,
        chr = region$chr,
        start_bp = region$flank_start,
        end_bp = region$flank_end,
        trait1_name = "CRC",
        trait2_name = "Obesity",
        region_id = region$region_id
      )
      
      if(!is.null(result)) {
        # Save plot
        ggsave(
          file.path(outdir, paste0("Regional_CRC_Obesity_", region$region_id, ".png")),
          result$plot,
          width = 12, height = 10, dpi = 300
        )
        
        # Analyze top SNPs
        top_snps <- analyze_region_snps(result$data)
        message("\nTop SNPs in ", region$region_id, ":")
        print(top_snps)
        
        regional_results[[paste0("Obesity_", region$region_id)]] <- result
      }
    }
  }
  
  # Plot top CRC-IBD regions  
  if(!is.null(ibd_regions) && nrow(ibd_regions) > 0) {
    for(i in 1:min(3, nrow(ibd_regions))) {  # Top 3 regions
      region <- ibd_regions[i, ]
      
      result <- create_regional_plot(
        dt = dt_ibd,
        chr = region$chr,
        start_bp = region$flank_start,
        end_bp = region$flank_end,
        trait1_name = "CRC",
        trait2_name = "IBD",
        region_id = region$region_id
      )
      
      if(!is.null(result)) {
        # Save plot
        ggsave(
          file.path(outdir, paste0("Regional_CRC_IBD_", region$region_id, ".png")),
          result$plot,
          width = 12, height = 10, dpi = 300
        )
        
        # Analyze top SNPs
        top_snps <- analyze_region_snps(result$data)
        message("\nTop SNPs in ", region$region_id, ":")
        print(top_snps)
        
        regional_results[[paste0("IBD_", region$region_id)]] <- result
      }
    }
  }
  
  return(regional_results)
}

# results and analysis!!
results <- main_regional_analysis()

library(data.table); library(dplyr)

dt_obesity <- fread("Fig1/output/Fig4_data_aID_less_strict1_newlfdr.txt")
dt_ibd     <- fread("Fig1/output/Fig4_data_aID_less_strict2_newlfdr.txt")

obesity_regions <- identify_top_regions(dt_obesity, "CRC-Obesity", 0.001)
ibd_regions     <- identify_top_regions(dt_ibd,     "CRC-IBD",     0.001)

#sanity check
#plot for lfdr < 0.001, the taller the points, the stronger evidence that the SNP is pleiotropic. 
library(ggplot2)
sig_both <- rbind(
  mutate(dt_obesity[lfdr_new < 0.001], pair = "CRC-Obesity"),
  mutate(dt_ibd    [lfdr_new < 0.001], pair = "CRC-IBD")
)

ggplot(sig_both, aes(x = as.numeric(chr) + bp/1e8, y = -log10(lfdr_new),
                     colour = pair)) +
  geom_point(alpha = 0.6, size = 0.4) +
  scale_x_continuous(breaks = 1:22, labels = 1:22) +
  labs(x = "Chr", y = "-log10(lfdr)", colour = NULL) +
  theme_minimal()

#The plot here did not discard duplicates, so here have more dots.




#!/usr/bin/env Rscript
# Regional analysis for pleiotropy signals - CRC vs IBD only

# - Panel 1: Pleiotropy signal (-log10 lfdr)
# - Panel 2: CRC Z-scores  
# - Panel 3: IBD Z-scores
.libPaths("/home/swang25/R/ubuntu/4.4.1")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(here)
  library(viridis)
  library(gridExtra)
})

# Function to identify top pleiotropic regions
identify_top_regions <- function(dt, trait_name, lfdr_threshold = 0.001, min_distance = 1e6) {
  # significant SNPs by lfdr below 0.001 here
  sig_snps <- dt[lfdr_new < lfdr_threshold][order(lfdr_new)]
  
  if(nrow(sig_snps) == 0) {
    message("No significant SNPs found for ", trait_name, " at threshold ", lfdr_threshold)
    return(NULL)
  }
  
  message("Found ", nrow(sig_snps), " significant SNPs for ", trait_name)
  
  # Group nearby SNPs into regions
  regions <- sig_snps %>%
    arrange(chr, bp) %>% #Sort significant snps by chromosome & position
    group_by(chr) %>%
    mutate(
      # Create region groups based on distance
      region_group = cumsum(c(1, diff(bp) > min_distance)) # Cluster nearby hits into one locus, count as next for gap exceeds min_distance
    ) %>%
    group_by(chr, region_group) %>%
    summarise(
      start_bp = min(bp),
      end_bp = max(bp),
      center_bp = round(mean(bp)),
      n_snps = n(),# number of significant SNPs that each locus contain
      min_lfdr = min(lfdr_new),# use it for the best SNP in the region ranking
      max_neg_log_lfdr = max(-log10(lfdr_new)),
      lead_snp_bp = bp[which.min(lfdr_new)][1],
      lead_snp_lfdr = min(lfdr_new),
      .groups = 'drop'
    ) %>%
    arrange(min_lfdr) %>% # smallest lfdr is placed first!
    mutate(
      region_id = paste0("Chr", chr, "_", round(center_bp/1e6, 1), "Mb"),
      flank_start = pmax(1, start_bp - 500000),  # 500kb flanking
      flank_end = end_bp + 500000
    )
  
  return(regions)
}

# Function to create regional plot with P-values and rs numbers
create_regional_plot <- function(dt, chr, start_bp, end_bp, trait1_name, trait2_name, 
                                 region_id, lfdr_threshold = 0.01) {
  
  # Extract regional data
  region_data <- dt %>%
    filter(chr == !!chr, bp >= start_bp, bp <= end_bp) %>%
    mutate(
      position_mb = bp / 1e6,
      significant = lfdr_new < lfdr_threshold,
      neg_log_lfdr = -log10(lfdr_new),
      # Calculate -log10 P-values for both traits
      neg_log_p1 = -log10(pnorm(abs(Z1), lower.tail = FALSE) * 2), # Two-tailed p-value for trait 1
      neg_log_p2 = -log10(pnorm(abs(Z2), lower.tail = FALSE) * 2)  # Two-tailed p-value for trait 2
    )
  
  if(nrow(region_data) == 0) {
    message("No data in region Chr", chr, ":", start_bp, "-", end_bp)
    return(NULL)
  }
  
  message("Creating regional plot for ", region_id, " with ", nrow(region_data), " SNPs")
  
  # Calculate plot limits
  x_min <- min(region_data$position_mb)
  x_max <- max(region_data$position_mb)
  
  # Find most significant SNPs for labeling
  top_pleiotropy_snp <- region_data %>% slice_min(lfdr_new, n = 1)
  top_trait1_snp <- region_data %>% slice_max(neg_log_p1, n = 1)
  top_trait2_snp <- region_data %>% slice_max(neg_log_p2, n = 1)
  
  # Plot 1: Pleiotropy signal (-log10 lfdr)
  p1 <- ggplot(region_data, aes(x = position_mb, y = neg_log_lfdr)) +
    geom_point(aes(color = significant), alpha = 0.7, size = 1) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red"),
                       name = paste0("lfdr < ", lfdr_threshold)) +
    geom_hline(yintercept = -log10(lfdr_threshold), linetype = "dashed", color = "red", alpha = 0.7) +
    # Add label for most significant pleiotropy SNP
    {if(nrow(top_pleiotropy_snp) > 0 && !is.na(top_pleiotropy_snp$snp[1])) {
      geom_text(data = top_pleiotropy_snp, 
                aes(x = position_mb, y = neg_log_lfdr, label = snp),
                vjust = -0.5, hjust = 0.5, size = 3, color = "black")
    }} +
    labs(
      title = paste0("Pleiotropy Signal: ", trait1_name, " vs ", trait2_name),
      subtitle = region_id,
      x = "",  # No x-axis label for top plot
      y = "-log10(lfdr)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_blank(),  # Remove x-axis labels
      legend.position = "right"
    ) +
    xlim(x_min, x_max)
  
  # Plot 2: Individual trait 1 P-values
  p2 <- ggplot(region_data, aes(x = position_mb, y = neg_log_p1)) +
    geom_point(aes(color = Z1 > 0), alpha = 0.7, size = 1) +
    scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red"),
                       name = paste0(trait1_name, " effect"),
                       labels = c("Negative", "Positive")) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", alpha = 0.5, color = "grey50") +  # GWAS significance
    # Add label for most significant trait 1 SNP
    {if(nrow(top_trait1_snp) > 0 && !is.na(top_trait1_snp$snp[1])) {
      geom_text(data = top_trait1_snp, 
                aes(x = position_mb, y = neg_log_p1, label = snp),
                vjust = -0.5, hjust = 0.5, size = 3, color = "black")
    }} +
    labs(
      x = "",
      y = paste0("-log10(P) ", trait1_name)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      legend.position = "right"
    ) +
    xlim(x_min, x_max)
  
  # Plot 3: Individual trait 2 P-values
  p3 <- ggplot(region_data, aes(x = position_mb, y = neg_log_p2)) +
    geom_point(aes(color = Z2 > 0), alpha = 0.7, size = 1) +
    scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red"),
                       name = paste0(trait2_name, " effect"),
                       labels = c("Negative", "Positive")) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", alpha = 0.5, color = "grey50") +  # GWAS significance
    # Add label for most significant trait 2 SNP
    {if(nrow(top_trait2_snp) > 0 && !is.na(top_trait2_snp$snp[1])) {
      geom_text(data = top_trait2_snp, 
                aes(x = position_mb, y = neg_log_p2, label = snp),
                vjust = -0.5, hjust = 0.5, size = 3, color = "black")
    }} +
    labs(
      x = paste0("Position (Mb) - Chromosome ", chr),
      y = paste0("-log10(P) ", trait2_name)
    ) +
    theme_minimal() +
    theme(legend.position = "right") +
    xlim(x_min, x_max)
  
  # Combine plots
  combined_plot <- grid.arrange(p1, p2, p3, ncol = 1, heights = c(1, 0.8, 0.8))
  
  return(list(
    plot = combined_plot,
    data = region_data,
    summary = list(
      region = region_id,
      chr = chr,
      start_mb = x_min,
      end_mb = x_max,
      n_snps = nrow(region_data),
      n_significant = sum(region_data$significant),
      max_neg_log_lfdr = max(region_data$neg_log_lfdr),
      trait1_cor = cor(region_data$Z1, region_data$Z2, use = "complete.obs"),
      top_pleiotropy_snp = if(nrow(top_pleiotropy_snp) > 0) top_pleiotropy_snp$snp[1] else NA,
      top_trait1_snp = if(nrow(top_trait1_snp) > 0) top_trait1_snp$snp[1] else NA,
      top_trait2_snp = if(nrow(top_trait2_snp) > 0) top_trait2_snp$snp[1] else NA
    )
  ))
}


# Function to analyze top SNPs in region
analyze_region_snps <- function(region_data, top_n = 10) {
  top_snps <- region_data %>%
    arrange(lfdr_new) %>%
    head(top_n) %>%
    select(chr, bp, snp, Z1, Z2, lfdr_new, neg_log_p1, neg_log_p2) %>%  # Changed rsid to snp
    mutate(
      position_mb = round(bp / 1e6, 2),
      neg_log_lfdr = round(-log10(lfdr_new), 2),
      Z1 = round(Z1, 3),
      Z2 = round(Z2, 3),
      neg_log_p1 = round(neg_log_p1, 2),
      neg_log_p2 = round(neg_log_p2, 2),
      effect_direction = case_when(
        Z1 > 0 & Z2 > 0 ~ "Same (+/+)",
        Z1 < 0 & Z2 < 0 ~ "Same (-/-)", 
        Z1 > 0 & Z2 < 0 ~ "Opposite (+/-)",
        Z1 < 0 & Z2 > 0 ~ "Opposite (-/+)",
        TRUE ~ "Mixed"
      )
    ) %>%
    select(chr, position_mb, snp, Z1, Z2, neg_log_p1, neg_log_p2, neg_log_lfdr, effect_direction)  # Changed rsid to snp
  
  return(top_snps)
}

# Main analysis - CRC-IBD regional plots for ALL chromosomes
main_regional_analysis_all_chr <- function() {
  
  outdir <- here("Fig1", "output")
  
  # Load pleiotropy results - CRC-IBD only
  message("Loading CRC-IBD pleiotropy results...")
  
  ibd_file <- file.path(outdir, "Fig4_data_aID_less_strict2_newlfdr.txt")
  if(file.exists(ibd_file)) {
    dt_ibd <- fread(ibd_file)
    message("Loaded CRC-IBD data: ", nrow(dt_ibd), " SNPs")
  } else {
    stop("CRC-IBD file not found: ", ibd_file)
  }
  
  # Check if snp column exists, if not create placeholder
  if(!"snp" %in% colnames(dt_ibd)) {
    message("Warning: 'snp' column not found. Creating placeholder identifiers.")
    dt_ibd$snp <- paste0("Chr", dt_ibd$chr, ":", dt_ibd$bp)
  }
  
  # Get all available chromosomes in the data
  available_chrs <- sort(unique(dt_ibd$chr))
  available_chrs <- available_chrs[available_chrs %in% c(1:22, "X", "Y")]  # Only standard chromosomes
  
  message("\n=== CREATING REGIONAL PLOTS FOR ALL CHROMOSOMES ===")
  message("Available chromosomes in data: ", paste(available_chrs, collapse = ", "))
  
  # Identify significant regions for reference
  ibd_regions <- identify_top_regions(dt_ibd, "CRC-IBD", lfdr_threshold = 0.001)
  
  if(!is.null(ibd_regions)) {
    message("Found significant pleiotropy regions on chromosomes: ", 
            paste(unique(ibd_regions$chr), collapse = ", "))
  }
  
  regional_results <- list()
  
  # Create regional plots for ALL chromosomes
  for(chr in available_chrs) {  # use c(1) for testing, available_chrs for all final plots
    message("\nProcessing chromosome ", chr, "...")
    
    # Get chromosome data
    chr_data <- dt_ibd %>% filter(chr == !!chr)
    
    if(nrow(chr_data) == 0) {
      message("No data for chromosome ", chr)
      next
    }
    
    # Check if this chromosome has significant pleiotropy signals
    chr_regions <- NULL
    if(!is.null(ibd_regions)) {
      chr_regions <- ibd_regions %>% filter(chr == !!chr)
    }
    
    if(!is.null(chr_regions) && nrow(chr_regions) > 0) {
      # Use the most significant region for this chromosome
      region <- chr_regions %>% slice_min(min_lfdr, n = 1)
      start_bp <- region$flank_start[1]
      end_bp <- region$flank_end[1]
      region_id <- region$region_id[1]
      message("Using significant region: ", region_id)
    } else {
      # For chromosomes without significant pleiotropy, create a region around the most significant individual SNP
      most_sig_snp <- chr_data %>% slice_min(lfdr_new, n = 1)
      center_bp <- most_sig_snp$bp[1]
      start_bp <- max(1, center_bp - 1e6)  # 1Mb flanking
      end_bp <- center_bp + 1e6
      region_id <- paste0("Chr", chr, "_", round(center_bp/1e6, 1), "Mb")
      message("No significant pleiotropy found. Using region around most significant SNP: ", region_id)
    }
    
    # Create regional plot
    result <- create_regional_plot(
      dt = dt_ibd,
      chr = chr,
      start_bp = start_bp,
      end_bp = end_bp,
      trait1_name = "CRC",
      trait2_name = "IBD",
      region_id = region_id
    )
    
    if(!is.null(result)) {
      # Save plot with chromosome number in filename
      filename <- paste0("Regional_CRC_IBD_Chr", sprintf("%02d", as.numeric(chr)), "_", region_id, ".png")
      ggsave(
        file.path(outdir, filename),
        result$plot,
        width = 12, height = 10, dpi = 300
      )
      
      # Analyze top SNPs
      top_snps <- analyze_region_snps(result$data, top_n = 3)
      message("Top 3 SNPs in ", region_id, ":")
      print(top_snps)
      
      regional_results[[paste0("Chr", chr, "_", region_id)]] <- result
      
      message("Saved: ", filename)
    }
  }
  
  # Final summary
  message("\n=== FINAL SUMMARY ===")
  message("Created regional plots for ", length(regional_results), " chromosomes")
  
  # Create summary table
  summary_list <- lapply(regional_results, function(x) {
    data.frame(
      chromosome = x$summary$chr,
      region = x$summary$region,
      n_snps = x$summary$n_snps,
      n_significant = x$summary$n_significant,
      max_pleiotropy = round(x$summary$max_neg_log_lfdr, 2),
      trait_correlation = round(x$summary$trait1_cor, 3),
      top_pleiotropy_snp = x$summary$top_pleiotropy_snp %||% "None",
      stringsAsFactors = FALSE
    )
  })
  
  summary_table <- do.call(rbind, summary_list)
  summary_table <- summary_table[order(summary_table$chromosome), ]
  
  print(summary_table)
  
  return(regional_results)
}
# Run the analysis
results <- main_regional_analysis_all_chr()

# Load data for sanity check
library(data.table); library(dplyr)
dt_ibd <- fread("Fig1/output/Fig4_data_aID_less_strict2_newlfdr.txt")
ibd_regions <- identify_top_regions(dt_ibd, "CRC-IBD", 0.001)

# Sanity check plot - CRC-IBD only
message("\n=== CREATING SANITY CHECK MANHATTAN PLOT ===")
sig_ibd <- dt_ibd[lfdr_new < 0.001]

library(ggplot2)
manhattan_check <- ggplot(sig_ibd, aes(x = as.numeric(chr) + bp/1e8, y = -log10(lfdr_new))) +
  geom_point(alpha = 0.6, size = 0.4, color = "darkblue") +
  scale_x_continuous(breaks = 1:22, labels = 1:22) +
  labs(title = "CRC-IBD Pleiotropy Signals (lfdr < 0.001)",
       x = "Chromosome", y = "-log10(lfdr)") +
  theme_minimal()

ggsave("Fig1/output/Manhattan_CRC_IBD_sanity_check.png", manhattan_check, 
       width = 14, height = 6, dpi = 300)

message("Sanity check Manhattan plot saved as: Manhattan_CRC_IBD_sanity_check.png")




