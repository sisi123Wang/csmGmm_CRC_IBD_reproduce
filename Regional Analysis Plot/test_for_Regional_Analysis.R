#!/usr/bin/env Rscript
# Regional pleiotropy analysis for CRC vs IBD
# Purpose:
# 1. Identify top pleiotropic regions using lfdr
# 2. Plot one regional figure per chromosome
# 3. Summarize the top SNPs in each plotted region

.libPaths("/home/swang25/R/ubuntu/4.4.1")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(here)
  library(gridExtra)
})

# ---------------------------
# Helper: find top regions
# ---------------------------
# Significant SNPs are grouped into loci if they are within min_distance bp.
identify_top_regions <- function(dt, lfdr_threshold = 0.001, min_distance = 1e6) {
  sig_snps <- dt %>%
    filter(lfdr_new < lfdr_threshold) %>%
    arrange(chr, bp, lfdr_new)

  if (nrow(sig_snps) == 0) {
    message("No significant SNPs found at lfdr < ", lfdr_threshold)
    return(NULL)
  }

  regions <- sig_snps %>%
    group_by(chr) %>%
    mutate(region_group = cumsum(c(1, diff(bp) > min_distance))) %>%
    group_by(chr, region_group) %>%
    summarise(
      start_bp   = min(bp),
      end_bp     = max(bp),
      center_bp  = round(mean(bp)),
      n_snps     = n(),
      min_lfdr   = min(lfdr_new),
      lead_snp   = snp[which.min(lfdr_new)][1],
      lead_bp    = bp[which.min(lfdr_new)][1],
      .groups    = "drop"
    ) %>%
    arrange(min_lfdr) %>%
    mutate(
      region_id   = paste0("Chr", chr, "_", round(center_bp / 1e6, 1), "Mb"),
      flank_start = pmax(1, start_bp - 5e5),
      flank_end   = end_bp + 5e5
    )

  return(regions)
}

# ---------------------------
# Helper: prepare one region
# ---------------------------
prepare_region_data <- function(dt, chr, start_bp, end_bp, lfdr_threshold = 0.01) {
  region_data <- dt %>%
    filter(chr == !!chr, bp >= start_bp, bp <= end_bp) %>%
    mutate(
      position_mb   = bp / 1e6,
      significant   = lfdr_new < lfdr_threshold,
      neg_log_lfdr  = -log10(lfdr_new),
      neg_log_p_crc = -log10(2 * pnorm(abs(Z1), lower.tail = FALSE)),
      neg_log_p_ibd = -log10(2 * pnorm(abs(Z2), lower.tail = FALSE))
    )

  if (nrow(region_data) == 0) return(NULL)
  region_data
}

# ---------------------------
# Helper: make 3-panel plot
# ---------------------------
create_regional_plot <- function(region_data, chr, region_id,
                                 trait1_name = "CRC",
                                 trait2_name = "IBD",
                                 lfdr_threshold = 0.01) {
  x_min <- min(region_data$position_mb)
  x_max <- max(region_data$position_mb)

  top_pleio  <- region_data %>% slice_min(lfdr_new, n = 1)
  top_crc    <- region_data %>% slice_max(neg_log_p_crc, n = 1)
  top_ibd    <- region_data %>% slice_max(neg_log_p_ibd, n = 1)

  p1 <- ggplot(region_data, aes(position_mb, neg_log_lfdr)) +
    geom_point(aes(color = significant), alpha = 0.7, size = 1) +
    geom_hline(yintercept = -log10(lfdr_threshold), linetype = "dashed", color = "red") +
    geom_text(
      data = top_pleio,
      aes(label = snp),
      vjust = -0.5, size = 3, color = "black"
    ) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
    labs(
      title = paste0("Pleiotropy Signal: ", trait1_name, " vs ", trait2_name),
      subtitle = region_id,
      x = NULL,
      y = "-log10(lfdr)"
    ) +
    coord_cartesian(xlim = c(x_min, x_max)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), legend.position = "right")

  p2 <- ggplot(region_data, aes(position_mb, neg_log_p_crc)) +
    geom_point(aes(color = Z1 > 0), alpha = 0.7, size = 1) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "grey50") +
    geom_text(
      data = top_crc,
      aes(label = snp),
      vjust = -0.5, size = 3, color = "black"
    ) +
    scale_color_manual(
      values = c("FALSE" = "blue", "TRUE" = "red"),
      labels = c("Negative", "Positive")
    ) +
    labs(x = NULL, y = paste0("-log10(P) ", trait1_name)) +
    coord_cartesian(xlim = c(x_min, x_max)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), legend.position = "right")

  p3 <- ggplot(region_data, aes(position_mb, neg_log_p_ibd)) +
    geom_point(aes(color = Z2 > 0), alpha = 0.7, size = 1) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "grey50") +
    geom_text(
      data = top_ibd,
      aes(label = snp),
      vjust = -0.5, size = 3, color = "black"
    ) +
    scale_color_manual(
      values = c("FALSE" = "blue", "TRUE" = "red"),
      labels = c("Negative", "Positive")
    ) +
    labs(
      x = paste0("Position (Mb) - Chromosome ", chr),
      y = paste0("-log10(P) ", trait2_name)
    ) +
    coord_cartesian(xlim = c(x_min, x_max)) +
    theme_minimal() +
    theme(legend.position = "right")

  combined_plot <- grid.arrange(p1, p2, p3, ncol = 1, heights = c(1, 0.8, 0.8))

  summary_info <- list(
    chr = chr,
    region = region_id,
    n_snps = nrow(region_data),
    n_significant = sum(region_data$significant),
    max_neg_log_lfdr = max(region_data$neg_log_lfdr),
    cor_z = cor(region_data$Z1, region_data$Z2, use = "complete.obs"),
    top_pleio_snp = top_pleio$snp[1]
  )

  list(plot = combined_plot, data = region_data, summary = summary_info)
}

# ---------------------------
# Helper: summarize top SNPs
# ---------------------------
summarize_top_snps <- function(region_data, top_n = 3) {
  region_data %>%
    arrange(lfdr_new) %>%
    slice_head(n = top_n) %>%
    transmute(
      chr,
      position_mb = round(bp / 1e6, 2),
      snp,
      Z1 = round(Z1, 3),
      Z2 = round(Z2, 3),
      neg_log_p_crc = round(neg_log_p_crc, 2),
      neg_log_p_ibd = round(neg_log_p_ibd, 2),
      neg_log_lfdr = round(-log10(lfdr_new), 2),
      effect_direction = case_when(
        Z1 > 0 & Z2 > 0 ~ "Same (+/+)",
        Z1 < 0 & Z2 < 0 ~ "Same (-/-)",
        Z1 > 0 & Z2 < 0 ~ "Opposite (+/-)",
        Z1 < 0 & Z2 > 0 ~ "Opposite (-/+)",
        TRUE ~ "Mixed"
      )
    )
}

# ---------------------------
# Main analysis
# ---------------------------
main <- function() {
  outdir <- here("Fig1", "output")
  infile <- file.path(outdir, "Fig4_data_aID_less_strict2_newlfdr.txt")

  if (!file.exists(infile)) {
    stop("Input file not found: ", infile)
  }

  dt <- fread(infile)

  if (!"snp" %in% names(dt)) {
    dt[, snp := paste0("Chr", chr, ":", bp)]
  }

  regions <- identify_top_regions(dt, lfdr_threshold = 0.001)

  available_chrs <- sort(unique(dt$chr))
  available_chrs <- available_chrs[available_chrs %in% c(1:22, "X", "Y")]

  results <- list()

  for (chr_i in available_chrs) {
    message("Processing chromosome ", chr_i, "...")

    chr_data <- dt %>% filter(chr == !!chr_i)
    if (nrow(chr_data) == 0) next

    chr_regions <- NULL
    if (!is.null(regions)) {
      chr_regions <- regions %>% filter(chr == !!chr_i)
    }

    if (!is.null(chr_regions) && nrow(chr_regions) > 0) {
      region_row <- chr_regions %>% slice_min(min_lfdr, n = 1)
      start_bp <- region_row$flank_start[1]
      end_bp   <- region_row$flank_end[1]
      region_id <- region_row$region_id[1]
    } else {
      best_snp <- chr_data %>% slice_min(lfdr_new, n = 1)
      center_bp <- best_snp$bp[1]
      start_bp <- max(1, center_bp - 1e6)
      end_bp   <- center_bp + 1e6
      region_id <- paste0("Chr", chr_i, "_", round(center_bp / 1e6, 1), "Mb")
    }

    region_data <- prepare_region_data(dt, chr_i, start_bp, end_bp)
    if (is.null(region_data)) next

    result <- create_regional_plot(region_data, chr = chr_i, region_id = region_id)

    file_name <- paste0("Regional_CRC_IBD_", region_id, ".png")

    ggsave(
      filename = file.path(outdir, file_name),
      plot = result$plot,
      width = 12,
      height = 10,
      dpi = 300
    )

    message("Saved: ", file_name)
    print(summarize_top_snps(result$data, top_n = 3))

    results[[paste0("Chr", chr_i)]] <- result
  }

  summary_table <- bind_rows(lapply(results, function(x) {
    data.frame(
      chromosome = x$summary$chr,
      region = x$summary$region,
      n_snps = x$summary$n_snps,
      n_significant = x$summary$n_significant,
      max_pleiotropy = round(x$summary$max_neg_log_lfdr, 2),
      trait_correlation = round(x$summary$cor_z, 3),
      top_pleiotropy_snp = x$summary$top_pleio_snp,
      stringsAsFactors = FALSE
    )
  }))

  summary_table <- summary_table %>% arrange(chromosome)
  print(summary_table)

  invisible(results)
}

results <- main()

