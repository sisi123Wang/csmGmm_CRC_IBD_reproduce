#!/usr/bin/env Rscript
# Fig1/make_fig4.R
# Two-way pleiotropy analysis for CRC trait pairs
#
# Usage:
#   Rscript make_fig4.R <aID 1-2> <Snum>
#
# aID = 1 : CRC–Obesity
# aID = 2 : CRC–IBD
#
# Note:
# Although cfg still records whether correlation was considered before,
# the current script runs the independence method for both pairs.

options(vroom.integer64 = "double")

.libPaths(c(
  "/home/swang25/R/ubuntu/4.4.1",
  Sys.getenv("R_LIBS_USER"),
  .libPaths()
))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(csmGmm)   # symm_fit_ind_EM / symm_fit_cor_EM
  library(locfdr)
  library(here)
})

here::i_am("Fig1/make_fig4.R") #change path.

# ---------------------------
# Command-line arguments
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript make_fig4.R <aID 1-2> <Snum>")
}

aID  <- as.integer(args[1])
Snum <- as.integer(args[2])

# ---------------------------
# Analysis configuration
# ---------------------------
cfg <- list(
  `1` = list(
    pair = "crc_ob",
    infile = here("Data", "crc_obesity2.txt"),
    pleio = TRUE,
    use_correlation = TRUE,
    beta1 = "es_trait1",
    p1 = "p_crc",
    beta2 = "es_trait2",
    p2 = "p_ob"
  ),
  `2` = list(
    pair = "crc_ibd",
    infile = here("Data", "crc_ibd.txt"),
    pleio = TRUE,
    use_correlation = FALSE,
    beta1 = "es_trait1",
    p1 = "p_crc",
    beta2 = "es_trait2",
    p2 = "p_ibd"
  )
)

item <- cfg[[as.character(aID)]]

if (is.null(item)) {
  stop("aID must be 1 or 2 (pleiotropy analysis only)")
}
if (!file.exists(item$infile)) {
  stop("Input file not found: ", item$infile)
}

message(
  "Running aID = ", aID,
  "; pleio? ", item$pleio,
  "; use_correlation? ", item$use_correlation,
  "; infile = ", item$infile
)

# ---------------------------
# Read input data
# ---------------------------
dt <- fread(item$infile, na.strings = c("NA", "#NA", ""))

message("   Columns in dt: ", paste(names(dt), collapse = ", "))
message("   Total rows: ", nrow(dt))

# ---------------------------
# Required columns
# ---------------------------
beta1 <- item$beta1
p1col <- item$p1
beta2 <- item$beta2
p2col <- item$p2

required <- c("chromosome", "base_pair_location", beta1, p1col, beta2, p2col)

if (!all(required %in% names(dt))) {
  stop("Raw summary file must contain: ", paste(required, collapse = ", "))
}

# ---------------------------
# Convert to numeric and compute Z-scores
# SE is reconstructed as |beta| / |qnorm(p/2)|
# Then Z = beta / SE
# ---------------------------
dt2 <- dt %>%
  mutate(
    beta1n = as.numeric(.data[[beta1]]),
    p1n    = as.numeric(.data[[p1col]]),
    beta2n = as.numeric(.data[[beta2]]),
    p2n    = as.numeric(.data[[p2col]])
  )

message(
  sprintf(
    "After coercion: NA %s = %d, NA %s = %d; NA %s = %d, NA %s = %d",
    beta1, sum(is.na(dt2$beta1n)),
    p1col, sum(is.na(dt2$p1n)),
    beta2, sum(is.na(dt2$beta2n)),
    p2col, sum(is.na(dt2$p2n))
  )
)

dt3 <- dt2 %>%
  mutate(
    se1 = abs(beta1n) / abs(qnorm(p1n / 2)),
    se2 = abs(beta2n) / abs(qnorm(p2n / 2)),
    Z1  = beta1n / se1,
    Z2  = beta2n / se2
  ) %>%
  filter(!is.na(Z1), !is.na(Z2), is.finite(Z1), is.finite(Z2))

message("   Rows after computing Z and filtering: ", nrow(dt3))

if (nrow(dt3) == 0) {
  stop("No valid rows after computing Z for pleiotropy.")
}

# ---------------------------
# Build analysis matrix
# ---------------------------
zmat <- as.matrix(dt3 %>% select(Z1, Z2))

# kept from your original script, although not used later
pA <- dt3[[p1col]]
pB <- dt3[[p2col]]

# Keep genomic information for downstream Manhattan plots
snp_col <- if ("rsid.x" %in% names(dt3)) {
  "rsid.x"
} else if ("rsid" %in% names(dt3)) {
  "rsid"
} else {
  NULL
}

if (!is.null(snp_col)) {
  dat0 <- dt3 %>%
    transmute(
      chr = chromosome,
      bp  = base_pair_location,
      snp = .data[[snp_col]],
      Z1, Z2
    )
} else {
  dat0 <- dt3 %>%
    transmute(
      chr = chromosome,
      bp  = base_pair_location,
      Z1, Z2
    )
}

# ---------------------------
# Clip extreme Z values
# This stabilizes very large values before EM fitting
# ---------------------------
zmat[zmat >  8.1] <-  8.1
zmat[zmat < -8.1] <- -8.1

message(
  "   zmat dims: ", paste(dim(zmat), collapse = " x "),
  "; anyNA? ", any(is.na(zmat)),
  "; anyInf? ", any(is.infinite(zmat))
)

# ---------------------------
# Output setup
# ---------------------------
outdir <- here("Fig1", "output")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

root <- file.path(outdir, paste0("Fig4_data_aID_ob2", aID))
# Alternative roots used before for less strict convergence option:
# root <- file.path(outdir, paste0("Fig4_data_aID_less_strict", aID))

# ---------------------------
# EM tuning parameters
# ---------------------------
newEps <- 1e-5
# newEps <- 1e-3   # less strict convergence option

initPi <- list(
  c(0.82),
  c(0.02, 0.02),
  c(0.02, 0.02),
  c(0.1)
)

initMu <- list(
  matrix(0, 2, 1),
  matrix(c(0, 3, 0, 6), 2),
  matrix(c(3, 0, 6, 0), 2),
  matrix(c(8, 8), 2)
)

# Alternative liberal initialization kept here for reference:
# initPi <- list(
#   c(0.60),
#   c(0.10, 0.10),
#   c(0.10, 0.10),
#   c(0.10)
# )
#
# initMu <- list(
#   matrix(0, 2, 1),
#   matrix(c(0, 1.5, 0, 3), 2),
#   matrix(c(1.5, 0, 3, 0), 2),
#   matrix(c(3, 3), 2)
# )

# ---------------------------
# Fit csmGmm model
# Current script uses independence method for both trait pairs
# ---------------------------
message("Using symm_fit_ind_EM on ", nrow(zmat), " rows")

newRes <- symm_fit_ind_EM(
  testStats  = zmat,
  initMuList = initMu,
  initPiList = initPi,
  eps        = newEps
)

# ---------------------------
# Save output for Manhattan plot
# ---------------------------
results_df <- cbind(dat0, lfdr_new = newRes$lfdrResults)

write.table(
  results_df,
  paste0(root, "_newlfdr.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Finished aID = ", aID, " | pleiotropy analysis | using independence method")
