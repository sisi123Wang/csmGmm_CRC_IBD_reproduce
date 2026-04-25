#!/usr/bin/env Rscript
# Fig1/make_fig4.R – two‑way pleiotropy & simple correlation for CRC pairs
# Usage: Rscript make_fig4.R <aID 1-2> <Snum>
# aID=1: CRC–Ob pleio; 2: CRC–IBD pleio

options(vroom.integer64 = "double")
.libPaths(c("/home/swang25/R/ubuntu/4.4.1",
            Sys.getenv("R_LIBS_USER"),
            .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(csmGmm)    # for symm_fit_ind_EM / symm_fit_cor_EM
  library(locfdr)    # for null_estimation
  library(here)
})

here::i_am("Fig1/make_fig4.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2) stop("Usage: Rscript make_fig4.R <aID 1-2> <Snum>")
aID  <- as.integer(args[1])
Snum <- as.integer(args[2])  # not used here but reserved

# PLEIOTROPY ANALYSIS
cfg <- list(
  `1` = list(pair="crc_ob",
             infile= here("Data","crc_obesity2.txt"),
             pleio = TRUE,
             use_correlation = TRUE,  # CRC-Obesity uses correlation method--MAYBE CHANGE TO FALSE--tried but not good.
             beta1="es_trait1", p1="p_crc",
             beta2="es_trait2", p2="p_ob"),
  `2` = list(pair="crc_ibd",
             infile= here("Data","crc_ibd.txt"),
             pleio = TRUE,
             use_correlation = FALSE,  # CRC-IBD uses independence method
             beta1="es_trait1", p1="p_crc",
             beta2="es_trait2", p2="p_ibd")
)
item <- cfg[[as.character(aID)]]
if(is.null(item)) stop("aID must be 1 or 2 (pleiotropy analysis only)")
if(!file.exists(item$infile)) stop("Input file not found: ", item$infile)

message("Running aID = ", aID, "; pleio? ", item$pleio, "; use_correlation? ", item$use_correlation, "; infile = ", item$infile)
dt <- fread(item$infile, na.strings=c("NA","#NA",""))
message("   Columns in dt: ", paste(names(dt), collapse=", "))
message("   Total rows: ", nrow(dt))

# Z-matrix and pA,pB - PLEIOTROPY ONLY
# PLEIOTROPY: raw summary must have chromosome, base_pair_location,
#   es_trait1, p_crc, es_trait2, p_ob or p_ibd.
beta1 <- item$beta1; p1col <- item$p1
beta2 <- item$beta2; p2col <- item$p2
required <- c("chromosome","base_pair_location", beta1, p1col, beta2, p2col)
if(!all(required %in% names(dt))) {
  stop("Raw summary file must contain: ", paste(required, collapse=", "))
}
# coerce numeric
dt2 <- dt %>%
  mutate(
    beta1n = as.numeric(.data[[beta1]]),
    p1n    = as.numeric(.data[[p1col]]),
    beta2n = as.numeric(.data[[beta2]]),
    p2n    = as.numeric(.data[[p2col]])
  )
message(sprintf("After coercion: NA %s = %d, NA %s = %d; NA %s = %d, NA %s = %d",
                beta1, sum(is.na(dt2$beta1n)),
                p1col, sum(is.na(dt2$p1n)),
                beta2, sum(is.na(dt2$beta2n)),
                p2col, sum(is.na(dt2$p2n))))
# fill missing SE via |β|/|qnorm(p/2)|, then Z = β/SE
dt3 <- dt2 %>%
  mutate(
    se1 = abs(beta1n) / abs(qnorm(p1n/2)),
    se2 = abs(beta2n) / abs(qnorm(p2n/2)),
    Z1  = beta1n / se1,
    Z2  = beta2n / se2
  ) %>%
  # filter valid
  filter(!is.na(Z1), !is.na(Z2), is.finite(Z1), is.finite(Z2))
message("   Rows after computing Z and filtering: ", nrow(dt3))
if(nrow(dt3)==0) stop("No valid rows after computing Z for pleiotropy.")
zmat <- as.matrix(dt3 %>% select(Z1,Z2))
pA <- dt3[[p1col]]
pB <- dt3[[p2col]]
# optionally, keep chr/bp/snp
snp_col <- if("rsid.x" %in% names(dt3)) "rsid.x" else if("rsid" %in% names(dt3)) "rsid" else NULL
if(!is.null(snp_col)) {
  dat0 <- dt3 %>% transmute(chr=chromosome, bp=base_pair_location, snp=.data[[snp_col]], Z1, Z2)
} else {
  dat0 <- dt3 %>% transmute(chr=chromosome, bp=base_pair_location, Z1, Z2)
}

# Clip extremes
zmat[zmat >  8.1] <-  8.1
zmat[zmat < -8.1] <- -8.1
message("   zmat dims: ", paste(dim(zmat), collapse=" x "),
        "; anyNA? ", any(is.na(zmat)), "; anyInf? ", any(is.infinite(zmat)))

# set up output
outdir <- here("Fig1","output")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
root   <- file.path(outdir, paste0("Fig4_data_aID_ob2", aID))
#root <- file.path(outdir, paste0("Fig4_data_aID_less_strict", aID))  # Less strict
#root <- file.path(outdir, paste0("Fig4_data_aID_reharm", aID))

newEps <- 1e-5
#newEps <- 1e-3  # 100x less strict convergence

# EM init
initPi <- list(c(0.82),
               c(0.02,0.02),
               c(0.02,0.02),
               c(0.1))
initMu <- list(matrix(0,2,1),
               matrix(c(0,3,0,6),2),
               matrix(c(3,0,6,0),2),
               matrix(c(8,8),2))

#EM init - LIBERAL PARAMETERS
#initPi <- list(c(0.60),        # Less null mass 82% to 60%
#               c(0.10,0.10),   # More signal mass (was 0.02,0.02)
#               c(0.10,0.10),   # More signal mass (was 0.02,0.02)
#               c(0.10))        # More pleiotropy mass (was 0.1)

#initMu <- list(matrix(0,2,1),
#               matrix(c(0,1.5,0,3),2),    # Lower thresholds (was 0,3,0,6)
#               matrix(c(1.5,0,3,0),2),    # Lower thresholds (was 3,0,6,0)
 #              matrix(c(3,3),2))          # Lower thresholds (was 8,8)

# methods
#if(item$use_correlation) {
  # Use correlation method for CRC-Obesity pairs (they share material)
#  message("Using symm_fit_cor_EM on ", nrow(zmat), " rows")
#  rho  <- cor(zmat[,1], zmat[,2])
#  cmat <- matrix(c(1,rho,rho,1), nrow=2, ncol=2)
#  corRes <- symm_fit_cor_EM(testStats   = zmat,
#                            corMat      = cmat,
#                            initMuList  = initMu,
#                            initPiList  = initPi,
#                            eps         = newEps)
  # Combine genomic info with lfdr results for Manhattan plots
#  results_df <- cbind(dat0, lfdr_new = corRes$lfdrResults)
#  write.table(results_df,
#              paste0(root,"_newlfdr.txt"),
#              sep="\t", quote=FALSE, row.names=FALSE)
#} else {
  # Use independence method for CRC-IBD pairs
#  message("Using symm_fit_ind_EM on ", nrow(zmat), " rows")
#  newRes <- symm_fit_ind_EM(testStats = zmat,
#                            initMuList  = initMu,
#                            initPiList  = initPi,
#                            eps         = newEps)
  # Combine genomic info with lfdr results for Manhattan plots
#  results_df <- cbind(dat0, lfdr_new = newRes$lfdrResults)
#  write.table(results_df,
#              paste0(root,"_newlfdr.txt"),
#              sep="\t", quote=FALSE, row.names=FALSE)
#}

#message("Finished aID = ", aID, " | pleiotropy analysis | use_correlation = ", item$use_correlation)

# methods--CRC-Obesity pairs now use independence method since they no longer share material
message("Using symm_fit_ind_EM on ", nrow(zmat), " rows")
newRes <- symm_fit_ind_EM(testStats = zmat,
                          initMuList  = initMu,
                          initPiList  = initPi,
                          eps         = newEps)

# Combine genomic info with lfdr results for Manhattan plots
results_df <- cbind(dat0, lfdr_new = newRes$lfdrResults)
write.table(results_df,
            paste0(root,"_newlfdr.txt"),
            sep="\t", quote=FALSE, row.names=FALSE)

message("Finished aID = ", aID, " | pleiotropy analysis | using independence method")



