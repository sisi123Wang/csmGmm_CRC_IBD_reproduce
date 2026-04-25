#!/usr/bin/env Rscript
# Fig4/plot_counts_and_manhattans.R
# Usage: Rscript plot_counts_and_manhattans.R <aID> <Snum>

.libPaths(c("/home/swang25/R/ubuntu/4.4.1", Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(here)
})

args <- commandArgs(trailingOnly=TRUE)
aID  <- as.integer(args[1])  # 1=crc_ob pleio, 2=crc_ibd pleio, 3=crc_ob corr, 4=crc_ibd corr
Snum <- as.integer(args[2])  # 1=strict (0.01), 2=relaxed (0.10)

# threshold
thr <- if (Snum == 1) 0.01 else 0.10

# which pair
pair <- if (aID %in% c(1,3)) "crc_ob" else "crc_ibd"

# output dir
outdir <- here("Fig4","output")

# helper: read one lfdr vector
read_lfdr <- function(method, suffix) {
  file <- file.path(outdir, paste0("Fig4_data_aID", aID, "_", suffix, ".txt"))
  if (!file.exists(file)) stop("Missing file: ", file)
  x <- read_tsv(file, col_names=TRUE, show_col_types=FALSE)[[1]]
  tibble(Method=method, lfdr=x)
}

if (aID %in% c(1,2)) {
  # ───────────────────────────────────────────────────────────────
  # Pleio: read all six methods
  methods <- tribble(
    ~method,   ~suffix,
    "New",     "newlfdr",
    "HDMT",    "hdmt",
    "DACT",    "DACTp",
    "Kernel",  "kernel",
    "Spline7", "df7",
    "Spline50","df50"
  )
  
  lfdr_tbl <- methods %>%
    mutate(data = map2(method, suffix, ~ read_lfdr(.x, .y))) %>%
    pull(data) %>%
    bind_cols() %>%      # binds Method1,lfdr1,Method2,lfdr2,...
    pivot_longer(everything(),
                 names_to=c(".value","set"),
                 names_pattern="(Method|lfdr)(.)") %>%
    select(-set)
  
  # counts
  counts <- lfdr_tbl %>%
    group_by(Method) %>%
    summarise(nReject = sum(lfdr <= thr), .groups="drop") %>%
    mutate(Pair=pair, Snum=Snum) %>%
    select(Pair, Method, nReject, Snum)
  
  write_tsv(counts, file.path(outdir, paste0(pair,"_counts_S",Snum,".txt")))
  message("Saved counts to ", file.path(outdir, paste0(pair,"_counts_S",Snum,".txt")))
  
  # prepare per‐SNP table for Manhattan
  # re‐read the SNP metadata + lfdr_new
  snp_dat <- read_tsv(
    file.path(outdir, paste0(pair,"_newlfdr.tsv")),
    show_col_types=FALSE
  ) %>%
    rename(lfdr_new = lfdr) %>%
    select(chr, bp, snp, lfdr_new)
  
  plot_tbl <- snp_dat %>%
    mutate(hit = lfdr_new <= thr,
           minuslog10 = -log10(lfdr_new))
  
  p <- ggplot(plot_tbl,
              aes(x=bp/1e6, y=minuslog10, colour=hit)) +
    geom_point(size=0.5, alpha=0.8) +
    facet_wrap(~ factor(chr, levels=1:22), nrow=2, scales="free_x") +
    geom_hline(yintercept=-log10(thr), linetype="dashed") +
    scale_colour_manual(values=c(`TRUE`="#d62728", `FALSE`="grey60"),
                        labels=c(`TRUE`="lfdr≤thr", `FALSE`="NS"),
                        name=NULL) +
    labs(x="Position (Mb)",
         y=expression(-log[10](lfdr[new])),
         title=paste("Two‑way pleiotropy –", toupper(pair))) +
    theme_bw(base_size=11) +
    theme(legend.position="top",
          panel.spacing.x=unit(0.2,"lines"))
  
  ggsave(file.path(outdir, paste0("Figure4_",pair,"_",Snum,".pdf")),
         p, width=10, height=5.5)
  message("Saved Manhattan plot to ", file.path(outdir, paste0("Figure4_",pair,"_",Snum,".pdf")))
  
} else {
  # ───────────────────────────────────────────────────────────────
  # Correlation: only one lfdr vector
  cor_dat <- read_tsv(
    file.path(outdir, paste0("Fig4_data_aID",aID,"_cor_lfdr.txt")),
    show_col_types=FALSE
  ) %>% rename(lfdr_cor = lfdr)
  
  # counts
  nrej <- sum(cor_dat$lfdr_cor <= thr)
  count_tbl <- tibble(
    Pair    = pair,
    Method  = "Correlation",
    nReject = nrej,
    Snum    = Snum
  )
  write_tsv(count_tbl,
            file.path(outdir, paste0(pair,"_corr_counts_S",Snum,".txt")))
  message("Saved correlation counts to ",
          file.path(outdir, paste0(pair,"_corr_counts_S",Snum,".txt")))
  
  # Manhattan—plot against Z1 from the original SNP file:
  meta <- fread(here("Data", ifelse(pair=="crc_ob","crc_obesity.txt","crc_ibd.txt")),
                select=c("chromosome","base_pair_location","rsid.x"))
  names(meta) <- c("chr","bp","snp")
  
  plot_dt <- meta %>%
    inner_join(cor_dat, by=c("snp")) %>%
    mutate(hit        = lfdr_cor <= thr,
           minuslog10 = -log10(lfdr_cor))
  
  p <- ggplot(plot_dt,
              aes(x=bp/1e6, y=minuslog10, colour=hit)) +
    geom_point(size=0.5, alpha=0.8) +
    facet_wrap(~ factor(chr, levels=1:22), nrow=2, scales="free_x") +
    geom_hline(yintercept=-log10(thr), linetype="dashed") +
    scale_colour_manual(values=c(`TRUE`="#d62728", `FALSE`="grey60"),
                        labels=c(`TRUE`="lfdr≤thr", `FALSE`="NS"),
                        name=NULL) +
    labs(x="Position (Mb)",
         y=expression(-log[10](lfdr[cor])),
         title=paste("Simple correlation –", toupper(pair))) +
    theme_bw(base_size=11) +
    theme(legend.position="top",
          panel.spacing.x=unit(0.2,"lines"))
  
  ggsave(file.path(outdir, paste0("Figure4_corr_",pair,"_",Snum,".pdf")),
         p, width=10, height=5.5)
  message("Saved correlation Manhattan to ",
          file.path(outdir, paste0("Figure4_corr_",pair,"_",Snum,".pdf")))
}
