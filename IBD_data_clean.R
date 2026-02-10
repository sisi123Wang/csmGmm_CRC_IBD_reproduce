---
title: "IBD_data_clean"
output: html_document
date: "2025-05-26"
---

```{r setup, message=FALSE, warning=FALSE, results='asis'}
library(data.table)
library(dplyr)
library(tidyr)
library(knitr)
library(scales)
library(readr)
#datasets and metadata
data_dir <- "/home/swang25/other_summary_statistics/"
files <- tibble::tribble(
  ~id,         ~path,                                                              ~build,
  "GCST90044155", "/home/swang25/other_summary_statistics/IBD/GCST90044155_buildGRCh37.tsv.gz", "GRCh37",
  "GCST90081487", "/home/swang25/other_summary_statistics/IBD/GCST90081487_buildGRCh38.tsv.gz", "GRCh38",
  "GCST90292538", "/home/swang25/other_summary_statistics/IBD/GCST90292538.tsv", "GRCh38",
  "GCST90476065", "/home/swang25/other_summary_statistics/IBD/GCST90476065.tsv.gz", "GRCh38",
  "GCST90476067", "/home/swang25/other_summary_statistics/IBD/GCST90476067.tsv.gz", "GRCh38",
  "ibd_build37_59957_20161107.txt.gz", "/home/swang25/other_summary_statistics/IBD/ibd_build37_59957_20161107.txt.gz", "GRCh37",
  "IBD_trans_ethnic_association_summ_stats_b37.txt.gz", "/home/swang25/other_summary_statistics/IBD/IBD_trans_ethnic_association_summ_stats_b37.txt.gz","GRCh37")


md <- knitr::kable(files, format = "markdown", caption = "Colorectal GWAS Files Metadata")
cat(md, "\n")

```

```{r iterate-all-datasets, include=TRUE, echo=FALSE, message=TRUE, warning=TRUE, results='asis'}
library(dplyr)
library(tidyr)
library(knitr)
library(scales)
library(readr)

setwd("~/other_summary_statistics/IBD/")

# helper that never writes a temp file
read_gwas <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE)
}

datasets <- tibble::tribble(
  ~id,            ~file,                                                            ~sample_size,
  "GCST90044155", "GCST90044155_buildGRCh37.tsv.gz",                               456348,
  "GCST90081487", "GCST90081487_buildGRCh38.tsv.gz",                               364036,
  "GCST90292538", "GCST90292538.tsv",                                              398668,
  "GCST90476065", "GCST90476065.tsv.gz",                                           446753,
  "GCST90476067", "GCST90476067.tsv.gz",                                           447331
) %>%
  mutate(file = file.path(getwd(), file))

# dataset head check
for (i in seq_len(nrow(datasets))) {
  ds <- datasets[i, ]
  cat("\n===== Dataset:", ds$id, "=====\n")
  cat("File:", ds$file, "\n")
  
  # Read only the first 5 rows for a fast preview
  dt_head <- fread(ds$file, nrows = 5)
  
  # Display column names and a preview
  cat("Column Names:\n")
  print(names(dt_head))
  cat("\nHead:\n")
  print(dt_head)
}

for (i in seq_len(nrow(datasets))) {
  ds <- datasets[i, ]
  cat("Dataset:", ds$id, "\n\n")
  cat("File:", ds$file, "\n\n")
  
  raw <- read_gwas(ds$file)
  cat("Original dimensions:", nrow(raw), "rows ×", ncol(raw), "columns\n\n")
  
  cat("Columns:\n")
  print(names(raw)); cat("\n")
  
  df <- raw %>% select(chromosome, base_pair_location, effect_allele, other_allele, p_value)
  cat("Preview (first 5 rows):\n")
  print(kable(head(df, 5), format = "markdown")); cat("\n")
  
  clean <- df %>% filter(!is.na(chromosome), !is.na(base_pair_location), !is.na(p_value))
  cat("Removed", nrow(df) - nrow(clean), "rows with missing; remaining:", nrow(clean), "\n\n")
  cat("Sample size (N):", ds$sample_size, "\n\n")
  
  snp_summary <- tibble(
    `Total SNPs`  = nrow(df),
    `Unique SNPs` = n_distinct(df$base_pair_location)
  )
  cat("SNP counts:\n")
  print(kable(snp_summary, format = "markdown")); cat("\n")
  
  pval_counts <- df %>%
    summarise(
      `<5e-8` = sum(p_value < 5e-8,  na.rm = TRUE),
      `<5e-6` = sum(p_value < 5e-6,  na.rm = TRUE),
      `<5e-4` = sum(p_value < 5e-4,  na.rm = TRUE),
      `<0.05` = sum(p_value < 0.05, na.rm = TRUE)
    ) %>%
    pivot_longer(everything(), names_to = "Threshold", values_to = "Count")
  cat("P‑value thresholds:\n")
  print(kable(pval_counts, format = "markdown")); cat("\n\n---\n\n")
}


```


```{R}
setwd("/home/swang25/other_summary_statistics/IBD/")
here::i_am("IBD_data_clean.Rmd")

clean_gwas <- function(id, path, build) {
  cat("Cleaning", id, "(assembly:", build, ")\n\n")
  
  raw <- read_tsv(path, show_col_types = FALSE)
  n0  <- nrow(raw)
  cat("Original rows:", n0, "Columns:", paste(names(raw), collapse=", "), "\n\n")
  if ("effect_allele" %in% names(raw)) {
    raw$effect_allele <- toupper(raw$effect_allele)
  }
  if ("other_allele" %in% names(raw)) {
    raw$other_allele <- toupper(raw$other_allele)
  }

# Drop strand-ambiguous SNPs
  if (all(c("effect_allele", "other_allele") %in% names(raw))) {
  is_ambiguous <- (raw$effect_allele == "A" & raw$other_allele == "T") |
                  (raw$effect_allele == "T" & raw$other_allele == "A") |
                  (raw$effect_allele == "C" & raw$other_allele == "G") |
                  (raw$effect_allele == "G" & raw$other_allele == "C")
  cat("Strand-ambiguous SNPs removed:", sum(is_ambiguous), "\n\n")
  raw <- raw[!is_ambiguous, ]
  }
  # Filter missing
  df <- raw %>%
    filter(!is.na(chromosome), !is.na(base_pair_location), !is.na(p_value))
  n1 <- nrow(df)
  cat("After dropping missing chr/BP/p_value:", n1, "(removed", n0-n1, ")\n\n")
  df <- df %>%
    mutate(
      chromosome = as.character(chromosome),
      base_pair_location = as.integer(base_pair_location)
    )
  #remove multi base
  df <- df %>%
    filter(nchar(effect_allele) == 1 & nchar(other_allele) == 1)
  cat("Removing multibase alleles:", nrow(df), "remain\n\n")
  # remove same and other allele
  same_allele <- df$effect_allele == df$other_allele
  df <- df[!same_allele, ]
  cat("SNPs with same effect and other alleles removed:", sum(same_allele), "\n")

  
  # Deduplicate
  df <- df %>%
    distinct(chromosome, base_pair_location, .keep_all = TRUE)
  n2 <- nrow(df)
  cat("After deduplication:", n2, "(removed", n1-n2, "dups)\n\n")
  #odd ratio and beta
  if ("odds_ratio" %in% names(df) & !"beta" %in% names(df)) {
    cat("Note: odds_ratio found, converting to beta using log(OR)\n")
    df <- df %>%
      mutate(beta = log(odds_ratio))
  } else if ("odds_ratio" %in% names(df) & "beta" %in% names(df)) {
    cat("Note: odds_ratio and beta both found; using existing beta column\n")
    df <- df %>%
      mutate(beta = ifelse(is.na(beta), log(odds_ratio), beta))
  }
  
  # sort
  chr_levels <- c(as.character(1:22), "X","Y","MT")
  df <- df %>%
    mutate(chromosome = factor(chromosome, levels=chr_levels, ordered=TRUE)) %>%
    arrange(chromosome, base_pair_location) %>%
    mutate(chromosome = as.character(chromosome))
  
  cat("Preview cleaned data:\n")
  print(kable(head(df,5), format="markdown"))
  cat("\n---\n\n")
  
  out_path <- sub("\\.tsv(\\.gz)?(\\.1)?$", "_cleaned.tsv.gz", path)
  write_tsv(df, out_path)
  cat("Cleaned file:", out_path, "\n\n")
  return(out_path)
}

cleaned_paths <- apply(files, 1, function(r) 
  clean_gwas(r["id"], r["path"], r["build"])
)

#clean_gwas("GCST90476067", "/home/swang25/other_summary_statistics/IBD/GCST90044155_buildGRCh38_c", "GRCh37")


```


```{r liftover-with-chr-prefix, echo=FALSE}
setwd("~/other_summary_statistics/IBD/")

library(readr)
library(dplyr)

# assemble GRCh37‐cleaned file then convert it to GRCh38-cleaned file 
df37 <- read_tsv("GCST90044155_buildGRCh37_cleaned.tsv.gz", show_col_types = FALSE)


# coords.bed with 'chr' prefix
coords <- df37 %>%
  mutate(
    chrom   = paste0("chr", chromosome),
    start0  = as.integer(base_pair_location) - 1L,
    end1    = as.integer(base_pair_location),
    idx     = row_number()
  ) %>%
  select(chrom, start0, end1, idx)

write.table(
  coords,
  file      = "coords.bed",
  sep       = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote     = FALSE
)

# liftOver
system("./liftOver coords.bed hg19ToHg38.over.chain.gz mapped.bed unmapped.bed")

# mapped.bed with correct types
mapped <- read_tsv("mapped.bed",
                   col_names = c("chrom","start0","end1","idx"),
                   col_types = cols(
                     chrom  = col_character(),
                     start0 = col_integer(),
                     end1   = col_integer(),
                     idx    = col_integer()
                   ))

# Drop the 'chr' prefix for consistency with other files
mapped <- mapped %>% mutate(chromosome = sub("^chr", "", chrom))

# update
df38 <- df37 %>%
  slice(mapped$idx) %>%
  mutate(
    chromosome         = mapped$chromosome,
    base_pair_location = mapped$end1,
    assembly           = "GRCh38"
  )
out38 <- "GCST90044155_buildGRCh38_cleaned.tsv.gz"
write_tsv(df38, out38)

cat("Lifted", nrow(df38), "of", nrow(df37),
    "variants to GRCh38 in file:", out38, "\n")
```
```{r liftover-check, message=FALSE, warning=FALSE}
library(data.table)
ckeck_liftover <- "/home/swang25/other_summary_statistics/IBD/GCST90044155_buildGRCh38_cleaned.tsv.gz"

cat("Previewing:", basename(ckeck_liftover))
cat("File:", ckeck_liftover, "\n")
dt_head <- fread(ckeck_liftover, nrows = 5)

# display
cat("Column Names:\n")
print(names(dt_head))
cat("\nHead:\n")
print(dt_head)

```


```{R}
setwd("~/other_summary_statistics/IBD/")

# SNPs
disc <- read_tsv("GCST90081487_buildGRCh38_cleaned.tsv.gz", show_col_types = FALSE)
cat("GCST90081487 columns:\n"); print(names(disc))
cat("\nGCST90081487 data preview:\n"); print(head(disc, 5))
#
rep1 <- read_tsv("GCST90292538_buildGRCh38_cleaned.tsv.gz", show_col_types = FALSE)
cat("GCST90292538 columns:\n"); print(names(rep1))
cat("\nGCST90292538 data preview:\n"); print(head(rep1, 5))
#
rep2 <- read_tsv("GCST90044155_buildGRCh38_cleaned.tsv.gz", show_col_types = FALSE)
cat("GCST90044155 columns:\n"); print(names(rep2))
cat("\nGCST90044155 data preview:\n"); print(head(rep2, 5))
#
rep3 <- read_tsv("GCST90476065_cleaned.tsv.gz", show_col_types = FALSE)
cat("GCST90476065 columns:\n"); print(names(rep3))
cat("\nGCST90476065 data preview:\n"); print(head(rep3, 5))
#
rep4 <- read_tsv("GCST90476067_cleaned.tsv.gz", show_col_types = FALSE)
cat("GCST90476067 columns:\n"); print(names(rep4))
cat("\nGCST90476067 data preview:\n"); print(head(rep4, 5))

```



```{r allele-harmonization-final, message=FALSE, warning=FALSE}
setwd("~/other_summary_statistics/IBD/")

library(readr)
library(dplyr)
library(knitr)
library(glue)
# to convert negative beta--es_disc here
N <- 10

# top N by p-value
disc <- read_tsv("GCST90081487_buildGRCh38_cleaned.tsv.gz", show_col_types=FALSE) %>%
  mutate(
    chromosome         = as.character(chromosome),
    base_pair_location = as.integer(base_pair_location)
  ) %>%
  arrange(p_value) %>%
  slice_head(n = N) %>%
  transmute(
    chromosome, base_pair_location,
    ea_disc   = effect_allele,
    ra_disc   = other_allele,
    es_disc   = log(odds_ratio),
    se_disc   = standard_error,
    p_disc    = p_value
  )

# beta as es
rep1 <- read_tsv("GCST90292538_cleaned.tsv.gz", show_col_types=FALSE) %>%
  mutate(
    chromosome         = as.character(chromosome),
    base_pair_location = as.integer(base_pair_location)
  ) %>%
  transmute(
    chromosome, base_pair_location,
    ea_rep1 = effect_allele,
    ra_rep1 = other_allele,
    es_rep1 = beta,
    se_rep1 = standard_error,
    p_rep1  = p_value
  )

# beta as es
rep2 <- read_tsv("GCST90044155_buildGRCh38_cleaned.tsv.gz", show_col_types=FALSE) %>%
  mutate(
    chromosome         = as.character(chromosome),
    base_pair_location = as.integer(base_pair_location)
  ) %>%
  transmute(
    chromosome, base_pair_location,
    ea_rep2 = effect_allele,
    ra_rep2 = other_allele,
    es_rep2 = beta,
    se_rep2 = standard_error,
    p_rep2  = p_value
  )

# beta as es
rep3 <- read_tsv("GCST90476065_cleaned.tsv.gz", show_col_types=FALSE) %>%
  mutate(
    chromosome         = as.character(chromosome),
    base_pair_location = as.integer(base_pair_location)
  ) %>%
  transmute(
    chromosome, base_pair_location,
    ea_rep3 = effect_allele,
    ra_rep3 = other_allele,
    es_rep3 = log(odds_ratio),
    se_rep3 = standard_error,
    p_rep3  = p_value
  )
# beta as es
rep4 <- read_tsv("GCST90476067_cleaned.tsv.gz", show_col_types=FALSE) %>%
  mutate(
    chromosome         = as.character(chromosome),
    base_pair_location = as.integer(base_pair_location)
  ) %>%
  transmute(
    chromosome, base_pair_location,
    ea_rep4 = effect_allele,
    ra_rep4 = other_allele,
    es_rep4 = log(odds_ratio),
    se_rep4 = standard_error,
    p_rep4  = p_value
  )

# Harmonize each replicate dataset with the discovery dataset
harm1 <- disc %>%
  inner_join(rep1, by = c("chromosome", "base_pair_location")) %>%
  mutate(
    flip1   = ea_rep1 != ea_disc,
    es1_adj = if_else(flip1, -es_rep1, es_rep1),
    ea1_adj = if_else(flip1, ra_rep1, ea_rep1),
    ra1_adj = if_else(flip1, ea_rep1, ra_rep1)
  )

harm2 <- disc %>%
  inner_join(rep2, by = c("chromosome", "base_pair_location")) %>%
  mutate(
    flip2   = ea_rep2 != ea_disc,
    es2_adj = if_else(flip2, -es_rep2, es_rep2),
    ea2_adj = if_else(flip2, ra_rep2, ea_rep2),
    ra2_adj = if_else(flip2, ea_rep2, ra_rep2)
  )

harm3 <- disc %>%
  inner_join(rep3, by = c("chromosome", "base_pair_location")) %>%
  mutate(
    flip3   = ea_rep3 != ea_disc,
    es3_adj = if_else(flip3, -es_rep3, es_rep3),
    ea3_adj = if_else(flip3, ra_rep3, ea_rep3),
    ra3_adj = if_else(flip3, ea_rep3, ra_rep3)
  )

harm4 <- disc %>%
  inner_join(rep4, by = c("chromosome", "base_pair_location")) %>%
  mutate(
    flip4   = ea_rep4 != ea_disc,
    es4_adj = if_else(flip4, -es_rep4, es_rep4),
    ea4_adj = if_else(flip4, ra_rep4, ea_rep4),
    ra4_adj = if_else(flip4, ea_rep4, ra_rep4)
  )

# Merge
harm <- harm1 %>%
  select(chromosome, base_pair_location,
         ea_disc, ra_disc, es_disc, se_disc, p_disc,
         ea_rep1, ra_rep1, es_rep1, se_rep1, p_rep1,
         flip1, ea1_adj, ra1_adj, es1_adj) %>%
  inner_join(
    harm2 %>% select(chromosome, base_pair_location,
                     ea_rep2, ra_rep2, es_rep2, se_rep2, p_rep2,
                     flip2, ea2_adj, ra2_adj, es2_adj),
    by = c("chromosome", "base_pair_location")
  ) %>%
  inner_join(
    harm3 %>% select(chromosome, base_pair_location,
                     ea_rep3, ra_rep3, es_rep3, se_rep3, p_rep3,
                     flip3, ea3_adj, ra3_adj, es3_adj),
    by = c("chromosome", "base_pair_location")
  ) %>%
  inner_join(
    harm4 %>% select(chromosome, base_pair_location,
                     ea_rep4, ra_rep4, es_rep4, se_rep4, p_rep4,
                     flip4, ea4_adj, ra4_adj, es4_adj),
    by = c("chromosome", "base_pair_location")
  )

harm %>%
  transmute(
    SNP = paste0(chromosome, ":", base_pair_location),
    ea_disc, ra_disc, es_disc, se_disc, p_disc,
    ea_rep1, ra_rep1, es_rep1, se_rep1, p_rep1, flip1, ea1_adj, ra1_adj, es1_adj,
    ea_rep2, ra_rep2, es_rep2, se_rep2, p_rep2, flip2, ea2_adj, ra2_adj, es2_adj,
    ea_rep3, ra_rep3, es_rep3, se_rep3, p_rep3, flip3, ea3_adj, ra3_adj, es3_adj,
    ea_rep4, ra_rep4, es_rep4, se_rep4, p_rep4, flip4, ea4_adj, ra4_adj, es4_adj
  ) %>%
  kable(
    caption = glue::glue("Allele Harmonization for Top {N} SNPs across 4 Replication Sets"),
    digits = 3
  )


```


