library(dplyr)
library(tidyr)
library(readr)
library(knitr)
library(scales)
library(purrr)
library(glue)

data_dir <- "/home/swang25/other_summary_statistics/colorectal_cancer"

files <- tibble::tribble(
  ~id,            ~path,                                                                              ~build,   ~sample_size,
  "GCST90013862", file.path(data_dir, "GCST90013862_buildGRCh38.tsv.gz"), "GRCh38",   407746,
  "GCST90255675", file.path(data_dir, "GCST90255675_buildGRCh37.tsv"),    "GRCh37",   185616,
  "GCST90475566", file.path(data_dir, "GCST90475566_buildGRCh38.tsv.gz"), "GRCh38",   446766
)

kable(files, caption = "Colorectal GWAS Files Metadata")

# preview and QC summaries
library(data.table); library(dplyr); library(knitr); library(tidyr); library(scales)

#can adjust based on your path
setwd("~/other_summary_statistics/colorectal_cancer/")

datasets <- tibble::tribble(
  ~id,            ~file,                                                                                 ~sample_size,
  "GCST90013862", "/home/swang25/other_summary_statistics/colorectal_cancer/GCST90013862_buildGRCh38.tsv.gz", 407746,
  "GCST90255675", "/home/swang25/other_summary_statistics/colorectal_cancer/GCST90255675_buildGRCh37.tsv",    185616,
  "GCST90475566", "/home/swang25/other_summary_statistics/colorectal_cancer/GCST90475566_buildGRCh38.tsv.gz", 446766
)

library(data.table)
library(dplyr)


# dataset head check
for (i in seq_len(nrow(datasets))) {
  ds <- datasets[i, ]
  cat("\n===== Dataset:", ds$id, "=====\n")
  cat("File:", ds$file, "\n")

  # first 5 rows for a fast preview
  dt_head <- fread(ds$file, nrows = 5)

  cat("Column Names:\n")
  print(names(dt_head))
  cat("\nHead:\n")
  print(dt_head)
}


for (i in seq_len(nrow(datasets))) {
  ds  <- datasets[i, ]
  cat("Dataset:", ds$id, "\n\n")
  cat("File:", ds$file, "\n\n")

  raw <- fread(ds$file)
  cat("Original dimensions:", nrow(raw), "rows ×", ncol(raw), "columns\n\n")
  cat("Columns:\n")
  print(names(raw)); cat("\n")

  df <- raw %>% select(chromosome, base_pair_location, effect_allele, other_allele, p_value)
  cat("Preview:\n")
  print(kable(head(df,5), format="markdown")); cat("\n")

  clean <- raw %>%
    filter(!is.na(chromosome), !is.na(base_pair_location), !is.na(p_value))
  cat("Removed", nrow(raw)-nrow(clean), "rows with missing", nrow(clean), "remain\n\n")
  cat("Sample size (N):", ds$sample_size, "\n\n")

  # SNP counts
  snp_summary <- tibble(`Total SNPs`=nrow(df), `Unique SNPs`=n_distinct(df$base_pair_location))

  # Unique SNPs ignores chromosome.
  # `Unique SNPs` = n_distinct(paste(df$chromosome, df$base_pair_location))

  cat("SNP counts:\n")
  print(kable(snp_summary, format="markdown")); cat("\n")

  # p-value thresholds
  pval_counts <- df %>%
    summarise(
      `<5e-8`=sum(p_value<5e-8,na.rm=TRUE),
      `<5e-6`=sum(p_value<5e-6,na.rm=TRUE),
      `<5e-4`=sum(p_value<5e-4,na.rm=TRUE),
      `<0.05`=sum(p_value<0.05,na.rm=TRUE)
    ) %>%
    pivot_longer(everything(),names_to="Threshold",values_to="Count")

  cat("P-value thresholds:\n")
  print(kable(pval_counts, format="markdown")); cat("\n\n---\n\n")
}

# Cleaning 
#can adjust path
setwd("/home/swang25/other_summary_statistics/colorectal_cancer/")

# here::i_am("Colorectal_data_clean.Rmd")

clean_gwas <- function(id, path, build) {
  cat("Cleaning", id, "(assembly:", build, ")\n\n")

  raw <- read_tsv(path, show_col_types = FALSE)
  n0  <- nrow(raw)
  cat("Original rows:", n0, "* Columns:*", paste(names(raw), collapse=", "), "\n\n")

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

  # remove multi base
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

  # include alleles to avoid collapsing
  # distinct(chromosome, base_pair_location, effect_allele, other_allele, .keep_all = TRUE)

  n2 <- nrow(df)
  cat("After deduplication:", n2, "(removed", n1-n2, "dups)\n\n")

  # odd ratio and beta
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


# LiftOver GRCh37 -> GRCh38 (GCST90255675)

setwd("~/other_summary_statistics/colorectal_cancer/")

library(readr)
library(dplyr)


# assemble GRCh37-cleaned file then convert it to GRCh38-cleaned file
df37 <- read_tsv("GCST90255675_buildGRCh37.tsv",
                 show_col_types = FALSE)

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

# liftOver (make sure! ./liftOver and hg19ToHg38.over.chain.gz are in this directory)
system("./liftOver coords.bed hg19ToHg38.over.chain.gz mapped.bed unmapped.bed")

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
  mutate(.idx = row_number()) %>%
  inner_join(mapped, by = c(".idx" = "idx")) %>%
  mutate(
    chromosome = sub("^chr", "", chrom),
    base_pair_location = end1,
    assembly = "GRCh38"
  ) %>%
  select(-.idx, -chrom, -start0, -end1)

out38 <- "GCST90255675_buildGRCh38_cleaned.tsv.gz"
readr::write_tsv(df38, gzfile(out38))

cat("Lifted", nrow(df38), "of", nrow(df37),
    "variants to GRCh38 in file:", out38, "\n")


setwd("~/other_summary_statistics/colorectal_cancer/")


disc <- read_tsv("GCST90013862_buildGRCh38_cleaned.tsv.gz", show_col_types = FALSE)
cat("GCST90013862 columns:\n"); print(names(disc))
cat("\nGCST90013862 data preview:\n"); print(head(disc, 5))

rep1 <- read_tsv("GCST90475566_buildGRCh38_cleaned.tsv.gz", show_col_types = FALSE)
cat("GCST90475566 columns:\n"); print(names(rep1))
cat("\nGCST90475566 data preview:\n"); print(head(rep1, 5))

rep2 <- read_tsv("GCST90255675_buildGRCh38_cleaned.tsv.gz", show_col_types = FALSE)
cat("GCST90255675 columns:\n"); print(names(rep2))
cat("\nGCST90255675 data preview:\n"); print(head(rep2, 5))

# Allele harmonization QC (Top N)

setwd("~/other_summary_statistics/colorectal_cancer/")

library(readr)
library(dplyr)
library(knitr)
library(glue)

N <- 10

disc <- read_tsv("GCST90013862_buildGRCh38_cleaned.tsv.gz", show_col_types=FALSE) %>%
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
    es_disc   = beta,
    se_disc   = standard_error,
    p_disc    = p_value
  )

# log(odds_ratio) as es
rep1 <- read_tsv("GCST90475566_buildGRCh38_cleaned.tsv.gz", show_col_types=FALSE) %>%
  mutate(
    chromosome         = as.character(chromosome),
    base_pair_location = as.integer(base_pair_location)
  ) %>%
  transmute(
    chromosome, base_pair_location,
    ea_rep1 = effect_allele,
    ra_rep1 = other_allele,
    es_rep1 = log(odds_ratio),
    se_rep1 = standard_error,
    p_rep1  = p_value
  )

# beta as es
rep2 <- read_tsv("GCST90255675_buildGRCh38_cleaned.tsv.gz", show_col_types=FALSE) %>%
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


harm1 <- disc %>%
  inner_join(rep1, by = c("chromosome","base_pair_location")) %>%
  mutate(
    flip1      = (ea_rep1 != ea_disc),
    es1_adj    = if_else(flip1, -es_rep1, es_rep1),
    ea1_adj    = if_else(flip1, ra_rep1, ea_rep1),
    ra1_adj    = if_else(flip1, ea_rep1, ra_rep1)
  )

harm2 <- disc %>%
  inner_join(rep2, by = c("chromosome","base_pair_location")) %>%
  mutate(
    flip2    = ea_rep2 != ea_disc,
    es2_adj  = if_else(flip2, -es_rep2, es_rep2),
    ea2_adj  = if_else(flip2, ra_rep2, ea_rep2),
    ra2_adj  = if_else(flip2, ea_rep2, ra_rep2)
  )

harm <- harm1 %>%
  select(chromosome, base_pair_location,
         ea_disc, ra_disc, es_disc, se_disc, p_disc,
         ea_rep1, ra_rep1, es_rep1, se_rep1, p_rep1,
         flip1, ea1_adj, ra1_adj, es1_adj) %>%
  inner_join(
    harm2 %>%
      select(chromosome, base_pair_location,
             ea_rep2, ra_rep2, es_rep2, se_rep2, p_rep2,
             flip2, ea2_adj, ra2_adj, es2_adj),
    by = c("chromosome","base_pair_location")
  )

harm %>%
  transmute(
    SNP              = paste0(chromosome, ":", base_pair_location),
    ea_disc, ra_disc, es_disc, se_disc, p_disc,
    ea_rep1, ra_rep1, es_rep1, se_rep1, p_rep1, flip1, ea1_adj, ra1_adj, es1_adj,
    ea_rep2, ra_rep2, es_rep2, se_rep2, p_rep2, flip2, ea2_adj, ra2_adj, es2_adj
  ) %>%
  kable(
    caption = glue("Allele Harmonization for Top {N} SNPs"),
    digits  = 3
  )

