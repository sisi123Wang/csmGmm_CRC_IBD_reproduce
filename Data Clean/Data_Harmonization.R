# - Please change the seetwd path to your lolcal path.
setwd("/home/swang25/other_summary_statistics/")

library(tibble)
library(dplyr)
library(readr)
library(stats)
library(rlang)

na_strings <- c("", "NA", "#NA", "N/A", "null", "NULL")

# cleaned dataset list

files_info <- tribble(
  ~type,        ~id,            ~file_path,                                                                                 ~sample_size,
  "colorectal", "GCST90013862", "/home/swang25/other_summary_statistics/colorectal_cancer/GCST90013862_buildGRCh38_cleaned.tsv.gz", 407746,
  "colorectal", "GCST90255675", "/home/swang25/other_summary_statistics/colorectal_cancer/GCST90255675_buildGRCh38_cleaned.tsv.gz", 185616,
  "colorectal", "GCST90475566", "/home/swang25/other_summary_statistics/colorectal_cancer/GCST90475566_buildGRCh38_cleaned.tsv.gz", 446766,
  "IBD",        "GCST90044155", "/home/swang25/other_summary_statistics/IBD/GCST90044155_buildGRCh38_cleaned.tsv.gz",              456348,
  "IBD",        "GCST90081487", "/home/swang25/other_summary_statistics/IBD/GCST90081487_buildGRCh38_cleaned.tsv.gz",              364036,
  "IBD",        "GCST90292538", "/home/swang25/other_summary_statistics/IBD/GCST90292538_cleaned.tsv.gz",                          398668,
  "IBD",        "GCST90476065", "/home/swang25/other_summary_statistics/IBD/GCST90476065_cleaned.tsv.gz",                          446753,
  "IBD",        "GCST90476067", "/home/swang25/other_summary_statistics/IBD/GCST90476067_cleaned.tsv.gz",                          447331,
  "obesity",    "GCST90081522", "/home/swang25/other_summary_statistics/obesity/GCST90081522_buildGRCh38_cleaned.tsv.gz",          246639,
  "obesity",    "GCST90475762", "/home/swang25/other_summary_statistics/obesity/GCST90475762_cleaned.tsv.gz",                      413854
)

# Harmonization function
harmonize_two_traits <- function(df1, df2,
                                 chr = "chromosome", pos = "base_pair_location",
                                 ea1 = "ea_trait1", ra1 = "ra_trait1", es1 = "es_trait1",
                                 ea2 = "ea_trait2", ra2 = "ra_trait2", es2 = "es_trait2") {

  merged <- inner_join(df1, df2, by = c(chr, pos)) %>%
    mutate(
      across(c(!!ea1, !!ra1, !!ea2, !!ra2), ~ toupper(.x)),
      es2_original = !!sym(es2)  # Save original
    )

  harmonized <- merged %>%
    mutate(
      flip = case_when(
        !!sym(ea1) == !!sym(ra2) & !!sym(ra1) == !!sym(ea2) ~ TRUE,
        !!sym(ea1) == !!sym(ea2) & !!sym(ra1) == !!sym(ra2) ~ FALSE,
        TRUE ~ NA
      )
    ) %>%
    filter(!is.na(flip)) %>%
    mutate(
      # Flip effect size
      es_trait2 = if_else(flip, -!!sym(es2), !!sym(es2)),
      # Use temp variables to avoid overwriting in-place
      temp_ea = if_else(flip, !!sym(ra2), !!sym(ea2)),
      temp_ra = if_else(flip, !!sym(ea2), !!sym(ra2))
    ) %>%
    mutate(
      !!ea2 := temp_ea,
      !!ra2 := temp_ra
    ) %>%
    select(-temp_ea, -temp_ra) %>%
    relocate(es2_original, .before = es_trait2)

  return(harmonized)
}

# CRC: load and rename with SE imputation, Z
df_crc_raw <- read_tsv(files_info$file_path[files_info$type == "colorectal"][2],
                       show_col_types = FALSE,
                       na = na_strings)

df_crc <- df_crc_raw %>%
  rename(
    es_crc = beta,
    se_crc = standard_error,
    p_crc  = p_value,
    ea_crc = effect_allele,
    ra_crc = other_allele
  ) %>%
  mutate(
    # SE values are missing
    se_crc_imputed = is.na(se_crc),

    # SE = |beta| / |qnorm(p/2)|
    se_crc = case_when(
      !is.na(se_crc) ~ se_crc,  # Keep existing SE
      is.na(se_crc) & !is.na(es_crc) & !is.na(p_crc) & p_crc > 0 & p_crc < 1 ~
        abs(es_crc) / abs(qnorm(p_crc / 2)),  # Impute SE
      TRUE ~ NA_real_  # Keep as NA if cannot impute
    ),

    # Calculate Z-score = beta / SE
    z_crc = case_when(
      !is.na(es_crc) & !is.na(se_crc) & se_crc != 0 ~ es_crc / se_crc,
      TRUE ~ NA_real_
    )
  )

cat(" SE missing:", sum(df_crc$se_crc_imputed, na.rm = TRUE), "\n")
cat("SE after imputation:", sum(!is.na(df_crc$se_crc)), "\n")
cat("Z-scores calculated:", sum(!is.na(df_crc$z_crc)), "\n\n")

# IBD: load + rename + SE imputation + Z
df_ibd_raw <- read_tsv(files_info$file_path[files_info$type == "IBD"][1],
                       show_col_types = FALSE,
                       na = na_strings)

df_ibd <- df_ibd_raw %>%
  rename(
    es_ibd = beta,
    se_ibd = standard_error,
    p_ibd  = p_value,
    ea_ibd = effect_allele,
    ra_ibd = other_allele
  ) %>%
  mutate(
    se_ibd_imputed = is.na(se_ibd),
    se_ibd = case_when(
      !is.na(se_ibd) ~ se_ibd,
      is.na(se_ibd) & !is.na(es_ibd) & !is.na(p_ibd) & p_ibd > 0 & p_ibd < 1 ~
        abs(es_ibd) / abs(qnorm(p_ibd / 2)),
      TRUE ~ NA_real_
    ),
    z_ibd = case_when(
      !is.na(es_ibd) & !is.na(se_ibd) & se_ibd != 0 ~ es_ibd / se_ibd,
      TRUE ~ NA_real_
    )
  )

cat("Total SNPs:", nrow(df_ibd), "\n")
cat("SE originally missing:", sum(df_ibd$se_ibd_imputed, na.rm = TRUE), "\n")
cat("Z-scores calculated:", sum(!is.na(df_ibd$z_ibd)), "\n\n")


# Obesity: load + beta conversion if needed + SE imputation + Z
ob_path <- files_info$file_path[files_info$type == "obesity"][1]
df_ob_raw <- read_tsv(ob_path, show_col_types = FALSE, na = na_strings)

df_ob <- df_ob_raw %>%
  # Convert odds ratio to beta
  {
    if (!"beta" %in% names(.) && "odds_ratio" %in% names(.)) {
      mutate(., beta = log(odds_ratio))
    } else {
      .
    }
  } %>%
  # Clean and standardize
  mutate(
    chromosome = as.character(chromosome),
    base_pair_location = as.integer(base_pair_location),
    effect_allele = toupper(effect_allele),
    other_allele = toupper(other_allele)
  ) %>%
  rename(
    es_ob = beta,
    se_ob = standard_error,
    p_ob  = p_value,
    ea_ob = effect_allele,
    ra_ob = other_allele
  ) %>%
  mutate(
    se_ob_imputed = is.na(se_ob),
    se_ob = case_when(
      !is.na(se_ob) ~ se_ob,
      is.na(se_ob) & !is.na(es_ob) & !is.na(p_ob) & p_ob > 0 & p_ob < 1 ~
        abs(es_ob) / abs(qnorm(p_ob / 2)),
      TRUE ~ NA_real_
    ),
    z_ob = case_when(
      !is.na(es_ob) & !is.na(se_ob) & se_ob != 0 ~ es_ob / se_ob,
      TRUE ~ NA_real_
    )
  )

cat("Obesity data loaded:\n")
cat("- Total SNPs:", nrow(df_ob), "\n")
cat("- SE originally missing:", sum(df_ob$se_ob_imputed, na.rm = TRUE), "\n")
cat("- SE after imputation (available):", sum(!is.na(df_ob$se_ob)), "\n")
cat("- Z-scores calculated:", sum(!is.na(df_ob$z_ob)), "\n\n")

# Prep for harmonization (rename to common trait1/trait2 schema)
df_crc_tmp <- df_crc %>%
  mutate(chromosome = as.character(chromosome)) %>%
  rename(
    ea_trait1 = ea_crc,
    ra_trait1 = ra_crc,
    es_trait1 = es_crc,
    se_trait1 = se_crc
  )

df_ibd_tmp <- df_ibd %>%
  mutate(chromosome = as.character(chromosome)) %>%
  rename(
    ea_trait2 = ea_ibd,
    ra_trait2 = ra_ibd,
    es_trait2 = es_ibd,
    se_trait2 = se_ibd
  )

df_ob_tmp <- df_ob %>%
  mutate(chromosome = as.character(chromosome)) %>%
  rename(
    ea_trait2 = ea_ob,
    ra_trait2 = ra_ob,
    es_trait2 = es_ob,
    se_trait2 = se_ob
  )

# Harmonize CRC and IBD use this first. if not working, then use harmonize_two_traits_final function below
crc_ob_harmonized <- harmonize_two_traits(df_crc_tmp, df_ibd_tmp)
cat("CRC-Obesity harmonized SNPs:", nrow(crc_ob_harmonized), "\n")
cat("Initial dataset sizes:\n")
cat("- CRC SNPs:", nrow(df_crc_tmp), "\n")
cat("- IBD SNPs:", nrow(df_ibd_tmp), "\n\n")

step1_merged <- inner_join(df_crc_tmp, df_ibd_tmp, by = c("chromosome", "base_pair_location"))
cat("After inner join (same chr+pos):", nrow(step1_merged), "\n")

step2_alleles <- step1_merged %>%
  mutate(across(all_of(c("ea_trait1", "ra_trait1", "ea_trait2", "ra_trait2")), toupper))

step3_harmonized <- step2_alleles %>%
  mutate(
    is_flip = case_when(
      ea_trait1 == ra_trait2 & ra_trait1 == ea_trait2 ~ TRUE,   # flip needed
      ea_trait1 == ea_trait2 & ra_trait1 == ra_trait2 ~ FALSE,  # same direction
      TRUE                                            ~ NA      # allele mismatch
    )
  )

# cat("Allele mismatch (dropped):", sum(is.na(step3_harmonized$is_flip)), "\n")


# Make output folder exists, common out once the folder exist
dir.create("Data", showWarnings = FALSE, recursive = TRUE)
# Harmonize CRC and IBD
crc_ibd_harmonized <- harmonize_two_traits(df_crc_tmp, df_ibd_tmp)
write_tsv(crc_ibd_harmonized, "Data/crc_ibd.txt")
readr::write_tsv(crc_ibd_harmonized, "Data/crc_ibd_harmonized.txt")

# BACKUP harmonization function, if function harmonize_two_traits works, no need to run below.
# ─────────────────────────────────────────────────────────────────────────────
harmonize_two_traits_final <- function(df1, df2,
                                       chr = "chromosome",
                                       pos = "base_pair_location",
                                       ea1 = "ea_trait1", ra1 = "ra_trait1", es1 = "es_trait1",
                                       ea2 = "ea_trait2", ra2 = "ra_trait2", es2 = "es_trait2") {

  if (!all(c(chr, pos, ea1, ra1, es1) %in% names(df1))) stop("Missing required columns in df1")
  if (!all(c(chr, pos, ea2, ra2, es2) %in% names(df2))) stop("Missing required columns in df2")

  cat("=== Harmonization Process ===\n")
  cat("Dataset 1 SNPs:", nrow(df1), "\n")
  cat("Dataset 2 SNPs:", nrow(df2), "\n")

  merged <- inner_join(df1, df2, by = c(chr, pos))
  cat("After position matching:", nrow(merged), "SNPs\n")

  harmonized <- merged %>%
    mutate(across(all_of(c(ea1, ra1, ea2, ra2)), toupper)) %>%
    mutate(
      is_flip = case_when(
        .data[[ea1]] == .data[[ra2]] & .data[[ra1]] == .data[[ea2]] ~ TRUE,
        .data[[ea1]] == .data[[ea2]] & .data[[ra1]] == .data[[ra2]] ~ FALSE,
        TRUE ~ NA
      )
    )

  n_same     <- sum(harmonized$is_flip == FALSE, na.rm = TRUE)
  n_flip     <- sum(harmonized$is_flip == TRUE,  na.rm = TRUE)
  n_mismatch <- sum(is.na(harmonized$is_flip))

  cat("Allele matching results:\n")
  cat("  Same direction:", n_same, "SNPs\n")
  cat("  Need flipping:", n_flip, "SNPs\n")
  cat("  Mismatched (dropped):", n_mismatch, "SNPs\n")

  harmonized <- harmonized %>%
    filter(!is.na(is_flip)) %>%
    mutate(
      ea2_aligned  = if_else(is_flip, .data[[ra2]], .data[[ea2]]),
      ra2_aligned  = if_else(is_flip, .data[[ea2]], .data[[ra2]]),
      es2_aligned  = if_else(is_flip, -.data[[es2]], .data[[es2]]),
      es1_original = .data[[es1]],
      es2_original = .data[[es2]],
      flipped      = is_flip
    )

  cat("Final harmonized SNPs:", nrow(harmonized), "\n")
  cat("Harmonization complete!\n\n")

  return(harmonized)
}

#save outputs
crc_ibd_harmonized <- harmonize_two_traits_final(
  df_crc_tmp, df_ibd_tmp,
  ea1 = "ea_trait1", ra1 = "ra_trait1", es1 = "es_trait1",
  ea2 = "ea_trait2", ra2 = "ra_trait2", es2 = "es_trait2"
)
write_tsv(crc_ibd_harmonized, "Data/crc_ibd_final.txt")


