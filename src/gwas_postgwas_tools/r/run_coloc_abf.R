args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 14) {
  stop(
    paste(
      "Expected 14 arguments:",
      "input_path output_path gene_id cell_type eqtl_trait_type gwas_trait_type",
      "pp4_threshold p1 p2 p12 eqtl_sd_y eqtl_sample_size gwas_sample_size gwas_case_prop"
    )
  )
}

input_path <- args[[1]]
output_path <- args[[2]]
gene_id <- args[[3]]
cell_type <- args[[4]]
eqtl_trait_type <- args[[5]]
gwas_trait_type <- args[[6]]
pp4_threshold <- as.numeric(args[[7]])
p1 <- as.numeric(args[[8]])
p2 <- as.numeric(args[[9]])
p12 <- as.numeric(args[[10]])
eqtl_sd_y <- if (args[[11]] == "NA") NA_real_ else as.numeric(args[[11]])
eqtl_sample_size <- if (args[[12]] == "NA") NA_real_ else as.numeric(args[[12]])
gwas_sample_size <- if (args[[13]] == "NA") NA_real_ else as.numeric(args[[13]])
gwas_case_prop <- if (args[[14]] == "NA") NA_real_ else as.numeric(args[[14]])

suppressPackageStartupMessages({
  library(coloc)
})

input_df <- read.delim(input_path, sep = "\t", header = TRUE, check.names = FALSE)

if (nrow(input_df) == 0) {
  stop("No harmonized variants available for coloc.")
}

required_columns <- c(
  "variant_id",
  "beta_eqtl",
  "se_eqtl",
  "beta_gwas",
  "gwas_beta_aligned",
  "se_gwas"
)
missing_columns <- setdiff(required_columns, colnames(input_df))
if (length(missing_columns) > 0) {
  stop(paste("Missing required harmonized columns:", paste(missing_columns, collapse = ", ")))
}

build_eqtl_dataset <- function(df, trait_type, sd_y, sample_size) {
  dataset <- list(
    snp = df$variant_id,
    beta = df$beta_eqtl,
    varbeta = (df$se_eqtl)^2,
    type = trait_type
  )

  if (!is.na(sd_y)) {
    dataset$sdY <- sd_y
  }
  if (!is.na(sample_size)) {
    dataset$N <- sample_size
  }
  dataset
}

build_gwas_dataset <- function(df, trait_type, sample_size, case_prop) {
  dataset <- list(
    snp = df$variant_id,
    beta = df$gwas_beta_aligned,
    varbeta = (df$se_gwas)^2,
    type = trait_type
  )

  if (!is.na(sample_size)) {
    dataset$N <- sample_size
  }
  if (trait_type == "cc" && !is.na(case_prop)) {
    dataset$s <- case_prop
  }
  dataset
}

eqtl_dataset <- build_eqtl_dataset(input_df, eqtl_trait_type, eqtl_sd_y, eqtl_sample_size)
gwas_dataset <- build_gwas_dataset(input_df, gwas_trait_type, gwas_sample_size, gwas_case_prop)

coloc_fit <- coloc.abf(
  dataset1 = eqtl_dataset,
  dataset2 = gwas_dataset,
  p1 = p1,
  p2 = p2,
  p12 = p12
)

pp4 <- unname(coloc_fit$summary[["PP.H4.abf"]])
n_overlap_snps <- nrow(input_df)

output_df <- data.frame(
  gene_id = gene_id,
  cell_type = cell_type,
  n_overlap_snps = n_overlap_snps,
  pp4 = pp4,
  passed_threshold = pp4 >= pp4_threshold
)

write.table(output_df, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
