args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop("Expected 5 arguments: input_path output_path gene_id cell_type method")
}

input_path <- args[[1]]
output_path <- args[[2]]
gene_id <- args[[3]]
cell_type <- args[[4]]
method <- args[[5]]

suppressPackageStartupMessages({
  library(MendelianRandomization)
})

input_df <- read.delim(input_path, sep = "\t", header = TRUE, check.names = FALSE)

if (nrow(input_df) == 0) {
  stop("No harmonized variants available for MR.")
}

required_columns <- c("beta_eqtl", "se_eqtl", "gwas_beta_aligned", "se_gwas", "p_value_eqtl")
missing_columns <- setdiff(required_columns, colnames(input_df))
if (length(missing_columns) > 0) {
  stop(paste("Missing required harmonized columns:", paste(missing_columns, collapse = ", ")))
}

if (method != "wald_ratio") {
  stop(paste("Unsupported MR method:", method))
}

lead_row <- input_df[order(input_df$p_value_eqtl, decreasing = FALSE), ][1, ]
mr_input <- mr_input(
  bx = lead_row$beta_eqtl,
  bxse = lead_row$se_eqtl,
  by = lead_row$gwas_beta_aligned,
  byse = lead_row$se_gwas
)
mr_fit <- mr_wald_ratio(mr_input)

output_df <- data.frame(
  gene_id = gene_id,
  cell_type = cell_type,
  method = method,
  beta = mr_fit@Estimate,
  se = mr_fit@StdError,
  p_value = mr_fit@Pvalue
)

write.table(output_df, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
