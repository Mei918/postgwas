from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a standardized GWAS table from a mapped GWAS summary file."
    )
    parser.add_argument("--gwas", required=True, help="Path to mapped GWAS file.")
    parser.add_argument("--output", required=True, help="Output TSV or TSV.GZ path.")
    return parser


def main() -> int:
    args = build_parser().parse_args()

    gwas_path = Path(args.gwas).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()

    gwas_df = pd.read_csv(gwas_path, sep=r"\s+")
    standardized = standardize_gwas(gwas_df)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if output_path.suffix == ".gz" else None
    standardized.to_csv(output_path, sep="\t", index=False, compression=compression)

    print(f"rows\t{len(standardized)}")
    print(f"output\t{output_path}")
    return 0


def standardize_gwas(gwas_df: pd.DataFrame) -> pd.DataFrame:
    required = {
        "SNP_id_hg38",
        "beta",
        "stderr_beta",
        "ref",
        "alt",
        "MAF",
    }
    missing = required.difference(gwas_df.columns)
    if missing:
        raise ValueError(f"GWAS input is missing required columns: {sorted(missing)}")

    coord = gwas_df["SNP_id_hg38"].astype(str).str.replace("^chr", "", regex=True).str.split(":", expand=True)
    standardized = pd.DataFrame(
        {
            "chrom": coord[0],
            "position": coord[1].astype(int),
            "effect_allele": gwas_df["alt"].astype(str).str.upper(),
            "other_allele": gwas_df["ref"].astype(str).str.upper(),
            "beta": gwas_df["beta"].astype(float),
            "se": gwas_df["stderr_beta"].astype(float),
            "maf": gwas_df["MAF"].astype(float),
            "snp_id_hg38": gwas_df["SNP_id_hg38"].astype(str),
        }
    )

    if "p" in gwas_df.columns:
        standardized["p_value"] = gwas_df["p"].astype(float)
    elif "neg_log_pvalue" in gwas_df.columns:
        standardized["p_value"] = 10 ** (-gwas_df["neg_log_pvalue"].astype(float))
    else:
        raise ValueError("GWAS input must contain either 'p' or 'neg_log_pvalue'.")

    if "rsid" in gwas_df.columns:
        standardized["rsid"] = gwas_df["rsid"].astype(str)

    return standardized.sort_values(["chrom", "position"]).reset_index(drop=True)


if __name__ == "__main__":
    raise SystemExit(main())
