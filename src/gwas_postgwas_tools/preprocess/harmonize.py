from __future__ import annotations

import pandas as pd


AMBIGUOUS_PAIRS = {
    ("A", "T"),
    ("T", "A"),
    ("C", "G"),
    ("G", "C"),
}


def _normalize_alleles(df: pd.DataFrame) -> pd.DataFrame:
    normalized = df.copy()
    normalized["chrom"] = normalized["chrom"].astype(str).str.replace("^chr", "", regex=True)
    normalized["effect_allele"] = normalized["effect_allele"].str.upper()
    normalized["other_allele"] = normalized["other_allele"].str.upper()
    normalized["variant_id"] = (
        normalized["chrom"].astype(str) + ":" + normalized["position"].astype(str)
    )
    return normalized


def prepare_gwas_index(gwas_df: pd.DataFrame) -> pd.DataFrame:
    gwas = _normalize_alleles(gwas_df)
    return gwas.drop_duplicates(subset=["variant_id"]).set_index("variant_id", drop=False)


def subset_gwas_for_eqtl(eqtl_df: pd.DataFrame, gwas_indexed_df: pd.DataFrame) -> pd.DataFrame:
    eqtl = _normalize_alleles(eqtl_df)
    variant_ids = eqtl["variant_id"].drop_duplicates()
    matched = gwas_indexed_df.reindex(variant_ids)
    matched = matched.dropna(subset=["variant_id"]).reset_index(drop=True)
    return matched


def harmonize_variants(
    eqtl_df: pd.DataFrame,
    gwas_df: pd.DataFrame | None = None,
    gwas_subset_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    eqtl = _normalize_alleles(eqtl_df)
    if gwas_subset_df is not None:
        gwas = _normalize_alleles(gwas_subset_df)
    elif gwas_df is not None:
        gwas = _normalize_alleles(gwas_df)
    else:
        raise ValueError("harmonize_variants requires either gwas_df or gwas_subset_df.")

    merged = eqtl.merge(
        gwas,
        on="variant_id",
        suffixes=("_eqtl", "_gwas"),
        how="inner",
    )

    same_orientation = (
        (merged["effect_allele_eqtl"] == merged["effect_allele_gwas"])
        & (merged["other_allele_eqtl"] == merged["other_allele_gwas"])
    )
    flipped_orientation = (
        (merged["effect_allele_eqtl"] == merged["other_allele_gwas"])
        & (merged["other_allele_eqtl"] == merged["effect_allele_gwas"])
    )

    aligned = merged.loc[same_orientation | flipped_orientation].copy()
    aligned["gwas_beta_aligned"] = aligned["beta_gwas"]
    aligned.loc[flipped_orientation.loc[aligned.index], "gwas_beta_aligned"] *= -1

    ambiguous = aligned.apply(
        lambda row: (row["effect_allele_eqtl"], row["other_allele_eqtl"]) in AMBIGUOUS_PAIRS,
        axis=1,
    )
    return aligned.loc[~ambiguous].reset_index(drop=True)
