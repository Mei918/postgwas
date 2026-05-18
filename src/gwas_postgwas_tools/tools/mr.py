from __future__ import annotations

from dataclasses import asdict
from statistics import NormalDist
from typing import Optional, Union

import pandas as pd

from gwas_postgwas_tools.models.result import MRResult


def run_single_variant_mr(
    harmonized_df: pd.DataFrame,
    gene_id: str,
    cell_type: str,
    rscript_executable: str = "Rscript",
    method: str = "wald_ratio",
    working_dir: Optional[Union[str, Path]] = None,
) -> MRResult:
    if harmonized_df.empty:
        return MRResult(
            gene_id=gene_id,
            cell_type=cell_type,
            method=method,
            beta=0.0,
            se=0.0,
            p_value=1.0,
        )

    if method != "wald_ratio":
        raise ValueError(f"Unsupported MR method: {method}")

    required_columns = {"beta_eqtl", "se_eqtl", "gwas_beta_aligned", "se_gwas", "p_value_eqtl"}
    missing_columns = required_columns.difference(harmonized_df.columns)
    if missing_columns:
        raise ValueError(
            "Missing required harmonized columns: " + ", ".join(sorted(missing_columns))
        )

    lead_row = (
        harmonized_df.sort_values("p_value_eqtl", ascending=True)
        .iloc[0]
    )

    bx = float(lead_row["beta_eqtl"])
    bxse = float(lead_row["se_eqtl"])
    by = float(lead_row["gwas_beta_aligned"])
    byse = float(lead_row["se_gwas"])

    if bx == 0.0:
        raise ValueError("Lead eQTL beta is zero; cannot compute Wald ratio.")

    beta = by / bx
    se = ((byse ** 2) / (bx ** 2) + ((by ** 2) * (bxse ** 2)) / (bx ** 4)) ** 0.5
    z_score = 0.0 if se == 0.0 else beta / se
    p_value = 2 * NormalDist().cdf(-abs(z_score))

    return MRResult(
        gene_id=gene_id,
        cell_type=cell_type,
        method=method,
        beta=float(beta),
        se=float(se),
        p_value=float(p_value),
    )


def mr_result_to_dict(result: MRResult) -> dict[str, str | float]:
    return asdict(result)
