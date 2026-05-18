import math

import pandas as pd
import pytest

from gwas_postgwas_tools.tools.mr import run_single_variant_mr


def test_run_single_variant_mr_returns_default_for_empty_input():
    result = run_single_variant_mr(
        harmonized_df=pd.DataFrame(),
        gene_id="GENE1",
        cell_type="microglia",
    )

    assert result.gene_id == "GENE1"
    assert result.cell_type == "microglia"
    assert result.method == "wald_ratio"
    assert result.beta == 0.0
    assert result.se == 0.0
    assert result.p_value == 1.0


def test_run_single_variant_mr_computes_python_wald_ratio():
    harmonized_df = pd.DataFrame(
        [
            {
                "beta_eqtl": 0.2,
                "se_eqtl": 0.05,
                "gwas_beta_aligned": 0.1,
                "se_gwas": 0.02,
                "p_value_eqtl": 0.01,
            },
            {
                "beta_eqtl": 0.3,
                "se_eqtl": 0.04,
                "gwas_beta_aligned": 0.12,
                "se_gwas": 0.03,
                "p_value_eqtl": 0.02,
            },
        ]
    )

    result = run_single_variant_mr(
        harmonized_df=harmonized_df,
        gene_id="GENE1",
        cell_type="microglia",
    )

    expected_beta = 0.1 / 0.2
    expected_se = math.sqrt((0.02 ** 2) / (0.2 ** 2) + ((0.1 ** 2) * (0.05 ** 2)) / (0.2 ** 4))

    assert result.gene_id == "GENE1"
    assert result.cell_type == "microglia"
    assert result.method == "wald_ratio"
    assert result.beta == pytest.approx(expected_beta)
    assert result.se == pytest.approx(expected_se)
    assert 0.0 <= result.p_value <= 1.0


def test_run_single_variant_mr_rejects_unsupported_method():
    harmonized_df = pd.DataFrame(
        [
            {
                "beta_eqtl": 0.2,
                "se_eqtl": 0.05,
                "gwas_beta_aligned": 0.1,
                "se_gwas": 0.02,
                "p_value_eqtl": 0.01,
            }
        ]
    )

    with pytest.raises(ValueError, match="Unsupported MR method"):
        run_single_variant_mr(
            harmonized_df=harmonized_df,
            gene_id="GENE1",
            cell_type="microglia",
            method="ivw",
        )
