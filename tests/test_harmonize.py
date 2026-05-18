from gwas_postgwas_tools.io.sumstats import read_sumstats
from gwas_postgwas_tools.preprocess.harmonize import (
    harmonize_variants,
    prepare_gwas_index,
    subset_gwas_for_eqtl,
)


def test_harmonize_variants_flips_gwas_beta_when_alleles_are_reversed():
    eqtl = read_sumstats("examples/demo_eqtl.tsv")
    gwas = read_sumstats("examples/demo_gwas.tsv")

    result = harmonize_variants(eqtl, gwas)
    flipped = result.loc[result["variant_id"] == "1:1010"].iloc[0]

    assert flipped["gwas_beta_aligned"] == -0.06


def test_subset_gwas_for_eqtl_limits_to_requested_variants():
    eqtl = read_sumstats("examples/demo_eqtl.tsv")
    gwas = read_sumstats("examples/demo_gwas.tsv")
    gwas_indexed = prepare_gwas_index(gwas)

    subset = subset_gwas_for_eqtl(eqtl.iloc[[0]], gwas_indexed)

    assert len(subset) == 1
    assert subset.iloc[0]["position"] == 1000
