from gwas_postgwas_tools.io.sumstats import read_sumstats
from gwas_postgwas_tools.config import InputSpec


def test_read_sumstats_maps_source_columns_to_standard_names():
    spec = InputSpec(
        path="examples/demo_eqtl.tsv",
        columns={
            "chrom": "chrom",
            "position": "position",
            "effect_allele": "effect_allele",
            "other_allele": "other_allele",
            "beta": "beta",
            "se": "se",
            "p_value": "p_value",
            "gene_id": "gene_id",
            "cell_type": "cell_type",
        },
    )

    df = read_sumstats(spec)

    assert "gene_id" in df.columns
    assert "cell_type" in df.columns
