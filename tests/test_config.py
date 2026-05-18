from pathlib import Path

from gwas_postgwas_tools.config import load_config


def test_load_config_supports_nested_input_specs():
    config = load_config(Path("examples/ad_ebv_config.yaml"))

    assert config.eqtl.columns["gene_id"] == "gene_id"
    assert config.gwas.columns["beta"] == "beta"
    assert config.batch_mode == "gene_cell_type"
    assert config.batch_chunk_size == 200000
    assert config.batch_progress_every == 10
    assert config.batch_min_harmonized_variants is None
    assert config.eqtl_pvalue_threshold is None
    assert config.gwas_pvalue_threshold is None
