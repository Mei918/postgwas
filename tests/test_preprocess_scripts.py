import pandas as pd

from scripts.build_standardized_eqtl import list_eqtl_files, standardize_eqtl, stream_standardized_eqtl
from scripts.build_standardized_gwas import standardize_gwas


def test_standardize_eqtl_derives_required_columns():
    eqtl_df = pd.DataFrame(
        {
            "gene_id": ["GENE1"],
            "snp_id": ["rs1"],
            "dist_tss": [1000],
            "p_value": [0.01],
            "beta": [0.2],
            "cell_type": ["microglia"],
        }
    )
    snp_map_df = pd.DataFrame(
        {
            "SNP": ["rs1"],
            "SNP_id_hg38": ["chr1:12345"],
            "effect_allele": ["A"],
            "other_allele": ["G"],
            "MAF": [0.2],
        }
    )

    standardized = standardize_eqtl(eqtl_df, snp_map_df)

    assert standardized.loc[0, "chrom"] == "1"
    assert standardized.loc[0, "position"] == 12345
    assert standardized.loc[0, "effect_allele"] == "A"
    assert standardized.loc[0, "cell_type"] == "microglia"


def test_stream_standardized_eqtl_writes_output_incrementally(tmp_path):
    eqtl_root = tmp_path / "eqtl"
    astro_dir = eqtl_root / "astrocytes"
    astro_dir.mkdir(parents=True)
    (astro_dir / "Astrocytes.1").write_text(
        "GENE1 rs1 1000 0.01 0.2\nGENE2 rs2 2000 0.02 -0.3\n"
    )

    snp_map_df = pd.DataFrame(
        {
            "SNP": ["rs1", "rs2"],
            "SNP_id_hg38": ["chr1:12345", "chr1:22222"],
            "effect_allele": ["A", "C"],
            "other_allele": ["G", "T"],
            "MAF": [0.2, 0.3],
        }
    )
    output_path = tmp_path / "out.tsv.gz"

    stats = stream_standardized_eqtl(
        eqtl_root=eqtl_root,
        snp_map_df=snp_map_df,
        output_path=output_path,
        selected_cell_types=["astrocytes"],
        chunk_size=1,
    )

    written = pd.read_csv(output_path, sep="\t", compression="gzip")
    assert stats["rows"] == 2
    assert stats["genes"] == 2
    assert stats["cell_types"] == 1
    assert len(written) == 2
    assert set(written["gene_id"]) == {"GENE1", "GENE2"}


def test_standardize_gwas_converts_neg_log_pvalue():
    gwas_df = pd.DataFrame(
        {
            "SNP_id_hg38": ["chr2:100"],
            "beta": [0.1],
            "stderr_beta": [0.02],
            "ref": ["C"],
            "alt": ["T"],
            "MAF": [0.4],
            "neg_log_pvalue": [3.0],
        }
    )

    standardized = standardize_gwas(gwas_df)

    assert standardized.loc[0, "chrom"] == "2"
    assert standardized.loc[0, "position"] == 100
    assert round(standardized.loc[0, "p_value"], 6) == 0.001
