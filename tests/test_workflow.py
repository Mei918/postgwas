from pathlib import Path
import gzip

from gwas_postgwas_tools.config import load_config
import gwas_postgwas_tools.workflows.postgwas as workflow_module


def test_workflow_groups_by_gene_and_cell_type(tmp_path):
    config = load_config(Path("examples/ad_ebv_config.yaml"))
    config.output_dir = str(tmp_path)
    config.batch_chunk_size = 2
    config.batch_progress_every = 1

    original_coloc = workflow_module.run_coloc_abf
    original_mr = workflow_module.run_single_variant_mr

    def fake_coloc(**kwargs):
        class Result:
            def __init__(self):
                self.gene_id = kwargs["gene_id"]
                self.cell_type = kwargs["cell_type"]
                self.n_overlap_snps = len(kwargs["harmonized_df"])
                self.pp4 = 0.9
                self.passed_threshold = True

        return Result()

    def fake_mr(**kwargs):
        class Result:
            def __init__(self):
                self.gene_id = kwargs["gene_id"]
                self.cell_type = kwargs["cell_type"]
                self.method = kwargs["method"]
                self.beta = 0.1
                self.se = 0.01
                self.p_value = 0.001

        return Result()

    workflow_module.run_coloc_abf = fake_coloc
    workflow_module.run_single_variant_mr = fake_mr

    try:
        outputs = workflow_module.run_postgwas_pipeline(config)
    finally:
        workflow_module.run_coloc_abf = original_coloc
        workflow_module.run_single_variant_mr = original_mr

    assert outputs["coloc"].exists()
    assert outputs["mr"].exists()
    assert outputs["group_summary"].exists()
    assert outputs["failures"].exists()
    assert outputs["coloc_hits"].exists()
    assert outputs["mr_hits"].exists()
    assert outputs["priority_hits"].exists()


def test_iter_eqtl_groups_handles_group_boundary_across_chunks(tmp_path):
    eqtl_path = tmp_path / "eqtl.tsv.gz"
    with gzip.open(eqtl_path, "wt") as handle:
        handle.write(
            "chrom\tposition\teffect_allele\tother_allele\tbeta\tse\tp_value\tgene_id\tcell_type\n"
            "1\t100\tA\tG\t0.1\t0.01\t0.001\tGENE1\tmicroglia\n"
            "1\t101\tA\tG\t0.2\t0.02\t0.002\tGENE1\tmicroglia\n"
            "1\t102\tA\tG\t0.3\t0.03\t0.003\tGENE1\tmicroglia\n"
            "1\t103\tA\tG\t0.4\t0.04\t0.004\tGENE2\tastrocytes\n"
        )

    groups = list(
        workflow_module.iter_eqtl_groups(
            spec=load_config(Path("examples/ad_ebv_config.yaml")).eqtl.__class__(
                path=str(eqtl_path),
                sep="\t",
                compression="gzip",
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
            ),
            chunk_size=2,
        )
    )

    assert len(groups) == 2
    assert groups[0][0] == "GENE1"
    assert groups[0][1] == "microglia"
    assert len(groups[0][2]) == 3


def test_workflow_respects_limit_groups(tmp_path):
    config = load_config(Path("examples/ad_ebv_config.yaml"))
    config.output_dir = str(tmp_path)
    config.batch_chunk_size = 2
    config.batch_limit_groups = 1
    config.batch_progress_every = 1

    original_coloc = workflow_module.run_coloc_abf
    original_mr = workflow_module.run_single_variant_mr

    def fake_coloc(**kwargs):
        class Result:
            def __init__(self):
                self.gene_id = kwargs["gene_id"]
                self.cell_type = kwargs["cell_type"]
                self.n_overlap_snps = len(kwargs["harmonized_df"])
                self.pp4 = 0.9
                self.passed_threshold = True
        return Result()

    def fake_mr(**kwargs):
        class Result:
            def __init__(self):
                self.gene_id = kwargs["gene_id"]
                self.cell_type = kwargs["cell_type"]
                self.method = kwargs["method"]
                self.beta = 0.1
                self.se = 0.01
                self.p_value = 0.001
        return Result()

    workflow_module.run_coloc_abf = fake_coloc
    workflow_module.run_single_variant_mr = fake_mr

    try:
        outputs = workflow_module.run_postgwas_pipeline(config)
    finally:
        workflow_module.run_coloc_abf = original_coloc
        workflow_module.run_single_variant_mr = original_mr

    group_summary = Path(outputs["group_summary"]).read_text().strip().splitlines()
    assert len(group_summary) == 2


def test_workflow_writes_hits_summaries(tmp_path):
    config = load_config(Path("examples/ad_ebv_config.yaml"))
    config.output_dir = str(tmp_path)
    config.batch_chunk_size = 2
    config.batch_limit_groups = 2
    config.batch_progress_every = 1

    original_coloc = workflow_module.run_coloc_abf
    original_mr = workflow_module.run_single_variant_mr

    def fake_coloc(**kwargs):
        class Result:
            def __init__(self):
                self.gene_id = kwargs["gene_id"]
                self.cell_type = kwargs["cell_type"]
                self.n_overlap_snps = len(kwargs["harmonized_df"])
                self.pp4 = 0.9
                self.passed_threshold = True
        return Result()

    def fake_mr(**kwargs):
        class Result:
            def __init__(self):
                self.gene_id = kwargs["gene_id"]
                self.cell_type = kwargs["cell_type"]
                self.method = kwargs["method"]
                self.beta = 0.1
                self.se = 0.01
                self.p_value = 0.001
        return Result()

    workflow_module.run_coloc_abf = fake_coloc
    workflow_module.run_single_variant_mr = fake_mr

    try:
        outputs = workflow_module.run_postgwas_pipeline(config)
    finally:
        workflow_module.run_coloc_abf = original_coloc
        workflow_module.run_single_variant_mr = original_mr

    coloc_hits = Path(outputs["coloc_hits"]).read_text().strip().splitlines()
    mr_hits = Path(outputs["mr_hits"]).read_text().strip().splitlines()
    priority_hits = Path(outputs["priority_hits"]).read_text().strip().splitlines()

    assert len(coloc_hits) >= 2
    assert len(mr_hits) >= 2
    assert len(priority_hits) >= 2


def test_workflow_prefilter_skips_weak_groups(tmp_path):
    config = load_config(Path("examples/ad_ebv_config.yaml"))
    config.output_dir = str(tmp_path)
    config.batch_chunk_size = 2
    config.batch_limit_groups = 1
    config.batch_progress_every = 1
    config.batch_min_harmonized_variants = 10

    original_coloc = workflow_module.run_coloc_abf
    original_mr = workflow_module.run_single_variant_mr

    def fail_if_called(**kwargs):
        raise AssertionError("coloc/MR should not run for skipped groups")

    workflow_module.run_coloc_abf = fail_if_called
    workflow_module.run_single_variant_mr = fail_if_called

    try:
        outputs = workflow_module.run_postgwas_pipeline(config)
    finally:
        workflow_module.run_coloc_abf = original_coloc
        workflow_module.run_single_variant_mr = original_mr

    group_summary = Path(outputs["group_summary"]).read_text().strip().splitlines()
    coloc_lines = Path(outputs["coloc"]).read_text().strip().splitlines()
    mr_lines = Path(outputs["mr"]).read_text().strip().splitlines()

    assert len(group_summary) == 2
    assert "skipped" in group_summary[1]
    assert "low_harmonized_variants" in group_summary[1]
    assert len(coloc_lines) == 2
    assert len(mr_lines) == 2
