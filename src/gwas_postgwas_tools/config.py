from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Union

import yaml


STANDARD_COLUMN_NAMES = {
    "chrom",
    "position",
    "effect_allele",
    "other_allele",
    "beta",
    "se",
    "p_value",
    "maf",
    "gene_id",
    "cell_type",
}

REQUIRED_EQTL_COLUMNS = {
    "chrom",
    "position",
    "effect_allele",
    "other_allele",
    "beta",
    "se",
    "p_value",
    "gene_id",
    "cell_type",
}

REQUIRED_GWAS_COLUMNS = {
    "chrom",
    "position",
    "effect_allele",
    "other_allele",
    "beta",
    "se",
    "p_value",
}


@dataclass
class InputSpec:
    path: str
    sep: str = "\t"
    compression: Optional[str] = None
    columns: dict[str, str] = field(default_factory=dict)


@dataclass
class StudyConfig:
    name: str
    genome_build: str
    eqtl: InputSpec
    gwas: InputSpec
    output_dir: str
    rscript_executable: str = "Rscript"
    coloc_min_overlap_snps: int = 20
    coloc_pp4_threshold: float = 0.8
    coloc_p1: float = 1e-4
    coloc_p2: float = 1e-4
    coloc_p12: float = 1e-5
    mr_method: str = "wald_ratio"
    eqtl_trait_type: str = "quant"
    gwas_trait_type: str = "quant"
    eqtl_sd_y: Optional[float] = None
    eqtl_sample_size: Optional[int] = None
    gwas_sample_size: Optional[int] = None
    gwas_case_prop: Optional[float] = None
    batch_mode: str = "gene_cell_type"
    write_group_harmonized: bool = False
    batch_chunk_size: int = 200_000
    batch_limit_groups: Optional[int] = None
    batch_progress_every: int = 10
    batch_min_harmonized_variants: Optional[int] = None
    eqtl_pvalue_threshold: Optional[float] = None
    gwas_pvalue_threshold: Optional[float] = None


def load_config(path: Union[str, Path]) -> StudyConfig:
    config_path = Path(path)
    payload = yaml.safe_load(config_path.read_text()) or {}

    study = payload.get("study", {})
    reference = payload.get("reference", {})
    inputs = payload.get("inputs", {})
    postgwas = payload.get("postgwas", {})
    coloc = postgwas.get("coloc", {})
    mr = postgwas.get("mr", {})
    runtime = payload.get("runtime", {})
    batch = payload.get("batch", {})

    eqtl_spec = _parse_input_spec(inputs, "eqtl")
    gwas_spec = _parse_input_spec(inputs, "gwas")

    config = StudyConfig(
        name=study["name"],
        genome_build=reference["genome_build"],
        eqtl=eqtl_spec,
        gwas=gwas_spec,
        output_dir=inputs["output_dir"],
        rscript_executable=runtime.get("rscript_executable", "Rscript"),
        coloc_min_overlap_snps=coloc.get("min_overlap_snps", 20),
        coloc_pp4_threshold=coloc.get("pp4_threshold", 0.8),
        coloc_p1=coloc.get("p1", 1e-4),
        coloc_p2=coloc.get("p2", 1e-4),
        coloc_p12=coloc.get("p12", 1e-5),
        mr_method=mr.get("method", "wald_ratio"),
        eqtl_trait_type=postgwas.get("eqtl_trait_type", "quant"),
        gwas_trait_type=postgwas.get("gwas_trait_type", postgwas.get("trait_type", "quant")),
        eqtl_sd_y=postgwas.get("eqtl_sd_y"),
        eqtl_sample_size=postgwas.get("eqtl_sample_size"),
        gwas_sample_size=postgwas.get("gwas_sample_size"),
        gwas_case_prop=postgwas.get("gwas_case_prop"),
        batch_mode=batch.get("mode", "gene_cell_type"),
        write_group_harmonized=batch.get("write_group_harmonized", False),
        batch_chunk_size=batch.get("chunk_size", 200_000),
        batch_limit_groups=batch.get("limit_groups"),
        batch_progress_every=batch.get("progress_every", 10),
        batch_min_harmonized_variants=batch.get("min_harmonized_variants"),
        eqtl_pvalue_threshold=postgwas.get("eqtl_pvalue_threshold"),
        gwas_pvalue_threshold=postgwas.get("gwas_pvalue_threshold"),
    )
    validate_config(config)
    return config


def validate_config(config: StudyConfig) -> None:
    if config.genome_build not in {"GRCh37", "GRCh38"}:
        raise ValueError(
            f"Unsupported genome_build '{config.genome_build}'. Expected GRCh37 or GRCh38."
        )

    if config.eqtl_trait_type not in {"quant"}:
        raise ValueError("eqtl_trait_type must currently be 'quant'.")

    if config.gwas_trait_type not in {"quant", "cc"}:
        raise ValueError("gwas_trait_type must be 'quant' or 'cc'.")

    if config.batch_mode not in {"gene_cell_type"}:
        raise ValueError("batch.mode must currently be 'gene_cell_type'.")

    if config.batch_chunk_size < 1:
        raise ValueError("batch.chunk_size must be >= 1.")

    if config.batch_limit_groups is not None and config.batch_limit_groups < 1:
        raise ValueError("batch.limit_groups must be >= 1 when provided.")

    if config.batch_progress_every < 1:
        raise ValueError("batch.progress_every must be >= 1.")

    if config.batch_min_harmonized_variants is not None and config.batch_min_harmonized_variants < 1:
        raise ValueError("batch.min_harmonized_variants must be >= 1 when provided.")

    if config.coloc_min_overlap_snps < 1:
        raise ValueError("postgwas.coloc.min_overlap_snps must be >= 1.")

    if not 0 <= config.coloc_pp4_threshold <= 1:
        raise ValueError("postgwas.coloc.pp4_threshold must be between 0 and 1.")

    if config.gwas_trait_type == "cc" and config.gwas_case_prop is None:
        raise ValueError("postgwas.gwas_case_prop is required when gwas_trait_type is 'cc'.")

    if config.eqtl_trait_type == "quant" and config.eqtl_sd_y is None:
        raise ValueError("postgwas.eqtl_sd_y is required for quant eQTL coloc input.")

    if config.eqtl_sample_size is not None and config.eqtl_sample_size <= 0:
        raise ValueError("postgwas.eqtl_sample_size must be > 0 when provided.")

    if config.gwas_sample_size is not None and config.gwas_sample_size <= 0:
        raise ValueError("postgwas.gwas_sample_size must be > 0 when provided.")

    if config.gwas_case_prop is not None and not 0 < config.gwas_case_prop < 1:
        raise ValueError("postgwas.gwas_case_prop must be between 0 and 1.")

    if config.eqtl_pvalue_threshold is not None and not 0 < config.eqtl_pvalue_threshold <= 1:
        raise ValueError("postgwas.eqtl_pvalue_threshold must be between 0 and 1.")

    if config.gwas_pvalue_threshold is not None and not 0 < config.gwas_pvalue_threshold <= 1:
        raise ValueError("postgwas.gwas_pvalue_threshold must be between 0 and 1.")

    _validate_input_spec(config.eqtl, REQUIRED_EQTL_COLUMNS, "inputs.eqtl")
    _validate_input_spec(config.gwas, REQUIRED_GWAS_COLUMNS, "inputs.gwas")


def config_to_dict(config: StudyConfig) -> dict[str, Any]:
    return {
        "study": {
            "name": config.name,
        },
        "reference": {
            "genome_build": config.genome_build,
        },
        "inputs": {
            "eqtl": _input_spec_to_dict(config.eqtl),
            "gwas": _input_spec_to_dict(config.gwas),
            "output_dir": config.output_dir,
        },
        "runtime": {
            "rscript_executable": config.rscript_executable,
        },
        "batch": {
            "mode": config.batch_mode,
            "write_group_harmonized": config.write_group_harmonized,
            "chunk_size": config.batch_chunk_size,
            "limit_groups": config.batch_limit_groups,
            "progress_every": config.batch_progress_every,
            "min_harmonized_variants": config.batch_min_harmonized_variants,
        },
        "postgwas": {
            "eqtl_trait_type": config.eqtl_trait_type,
            "gwas_trait_type": config.gwas_trait_type,
            "eqtl_sd_y": config.eqtl_sd_y,
            "eqtl_sample_size": config.eqtl_sample_size,
            "gwas_sample_size": config.gwas_sample_size,
            "gwas_case_prop": config.gwas_case_prop,
            "eqtl_pvalue_threshold": config.eqtl_pvalue_threshold,
            "gwas_pvalue_threshold": config.gwas_pvalue_threshold,
            "coloc": {
                "min_overlap_snps": config.coloc_min_overlap_snps,
                "pp4_threshold": config.coloc_pp4_threshold,
                "p1": config.coloc_p1,
                "p2": config.coloc_p2,
                "p12": config.coloc_p12,
            },
            "mr": {
                "method": config.mr_method,
            },
        },
    }


def _parse_input_spec(inputs: dict[str, Any], key: str) -> InputSpec:
    nested = inputs.get(key)
    if isinstance(nested, dict):
        return InputSpec(
            path=nested["path"],
            sep=nested.get("sep", "\t"),
            compression=nested.get("compression"),
            columns=nested.get("columns", {}),
        )

    legacy_path_key = f"{key}_path"
    if legacy_path_key in inputs:
        return InputSpec(path=inputs[legacy_path_key])

    raise ValueError(f"Missing required input section: inputs.{key}")


def _validate_input_spec(spec: InputSpec, required_columns: set[str], label: str) -> None:
    unknown_columns = set(spec.columns).difference(STANDARD_COLUMN_NAMES)
    if unknown_columns:
        raise ValueError(
            f"{label}.columns contains unsupported standard names: {sorted(unknown_columns)}"
        )

    missing_mappings = required_columns.difference(set(spec.columns))
    if missing_mappings:
        raise ValueError(
            f"{label}.columns is missing required mappings for: {sorted(missing_mappings)}"
        )


def _input_spec_to_dict(spec: InputSpec) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "path": spec.path,
        "sep": spec.sep,
        "columns": spec.columns,
    }
    if spec.compression is not None:
        payload["compression"] = spec.compression
    return payload
