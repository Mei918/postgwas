from __future__ import annotations

from dataclasses import asdict
from importlib import resources
from pathlib import Path
import tempfile
from typing import Optional, Union

import pandas as pd

from gwas_postgwas_tools.models.result import ColocResult
from gwas_postgwas_tools.tools.base import ExternalTool


def run_coloc_abf(
    harmonized_df: pd.DataFrame,
    gene_id: str,
    cell_type: str,
    rscript_executable: str = "Rscript",
    eqtl_trait_type: str = "quant",
    gwas_trait_type: str = "quant",
    min_overlap_snps: int = 20,
    pp4_threshold: float = 0.8,
    p1: float = 1e-4,
    p2: float = 1e-4,
    p12: float = 1e-5,
    eqtl_sd_y: Optional[float] = None,
    eqtl_sample_size: Optional[int] = None,
    gwas_sample_size: Optional[int] = None,
    gwas_case_prop: Optional[float] = None,
    working_dir: Optional[Union[str, Path]] = None,
) -> ColocResult:
    n_overlap = len(harmonized_df)
    if n_overlap < min_overlap_snps:
        return ColocResult(
            gene_id=gene_id,
            cell_type=cell_type,
            n_overlap_snps=n_overlap,
            pp4=0.0,
            passed_threshold=False,
        )

    script_path = _resolve_r_script_path("run_coloc_abf.R")
    tool = ExternalTool(executable=rscript_executable)

    temp_dir_root = Path(working_dir).resolve() if working_dir else None
    with tempfile.TemporaryDirectory(dir=temp_dir_root) as temp_dir:
        temp_path = Path(temp_dir).resolve()
        input_path = temp_path / "harmonized.tsv"
        output_path = temp_path / "coloc.tsv"
        harmonized_df.to_csv(input_path, sep="\t", index=False)

        args = [
            str(script_path),
            str(input_path),
            str(output_path),
            gene_id,
            cell_type,
            eqtl_trait_type,
            gwas_trait_type,
            str(pp4_threshold),
            str(p1),
            str(p2),
            str(p12),
            _optional_arg(eqtl_sd_y),
            _optional_arg(eqtl_sample_size),
            _optional_arg(gwas_sample_size),
            _optional_arg(gwas_case_prop),
        ]

        result = tool.run(args=args, cwd=temp_dir_root)
        if result.returncode != 0:
            raise RuntimeError(
                "coloc R wrapper failed with exit code "
                f"{result.returncode}: {result.stderr.strip() or result.stdout.strip()}"
            )

        output_df = pd.read_csv(output_path, sep="\t")
        row = output_df.iloc[0]
        return ColocResult(
            gene_id=str(row["gene_id"]),
            cell_type=str(row["cell_type"]),
            n_overlap_snps=int(row["n_overlap_snps"]),
            pp4=float(row["pp4"]),
            passed_threshold=bool(row["passed_threshold"]),
        )


def coloc_result_to_dict(result: ColocResult) -> dict[str, str | int | float | bool]:
    return asdict(result)


def _resolve_r_script_path(script_name: str) -> Path:
    package_root = resources.files("gwas_postgwas_tools")
    return Path(str(package_root.joinpath("r", script_name)))


def _optional_arg(value: object) -> str:
    return "NA" if value is None else str(value)
