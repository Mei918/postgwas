from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Iterator, Optional

import pandas as pd

from gwas_postgwas_tools.config import InputSpec, StudyConfig
from gwas_postgwas_tools.io.sumstats import read_sumstats
from gwas_postgwas_tools.models.result import ColocResult, MRResult
from gwas_postgwas_tools.preprocess.harmonize import (
    harmonize_variants,
    prepare_gwas_index,
    subset_gwas_for_eqtl,
)
from gwas_postgwas_tools.tools.coloc import run_coloc_abf
from gwas_postgwas_tools.tools.mr import run_single_variant_mr


def run_postgwas_pipeline(config: StudyConfig) -> dict[str, Path]:
    output_dir = Path(config.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    gwas_df = read_sumstats(config.gwas)
    gwas_indexed_df = prepare_gwas_index(gwas_df)
    del gwas_df

    coloc_rows = []
    mr_rows = []
    group_summary_rows = []
    failure_rows = []
    harmonized_paths = []

    processed_groups = 0
    for gene_id, cell_type, eqtl_group_df in iter_eqtl_groups(config.eqtl, chunk_size=config.batch_chunk_size):
        processed_groups += 1
        if config.batch_limit_groups is not None and processed_groups > config.batch_limit_groups:
            break

        try:
            gwas_subset_df = subset_gwas_for_eqtl(eqtl_df=eqtl_group_df, gwas_indexed_df=gwas_indexed_df)
            harmonized_df = harmonize_variants(eqtl_df=eqtl_group_df, gwas_subset_df=gwas_subset_df)

            skip_reason = _prefilter_skip_reason(config=config, harmonized_df=harmonized_df)
            if skip_reason is not None:
                coloc_rows.append(
                    _result_to_record(
                        ColocResult(
                            gene_id=gene_id,
                            cell_type=cell_type,
                            n_overlap_snps=len(harmonized_df),
                            pp4=0.0,
                            passed_threshold=False,
                        )
                    )
                )
                mr_rows.append(
                    _result_to_record(
                        MRResult(
                            gene_id=gene_id,
                            cell_type=cell_type,
                            method=config.mr_method,
                            beta=0.0,
                            se=0.0,
                            p_value=1.0,
                        )
                    )
                )
                group_summary_rows.append(
                    {
                        "gene_id": gene_id,
                        "cell_type": cell_type,
                        "eqtl_variants": len(eqtl_group_df),
                        "harmonized_variants": len(harmonized_df),
                        "status": "skipped",
                        "skip_reason": skip_reason,
                    }
                )
                continue

            if config.write_group_harmonized:
                harmonized_path = output_dir / f"{config.name}.{gene_id}.{cell_type}.harmonized.tsv"
                harmonized_df.to_csv(harmonized_path, sep="\t", index=False)
                harmonized_paths.append(harmonized_path)

            coloc_result = run_coloc_abf(
                harmonized_df=harmonized_df,
                gene_id=gene_id,
                cell_type=cell_type,
                rscript_executable=config.rscript_executable,
                eqtl_trait_type=config.eqtl_trait_type,
                gwas_trait_type=config.gwas_trait_type,
                min_overlap_snps=config.coloc_min_overlap_snps,
                pp4_threshold=config.coloc_pp4_threshold,
                p1=config.coloc_p1,
                p2=config.coloc_p2,
                p12=config.coloc_p12,
                eqtl_sd_y=config.eqtl_sd_y,
                eqtl_sample_size=config.eqtl_sample_size,
                gwas_sample_size=config.gwas_sample_size,
                gwas_case_prop=config.gwas_case_prop,
                working_dir=output_dir,
            )
            mr_result = run_single_variant_mr(
                harmonized_df=harmonized_df,
                gene_id=gene_id,
                cell_type=cell_type,
                rscript_executable=config.rscript_executable,
                method=config.mr_method,
                working_dir=output_dir,
            )

            coloc_rows.append(_result_to_record(coloc_result))
            mr_rows.append(_result_to_record(mr_result))
            group_summary_rows.append(
                {
                    "gene_id": gene_id,
                    "cell_type": cell_type,
                    "eqtl_variants": len(eqtl_group_df),
                    "harmonized_variants": len(harmonized_df),
                    "status": "ok",
                    "skip_reason": "",
                }
            )

            if processed_groups % config.batch_progress_every == 0:
                print(
                    f"[progress] groups={processed_groups} gene={gene_id} "
                    f"cell_type={cell_type} harmonized={len(harmonized_df)}",
                    flush=True,
                )
        except Exception as exc:
            failure_rows.append(
                {
                    "gene_id": gene_id,
                    "cell_type": cell_type,
                    "eqtl_variants": len(eqtl_group_df),
                    "error": str(exc),
                }
            )
            group_summary_rows.append(
                {
                    "gene_id": gene_id,
                    "cell_type": cell_type,
                    "eqtl_variants": len(eqtl_group_df),
                    "harmonized_variants": -1,
                    "status": "failed",
                    "skip_reason": "",
                }
            )
            print(
                f"[failure] groups={processed_groups} gene={gene_id} "
                f"cell_type={cell_type} error={exc}",
                flush=True,
            )

    coloc_df = _safe_sort_frame(pd.DataFrame(coloc_rows), ["gene_id", "cell_type"])
    mr_df = _safe_sort_frame(pd.DataFrame(mr_rows), ["gene_id", "cell_type"])
    summary_df = _safe_sort_frame(pd.DataFrame(group_summary_rows), ["gene_id", "cell_type"])
    failures_df = _safe_sort_frame(pd.DataFrame(failure_rows), ["gene_id", "cell_type"])

    coloc_path = output_dir / f"{config.name}.coloc.tsv"
    mr_path = output_dir / f"{config.name}.mr.tsv"
    summary_path = output_dir / f"{config.name}.groups.tsv"
    failures_path = output_dir / f"{config.name}.failures.tsv"
    coloc_df.to_csv(coloc_path, sep="\t", index=False)
    mr_df.to_csv(mr_path, sep="\t", index=False)
    summary_df.to_csv(summary_path, sep="\t", index=False)
    failures_df.to_csv(failures_path, sep="\t", index=False)

    coloc_hits_df = coloc_df[coloc_df["passed_threshold"] == True].copy() if not coloc_df.empty else coloc_df.copy()
    mr_hits_df = mr_df[mr_df["p_value"] < 0.05].copy() if not mr_df.empty else mr_df.copy()
    priority_hits_df = (
        coloc_hits_df.merge(
            mr_hits_df,
            on=["gene_id", "cell_type"],
            suffixes=("_coloc", "_mr"),
            how="inner",
        )
        if not coloc_hits_df.empty and not mr_hits_df.empty
        else pd.DataFrame()
    )

    coloc_hits_path = output_dir / f"{config.name}.coloc_hits.tsv"
    mr_hits_path = output_dir / f"{config.name}.mr_hits.tsv"
    priority_hits_path = output_dir / f"{config.name}.priority_hits.tsv"
    coloc_hits_df.to_csv(coloc_hits_path, sep="\t", index=False)
    mr_hits_df.to_csv(mr_hits_path, sep="\t", index=False)
    priority_hits_df.to_csv(priority_hits_path, sep="\t", index=False)

    outputs = {
        "group_summary": summary_path,
        "coloc": coloc_path,
        "mr": mr_path,
        "failures": failures_path,
        "coloc_hits": coloc_hits_path,
        "mr_hits": mr_hits_path,
        "priority_hits": priority_hits_path,
    }
    if harmonized_paths:
        outputs["harmonized_dir"] = output_dir
    return outputs


def iter_eqtl_groups(spec: InputSpec, chunk_size: int) -> Iterator[tuple[str, str, pd.DataFrame]]:
    path = Path(spec.path)
    compression = spec.compression
    if compression is None and path.suffix == ".gz":
        compression = "gzip"

    eqtl_chunks = pd.read_csv(
        path,
        sep=spec.sep,
        compression=compression,
        chunksize=chunk_size,
    )

    carryover: Optional[pd.DataFrame] = None
    for chunk in eqtl_chunks:
        if carryover is not None and not carryover.empty:
            chunk = pd.concat([carryover, chunk], ignore_index=True)
            carryover = None

        if chunk.empty:
            continue

        chunk["gene_id"] = chunk["gene_id"].astype(str)
        chunk["cell_type"] = chunk["cell_type"].astype(str)
        group_keys = list(zip(chunk["gene_id"], chunk["cell_type"]))
        if not group_keys:
            continue

        last_key = group_keys[-1]
        split_index = len(group_keys)
        for index in range(len(group_keys) - 2, -1, -1):
            if group_keys[index] != last_key:
                split_index = index + 1
                break
            if index == 0:
                split_index = 0

        complete = chunk.iloc[:split_index].copy() if split_index > 0 else pd.DataFrame(columns=chunk.columns)
        carryover = chunk.iloc[split_index:].copy()

        if not complete.empty:
            yield from _yield_grouped_frame(complete)

    if carryover is not None and not carryover.empty:
        yield from _yield_grouped_frame(carryover)


def _yield_grouped_frame(df: pd.DataFrame) -> Iterator[tuple[str, str, pd.DataFrame]]:
    for (gene_id, cell_type), group_df in df.groupby(["gene_id", "cell_type"], sort=False):
        yield str(gene_id), str(cell_type), group_df.reset_index(drop=True)


def _result_to_record(result):
    if hasattr(result, "__dataclass_fields__"):
        return asdict(result)
    return dict(vars(result))


def _safe_sort_frame(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    if df.empty:
        return df
    return df.sort_values(columns).reset_index(drop=True)


def _prefilter_skip_reason(config: StudyConfig, harmonized_df: pd.DataFrame) -> str | None:
    if harmonized_df.empty:
        return "empty_harmonized"

    if (
        config.batch_min_harmonized_variants is not None
        and len(harmonized_df) < config.batch_min_harmonized_variants
    ):
        return "low_harmonized_variants"

    if config.eqtl_pvalue_threshold is not None:
        min_eqtl_p = pd.to_numeric(harmonized_df["p_value_eqtl"], errors="coerce").min()
        if pd.isna(min_eqtl_p) or float(min_eqtl_p) >= config.eqtl_pvalue_threshold:
            return "weak_eqtl_signal"

    if config.gwas_pvalue_threshold is not None:
        min_gwas_p = pd.to_numeric(harmonized_df["p_value_gwas"], errors="coerce").min()
        if pd.isna(min_gwas_p) or float(min_gwas_p) >= config.gwas_pvalue_threshold:
            return "weak_gwas_signal"

    return None
