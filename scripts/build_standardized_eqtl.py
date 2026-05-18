from __future__ import annotations

import argparse
from pathlib import Path
import re
from statistics import NormalDist

import pandas as pd


CELL_TYPE_SPECS = {
    "astrocytes": {
        "dir_name": "astrocytes",
        "file_prefix": "Astrocytes",
        "cell_type": "astrocytes",
    },
    "microglia": {
        "dir_name": "microglia",
        "file_prefix": "Microglia",
        "cell_type": "microglia",
    },
    "endo": {
        "dir_name": "endo",
        "file_prefix": "Endothelial.cells",
        "cell_type": "endo",
    },
    "excitatory": {
        "dir_name": "excitatory",
        "file_prefix": "Excitatory.neurons",
        "cell_type": "excitatory",
    },
    "inh": {
        "dir_name": "inh",
        "file_prefix": "Inhibitory.neurons",
        "cell_type": "inh",
    },
    "oli": {
        "dir_name": "oli",
        "file_prefix": "Oligodendrocytes",
        "cell_type": "oli",
    },
    "opc": {
        "dir_name": "opc",
        "file_prefix": "OPCs",
        "cell_type": "opc",
    },
    "pb": {
        "dir_name": "pb",
        "file_prefix": "Periblasts",
        "cell_type": "pb",
    },
    "peri": {
        "dir_name": "peri",
        "file_prefix": "Pericytes",
        "cell_type": "peri",
    },
}

OUTPUT_COLUMNS = [
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
    "snp_id",
    "dist_tss",
    "SNP_id_hg38",
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a standardized eQTL table from per-cell-type chromosome split files."
    )
    parser.add_argument("--eqtl-root", required=True, help="Root directory containing cell-type folders.")
    parser.add_argument("--snp-map", required=True, help="Path to snp_pos.txt.")
    parser.add_argument("--output", required=True, help="Output TSV or TSV.GZ path.")
    parser.add_argument(
        "--cell-types",
        nargs="+",
        default=sorted(CELL_TYPE_SPECS),
        help="Cell types to include. Default includes all known cell types.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=1_000_000,
        help="Rows to read per chunk from each chromosome-split file.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()

    eqtl_root = Path(args.eqtl_root).expanduser().resolve()
    snp_map_path = Path(args.snp_map).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()

    snp_map_df = load_snp_map(snp_map_path)
    stats = stream_standardized_eqtl(
        eqtl_root=eqtl_root,
        snp_map_df=snp_map_df,
        output_path=output_path,
        selected_cell_types=args.cell_types,
        chunk_size=args.chunk_size,
    )

    print(f"rows\t{stats['rows']}")
    print(f"genes\t{stats['genes']}")
    print(f"cell_types\t{stats['cell_types']}")
    print(f"output\t{output_path}")
    return 0


def stream_standardized_eqtl(
    eqtl_root: Path,
    snp_map_df: pd.DataFrame,
    output_path: Path,
    selected_cell_types: list[str],
    chunk_size: int,
) -> dict[str, int]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        output_path.unlink()

    total_rows = 0
    unique_genes: set[str] = set()
    seen_cell_types: set[str] = set()
    header_written = False

    for cell_type_key in selected_cell_types:
        spec = CELL_TYPE_SPECS.get(cell_type_key)
        if spec is None:
            raise ValueError(f"Unsupported cell type '{cell_type_key}'.")

        files = list_eqtl_files(eqtl_root=eqtl_root, spec=spec)
        seen_cell_types.add(spec["cell_type"])

        for path in files:
            for chunk in pd.read_csv(
                path,
                sep=r"\s+",
                header=None,
                names=["gene_id", "snp_id", "dist_tss", "p_value", "beta"],
                chunksize=chunk_size,
            ):
                chunk["cell_type"] = spec["cell_type"]
                standardized = standardize_eqtl(eqtl_df=chunk, snp_map_df=snp_map_df)
                if standardized.empty:
                    continue

                standardized.to_csv(
                    output_path,
                    sep="\t",
                    index=False,
                    compression="gzip" if output_path.suffix == ".gz" else None,
                    mode="a",
                    header=not header_written,
                )
                header_written = True
                total_rows += len(standardized)
                unique_genes.update(standardized["gene_id"].astype(str).unique())

    if not header_written:
        pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(
            output_path,
            sep="\t",
            index=False,
            compression="gzip" if output_path.suffix == ".gz" else None,
        )

    return {
        "rows": total_rows,
        "genes": len(unique_genes),
        "cell_types": len(seen_cell_types),
    }


def list_eqtl_files(eqtl_root: Path, spec: dict[str, str]) -> list[Path]:
    cell_dir = eqtl_root / spec["dir_name"]
    if not cell_dir.exists():
        raise FileNotFoundError(f"Cell-type directory not found: {cell_dir}")

    files = sorted(
        path
        for path in cell_dir.iterdir()
        if path.is_file()
        and not path.name.startswith("._")
        and re.fullmatch(rf"{re.escape(spec['file_prefix'])}\.\d+", path.name)
    )
    if not files:
        raise FileNotFoundError(f"No chromosome-split eQTL files found in {cell_dir}")
    return files


def load_snp_map(path: Path) -> pd.DataFrame:
    snp_map = pd.read_csv(path, sep=r"\s+")
    required = {"SNP", "SNP_id_hg38", "effect_allele", "other_allele", "MAF"}
    missing = required.difference(snp_map.columns)
    if missing:
        raise ValueError(f"snp_pos.txt is missing required columns: {sorted(missing)}")

    trimmed = snp_map[["SNP", "SNP_id_hg38", "effect_allele", "other_allele", "MAF"]].copy()
    trimmed = trimmed.drop_duplicates(subset=["SNP"])
    return trimmed


def standardize_eqtl(eqtl_df: pd.DataFrame, snp_map_df: pd.DataFrame) -> pd.DataFrame:
    merged = eqtl_df.merge(
        snp_map_df,
        left_on="snp_id",
        right_on="SNP",
        how="inner",
    ).copy()
    if merged.empty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)

    coord = merged["SNP_id_hg38"].astype(str).str.replace("^chr", "", regex=True).str.split(":", expand=True)
    merged["chrom"] = coord[0]
    merged["position"] = pd.to_numeric(coord[1], errors="coerce")
    merged["effect_allele"] = merged["effect_allele"].astype(str).str.upper()
    merged["other_allele"] = merged["other_allele"].astype(str).str.upper()
    merged["maf"] = pd.to_numeric(merged["MAF"], errors="coerce")
    merged["p_value"] = pd.to_numeric(merged["p_value"], errors="coerce")
    merged["beta"] = pd.to_numeric(merged["beta"], errors="coerce")
    merged["dist_tss"] = pd.to_numeric(merged["dist_tss"], errors="coerce")
    merged["gene_id"] = merged["gene_id"].astype(str)

    merged = merged.dropna(subset=["position", "maf", "p_value", "beta"])
    if merged.empty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)

    merged["position"] = merged["position"].astype(int)
    merged["se"] = _compute_standard_error(merged["beta"], merged["p_value"])

    return merged[OUTPUT_COLUMNS].reset_index(drop=True)


def _compute_standard_error(beta: pd.Series, p_value: pd.Series) -> pd.Series:
    normal = NormalDist()
    clipped_p = pd.to_numeric(p_value, errors="coerce").fillna(1.0).clip(lower=1e-300, upper=1 - 1e-16)
    tail_probs = (1 - (clipped_p / 2)).clip(lower=1e-300, upper=1 - 1e-16)
    z_scores = tail_probs.map(lambda value: abs(normal.inv_cdf(float(value))))
    z_scores = pd.to_numeric(z_scores, errors="coerce").replace(0, pd.NA)
    se = pd.to_numeric(beta, errors="coerce").abs() / z_scores
    median_se = pd.to_numeric(se, errors="coerce").dropna().median()
    if pd.isna(median_se):
        median_se = 1.0
    return pd.to_numeric(se, errors="coerce").fillna(median_se).astype(float)


if __name__ == "__main__":
    raise SystemExit(main())
