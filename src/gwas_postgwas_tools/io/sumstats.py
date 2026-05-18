from __future__ import annotations

from pathlib import Path
from typing import Union

import pandas as pd

from gwas_postgwas_tools.config import InputSpec


def read_sumstats(spec: Union[InputSpec, str, Path]) -> pd.DataFrame:
    if isinstance(spec, (str, Path)):
        return _read_standard_sumstats(Path(spec))

    path = Path(spec.path)
    compression = spec.compression
    if compression is None and path.suffix == ".gz":
        compression = "gzip"

    df = pd.read_csv(path, sep=spec.sep, compression=compression)
    renamed = _rename_to_standard_columns(df, spec)
    _validate_loaded_columns(renamed, spec)
    return renamed


def _rename_to_standard_columns(df: pd.DataFrame, spec: InputSpec) -> pd.DataFrame:
    reverse_mapping = {source: target for target, source in spec.columns.items()}
    missing_source_columns = sorted(set(reverse_mapping).difference(df.columns))
    if missing_source_columns:
        raise ValueError(
            f"Input file {spec.path} is missing mapped source columns: {missing_source_columns}"
        )

    renamed = df.rename(columns=reverse_mapping).copy()
    return renamed.loc[:, ~renamed.columns.duplicated()]


def _validate_loaded_columns(df: pd.DataFrame, spec: InputSpec) -> None:
    required_columns = set(spec.columns)
    missing = sorted(required_columns.difference(df.columns))
    if missing:
        raise ValueError(f"Input file {spec.path} is missing standardized columns: {missing}")


def _read_standard_sumstats(path: Path) -> pd.DataFrame:
    compression = "gzip" if path.suffix == ".gz" else None
    df = pd.read_csv(path, sep="\t", compression=compression)
    return df
