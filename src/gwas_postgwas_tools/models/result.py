from __future__ import annotations

from dataclasses import dataclass


@dataclass
class ColocResult:
    gene_id: str
    cell_type: str
    n_overlap_snps: int
    pp4: float
    passed_threshold: bool


@dataclass
class MRResult:
    gene_id: str
    cell_type: str
    method: str
    beta: float
    se: float
    p_value: float
