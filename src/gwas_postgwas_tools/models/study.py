from __future__ import annotations

from dataclasses import dataclass


@dataclass
class DatasetSpec:
    name: str
    path: str
    genome_build: str
    trait_type: str
