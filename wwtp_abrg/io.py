"""I/O helpers for pipeline data."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd


def read_ko_table(path: str | Path) -> pd.DataFrame:
    """Read KO abundance table (rows KOs, columns samples)."""
    df = pd.read_csv(path, index_col=0)
    return df


def read_metadata(path: str | Path) -> pd.DataFrame:
    """Read sample metadata."""
    return pd.read_csv(path)


def read_annotations(path: str | Path) -> pd.DataFrame:
    """Read KO annotations (KO, gene_name, mechanism, antibiotic_class)."""
    return pd.read_csv(path)


def read_stratified(path: str | Path) -> pd.DataFrame:
    """Read stratified KO table (KO, Taxon, samples)."""
    return pd.read_csv(path)


def write_table(df: pd.DataFrame, path: str | Path) -> None:
    """Write DataFrame to CSV."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=True)


def write_tidy(df: pd.DataFrame, path: str | Path) -> None:
    """Write DataFrame to CSV without index."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def ensure_dirs(paths: Iterable[str | Path]) -> None:
    for path in paths:
        Path(path).mkdir(parents=True, exist_ok=True)
