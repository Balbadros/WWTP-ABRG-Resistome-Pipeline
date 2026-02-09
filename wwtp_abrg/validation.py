"""Validation checks for inputs."""

from __future__ import annotations

from typing import Iterable

import pandas as pd


def validate_samples(ko_table: pd.DataFrame, metadata: pd.DataFrame) -> None:
    """Ensure sample IDs match between KO table and metadata."""
    sample_cols = list(ko_table.columns)
    metadata_samples = metadata["sample_id"].tolist()
    missing_in_metadata = sorted(set(sample_cols) - set(metadata_samples))
    missing_in_table = sorted(set(metadata_samples) - set(sample_cols))
    if missing_in_metadata or missing_in_table:
        raise ValueError(
            "Sample ID mismatch."
            f" Missing in metadata: {missing_in_metadata}."
            f" Missing in KO table: {missing_in_table}."
        )


def validate_non_negative(df: pd.DataFrame, name: str) -> None:
    if (df < 0).any().any():
        raise ValueError(f"Negative values found in {name} table.")


def validate_no_missing(df: pd.DataFrame, name: str, required_cols: Iterable[str] | None = None) -> None:
    if df.isna().any().any():
        raise ValueError(f"Missing values found in {name} table.")
    if required_cols:
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns in {name}: {missing_cols}.")
