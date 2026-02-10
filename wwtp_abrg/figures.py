"""Figure generation."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_top_kos(top_kos: pd.DataFrame, output_path: str | Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4))
    top_kos.T.plot(kind="bar", ax=ax)
    ax.set_title("Top KOs Over Time")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Abundance")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def _scatter_by_group(
    coords: pd.DataFrame,
    metadata: pd.DataFrame | None,
    group_col: str | None,
    annotate_samples: bool,
    output_path: str | Path,
    title: str,
    x_label: str,
    y_label: str,
) -> None:
    fig, ax = plt.subplots(figsize=(5, 4))
    if metadata is not None and group_col and group_col in metadata.columns:
        merged = coords.merge(metadata[["sample_id", group_col]], left_index=True, right_on="sample_id", how="left")
        for group, group_df in merged.groupby(group_col):
            ax.scatter(group_df.iloc[:, 0], group_df.iloc[:, 1], label=str(group))
    else:
        ax.scatter(coords.iloc[:, 0], coords.iloc[:, 1])
    if annotate_samples:
        for sample_id, row in coords.iterrows():
            ax.text(row.iloc[0], row.iloc[1], sample_id, fontsize=7)
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc="best", fontsize=8, frameon=False)
    fig.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def plot_pca(
    coords: pd.DataFrame,
    metadata: pd.DataFrame | None,
    group_col: str | None,
    annotate_samples: bool,
    output_path: str | Path,
) -> None:
    _scatter_by_group(coords, metadata, group_col, annotate_samples, output_path, "PCA", "PC1", "PC2")


def plot_pcoa(
    coords: pd.DataFrame,
    metadata: pd.DataFrame | None,
    group_col: str | None,
    annotate_samples: bool,
    output_path: str | Path,
) -> None:
    _scatter_by_group(coords, metadata, group_col, annotate_samples, output_path, "PCoA", "Axis 1", "Axis 2")


def plot_top30_heatmap(
    relative_abundance: pd.DataFrame,
    metadata: pd.DataFrame,
    top30: pd.DataFrame,
    output_path: str | Path,
    label_style: str = "KO_gene",
    log1p: bool = True,
) -> None:
    ko_ids = top30.index.tolist()
    heatmap_data = relative_abundance.loc[ko_ids]
    sample_order = heatmap_data.columns.tolist()
    if "day" in metadata.columns:
        sample_order = (
            metadata.set_index("sample_id")
            .loc[sample_order]
            .sort_values("day")
            .index.tolist()
        )
    heatmap_data = heatmap_data[sample_order]
    if log1p:
        heatmap_values = np.log1p(heatmap_data.values)
    else:
        heatmap_values = heatmap_data.values

    if label_style == "KO_gene" and "gene_name" in top30.columns:
        labels = [f"{ko} {top30.loc[ko, 'gene_name']}" for ko in ko_ids]
    else:
        labels = ko_ids

    fig, ax = plt.subplots(figsize=(10, max(4, len(ko_ids) * 0.35)))
    im = ax.imshow(heatmap_values, aspect="auto", interpolation="nearest", cmap="viridis")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xticks(range(len(sample_order)))
    ax.set_xticklabels(sample_order, rotation=45, ha="right", fontsize=7)
    ax.set_title("Top-30 Major ABRGs Heatmap")
    ax.set_xlabel("Sample")
    ax.set_ylabel("KO")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="log1p(Relative Abundance)" if log1p else "Relative Abundance")
    fig.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def plot_richness_over_time(richness_df: pd.DataFrame, metadata: pd.DataFrame, output_path: str | Path) -> None:
    merged = richness_df.merge(metadata[["sample_id", "day"]], on="sample_id", how="left")
    merged = merged.dropna(subset=["day"]).sort_values("day")
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(merged["day"], merged["ko_richness"], marker="o")
    ax.set_title("KO Richness Over Time")
    ax.set_xlabel("Day")
    ax.set_ylabel("KO Richness")
    fig.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def plot_top_kos_time_series(
    tidy_df: pd.DataFrame,
    output_path: str | Path,
    rolling_window: int = 7,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 4))
    for ko, group_df in tidy_df.groupby("KO"):
        group_df = group_df.sort_values("day")
        rolling = group_df["abundance"].rolling(rolling_window, min_periods=1).mean()
        ax.plot(group_df["day"], rolling, label=ko)
    ax.set_title("Top KOs Time Series (Rolling Mean)")
    ax.set_xlabel("Day")
    ax.set_ylabel("Relative Abundance")
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=7, frameon=False)
    fig.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
