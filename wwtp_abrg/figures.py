"""Figure generation."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
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


def plot_pca(coords: pd.DataFrame, output_path: str | Path) -> None:
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(coords["PC1"], coords["PC2"], color="#4C78A8")
    for sample_id, row in coords.iterrows():
        ax.text(row["PC1"], row["PC2"], sample_id, fontsize=8)
    ax.set_title("PCA")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    fig.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def plot_pcoa(coords: pd.DataFrame, output_path: str | Path) -> None:
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(coords.iloc[:, 0], coords.iloc[:, 1], color="#F58518")
    for sample_id, row in coords.iterrows():
        ax.text(row.iloc[0], row.iloc[1], sample_id, fontsize=8)
    ax.set_title("PCoA")
    ax.set_xlabel("Axis 1")
    ax.set_ylabel("Axis 2")
    fig.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)
