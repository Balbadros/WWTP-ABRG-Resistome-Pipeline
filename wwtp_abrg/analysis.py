"""Analysis functions for ABRG pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Tuple

import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.cluster.hierarchy import linkage
from statsmodels.stats.multitest import multipletests
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import DistanceMatrix, permanova


@dataclass
class DiffResult:
    table: pd.DataFrame
    summary: pd.DataFrame


def compute_relative_abundance(ko_table: pd.DataFrame) -> pd.DataFrame:
    totals = ko_table.sum(axis=0)
    return ko_table.divide(totals, axis=1)


def compute_richness(ko_table: pd.DataFrame) -> pd.Series:
    return (ko_table > 0).sum(axis=0)


def top_kos_over_time(ko_table: pd.DataFrame, top_n: int) -> pd.DataFrame:
    mean_abundance = ko_table.mean(axis=1).sort_values(ascending=False)
    top_kos = mean_abundance.head(top_n).index
    return ko_table.loc[top_kos]


def differential_abundance(
    ko_table: pd.DataFrame,
    metadata: pd.DataFrame,
    group_col: str = "period",
    group_a: str = "early",
    group_b: str = "late",
    method: str = "ttest",
) -> DiffResult:
    group_a_samples = metadata.loc[metadata[group_col] == group_a, "sample_id"].tolist()
    group_b_samples = metadata.loc[metadata[group_col] == group_b, "sample_id"].tolist()
    results = []
    for ko, row in ko_table.iterrows():
        vals_a = row[group_a_samples].values
        vals_b = row[group_b_samples].values
        if method == "ttest":
            stat, p_value = stats.ttest_ind(vals_a, vals_b, equal_var=False)
        else:
            stat, p_value = stats.mannwhitneyu(vals_a, vals_b, alternative="two-sided")
        results.append((ko, stat, p_value, vals_a.mean(), vals_b.mean()))
    res_df = pd.DataFrame(
        results,
        columns=["KO", "stat", "p_value", "mean_early", "mean_late"],
    ).set_index("KO")
    res_df["q_value"] = multipletests(res_df["p_value"], method="fdr_bh")[1]
    summary = res_df.sort_values("q_value")
    return DiffResult(table=res_df, summary=summary)


def compute_pca(matrix: pd.DataFrame, n_components: int = 2) -> pd.DataFrame:
    data = np.ascontiguousarray(matrix.T.values)
    data = data - data.mean(axis=0)
    u, s, vh = np.linalg.svd(data, full_matrices=False)
    coords = u[:, :n_components] * s[:n_components]
    return pd.DataFrame(coords, index=matrix.columns, columns=[f"PC{i+1}" for i in range(n_components)])


def compute_pcoa(matrix: pd.DataFrame, metric: str = "braycurtis") -> Tuple[pd.DataFrame, pd.DataFrame]:
    data = np.ascontiguousarray(matrix.T.values)
    distance = beta_diversity(metric, data, ids=matrix.columns)
    ordination = pcoa(distance)
    coords = ordination.samples.iloc[:, :2]
    return coords, ordination.proportion_explained


def compute_bray_curtis(matrix: pd.DataFrame) -> pd.DataFrame:
    data = np.ascontiguousarray(matrix.T.values)
    distance = beta_diversity("braycurtis", data, ids=matrix.columns)
    return pd.DataFrame(distance.data, index=distance.ids, columns=distance.ids)


def compute_permanova(distance_df: pd.DataFrame, metadata: pd.DataFrame, group_col: str = "period") -> pd.DataFrame:
    distance = DistanceMatrix(distance_df.values, ids=distance_df.index)
    result = permanova(distance, metadata.set_index("sample_id"), column=group_col)
    return pd.DataFrame([result.to_dict()])


def compute_clustering(distance_df: pd.DataFrame, method: str = "average") -> pd.DataFrame:
    condensed = np.ascontiguousarray(distance_df.values)[
        np.triu_indices_from(distance_df.values, k=1)
    ]
    linkage_matrix = linkage(condensed, method=method)
    return pd.DataFrame(linkage_matrix, columns=["cluster1", "cluster2", "distance", "count"])


def compute_prevalence(ko_table: pd.DataFrame) -> pd.Series:
    return (ko_table > 0).sum(axis=1) / ko_table.shape[1]


def generate_top30(
    ko_table: pd.DataFrame,
    annotations: pd.DataFrame,
    top_n: int,
    prevalence_threshold: float,
    mechanism_filter: Iterable[str] | None = None,
) -> pd.DataFrame:
    mean_abundance = ko_table.mean(axis=1)
    prevalence = compute_prevalence(ko_table)
    merged = annotations.set_index("KO").copy()
    merged["mean_abundance"] = mean_abundance
    merged["prevalence"] = prevalence
    if mechanism_filter:
        merged = merged[merged["mechanism"].isin(mechanism_filter)]
    merged = merged[merged["prevalence"] >= prevalence_threshold]
    return merged.sort_values("mean_abundance", ascending=False).head(top_n)


def summarize_stratified(stratified: pd.DataFrame) -> pd.DataFrame:
    sample_cols = [col for col in stratified.columns if col not in {"KO", "Taxon"}]
    tidy = stratified.melt(id_vars=["KO", "Taxon"], value_vars=sample_cols, var_name="sample_id", value_name="abundance")
    summary = tidy.groupby(["KO", "Taxon"], as_index=False)["abundance"].sum()
    return summary


def stratified_to_ko_table(stratified: pd.DataFrame) -> pd.DataFrame:
    """Collapse stratified KO table to KO-by-sample abundance."""
    sample_cols = [col for col in stratified.columns if col not in {"KO", "Taxon"}]
    grouped = stratified.groupby("KO", as_index=False)[sample_cols].sum()
    return grouped.set_index("KO")[sample_cols]
