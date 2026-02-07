"""Co-occurrence network analysis."""

from __future__ import annotations

import pandas as pd
import scipy.stats as stats


def spearman_network(ko_table: pd.DataFrame, r_threshold: float, p_threshold: float) -> pd.DataFrame:
    kos = ko_table.index.tolist()
    edges = []
    for i, ko_a in enumerate(kos):
        for ko_b in kos[i + 1 :]:
            r, p = stats.spearmanr(ko_table.loc[ko_a], ko_table.loc[ko_b])
            if abs(r) >= r_threshold and p <= p_threshold:
                edges.append({"source": ko_a, "target": ko_b, "rho": r, "p_value": p})
    return pd.DataFrame(edges)
