from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from wwtp_abrg.pipeline import run_pipeline


def test_pipeline_runs(tmp_path: Path) -> None:
    ko_source = Path("data/raw/ko_abundance.csv")
    ko_df = pd.read_csv(ko_source, index_col=0)
    ko_contiguous = np.ascontiguousarray(ko_df.values)
    ko_df = pd.DataFrame(ko_contiguous, index=ko_df.index, columns=ko_df.columns)
    ko_table_path = tmp_path / "ko_abundance.csv"
    ko_df.to_csv(ko_table_path)

    config = {
        "project_name": "WWTP-ABRG-Resistome-Pipeline",
        "input": {
            "ko_table": str(ko_table_path),
            "ko_stratified": "data/raw/ko_stratified.csv",
            "metadata": "metadata/sample_metadata.csv",
            "ko_annotations": "metadata/ko_annotations.csv",
        },
        "output": {
            "processed_dir": str(tmp_path / "processed"),
            "tables_dir": str(tmp_path / "tables"),
            "figures_dir": str(tmp_path / "figures"),
            "networks_dir": str(tmp_path / "networks"),
            "manifest": str(tmp_path / "run_manifest.json"),
        },
        "parameters": {
            "top_n": 3,
            "top30_n": 3,
            "prevalence_threshold": 0.5,
            "correlation_r": 0.1,
            "p_value": 1.0,
            "diff_method": "ttest",
            "random_seed": 42,
            "pcoa_metric": "braycurtis",
            "clustering_method": "average",
        },
        "scenarios": {
            "mixed": {"mechanism": []},
            "otu_sum": {"use_stratified": True},
        },
    }

    run_pipeline(config)

    expected = [
        Path(config["output"]["processed_dir"]) / "ko_raw_abundance.csv",
        Path(config["output"]["processed_dir"]) / "ko_relative_abundance.csv",
        Path(config["output"]["tables_dir"]) / "ko_richness.csv",
        Path(config["output"]["tables_dir"]) / "top_kos_over_time.csv",
        Path(config["output"]["tables_dir"]) / "differential_abundance.csv",
        Path(config["output"]["networks_dir"]) / "spearman_edges.csv",
        Path(config["output"]["manifest"]),
    ]

    for path in expected:
        assert path.exists()
