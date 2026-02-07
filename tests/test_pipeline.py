from pathlib import Path

from wwtp_abrg.pipeline import run_pipeline


def test_pipeline_runs(tmp_path: Path) -> None:
    config = {
        "project_name": "WWTP-ABRG-Resistome-Pipeline",
        "input": {
            "ko_table": "data/raw/ko_abundance.csv",
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
