"""Pipeline execution logic."""

from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict

from wwtp_abrg import __version__, analysis, figures, io, network, validation


def run_pipeline(config: Dict[str, Any]) -> None:
    logger = logging.getLogger("wwtp_abrg")
    logger.info("Loading inputs")

    ko_table = io.read_ko_table(config["input"]["ko_table"])
    metadata = io.read_metadata(config["input"]["metadata"])
    annotations = io.read_annotations(config["input"]["ko_annotations"])

    validation.validate_samples(ko_table, metadata)
    validation.validate_non_negative(ko_table, "KO")
    validation.validate_no_missing(ko_table, "KO")
    validation.validate_no_missing(metadata, "metadata", required_cols=["sample_id"])

    io.ensure_dirs(
        [
            config["output"]["processed_dir"],
            config["output"]["tables_dir"],
            config["output"]["figures_dir"],
            config["output"]["networks_dir"],
        ]
    )

    logger.info("Computing relative abundance")
    rel_abundance = analysis.compute_relative_abundance(ko_table)
    io.write_table(rel_abundance, Path(config["output"]["processed_dir"]) / "ko_relative_abundance.csv")
    io.write_table(ko_table, Path(config["output"]["processed_dir"]) / "ko_raw_abundance.csv")

    logger.info("Computing KO richness")
    richness = analysis.compute_richness(ko_table)
    richness_df = richness.reset_index()
    richness_df.columns = ["sample_id", "ko_richness"]
    io.write_tidy(richness_df, Path(config["output"]["tables_dir"]) / "ko_richness.csv")

    logger.info("Computing top KOs over time")
    top_kos = analysis.top_kos_over_time(ko_table, config["parameters"]["top_n"])
    io.write_table(top_kos, Path(config["output"]["tables_dir"]) / "top_kos_over_time.csv")

    logger.info("Differential abundance")
    comparisons = config.get("comparisons")
    if not comparisons:
        if "period" in metadata.columns:
            comparisons = [
                {
                    "name": "early_vs_late",
                    "group_col": "period",
                    "group_a": "early",
                    "group_b": "late",
                    "method": config["parameters"].get("diff_method", "ttest"),
                }
            ]
        else:
            comparisons = []
            logger.warning("No comparisons configured and 'period' column missing; skipping differential abundance.")
    for comparison in comparisons:
        group_col = comparison["group_col"]
        if group_col not in metadata.columns:
            logger.warning("Skipping comparison '%s' (missing column '%s').", comparison["name"], group_col)
            continue
        diff = analysis.differential_abundance(
            ko_table,
            metadata,
            group_col=group_col,
            group_a=comparison["group_a"],
            group_b=comparison["group_b"],
            method=comparison.get("method", "ttest"),
        )
        io.write_table(
            diff.table,
            Path(config["output"]["tables_dir"]) / f"differential_abundance_{comparison['name']}.csv",
        )
        io.write_table(
            diff.summary,
            Path(config["output"]["tables_dir"]) / f"differential_abundance_{comparison['name']}_summary.csv",
        )

    logger.info("PCA and PCoA")
    pca_coords = analysis.compute_pca(rel_abundance)
    io.write_table(pca_coords, Path(config["output"]["tables_dir"]) / "pca_coordinates.csv")

    pcoa_coords, pcoa_variance = analysis.compute_pcoa(rel_abundance, metric=config["parameters"]["pcoa_metric"])
    io.write_table(pcoa_coords, Path(config["output"]["tables_dir"]) / "pcoa_coordinates.csv")
    io.write_table(pcoa_variance.to_frame("proportion"), Path(config["output"]["tables_dir"]) / "pcoa_variance.csv")

    bray_curtis = analysis.compute_bray_curtis(rel_abundance)
    io.write_table(bray_curtis, Path(config["output"]["tables_dir"]) / "bray_curtis_distance.csv")

    permanova_df = analysis.compute_permanova(bray_curtis, metadata)
    io.write_tidy(permanova_df, Path(config["output"]["tables_dir"]) / "permanova.csv")

    clustering = analysis.compute_clustering(bray_curtis, method=config["parameters"]["clustering_method"])
    io.write_table(clustering, Path(config["output"]["tables_dir"]) / "clustering.csv")

    logger.info("Network analysis")
    edges = network.spearman_network(
        rel_abundance,
        r_threshold=config["parameters"]["correlation_r"],
        p_threshold=config["parameters"]["p_value"],
    )
    io.write_tidy(edges, Path(config["output"]["networks_dir"]) / "spearman_edges.csv")

    logger.info("Taxon attribution")
    stratified_path = config["input"].get("ko_stratified")
    stratified = None
    if stratified_path:
        stratified = io.read_stratified(stratified_path)
        stratified_summary = analysis.summarize_stratified(stratified)
        io.write_tidy(
            stratified_summary,
            Path(config["output"]["tables_dir"]) / "taxon_contributions.csv",
        )

    logger.info("Top 30 ABRGs scenarios")
    scenarios = config.get("scenarios", {})
    top30_results = {}
    for name, settings in scenarios.items():
        mechanism = settings.get("mechanism")
        prevalence = settings.get("min_prevalence", config["parameters"]["prevalence_threshold"])
        use_stratified = settings.get("use_stratified", False)
        scenario_table = ko_table
        if use_stratified:
            if stratified is None:
                raise ValueError("use_stratified is True but no stratified KO table was provided.")
            scenario_table = analysis.stratified_to_ko_table(stratified)
        top30 = analysis.generate_top30(
            scenario_table,
            annotations,
            top_n=config["parameters"]["top30_n"],
            prevalence_threshold=prevalence,
            mechanism_filter=mechanism if mechanism else None,
        )
        top30_results[name] = top30
        io.write_table(top30, Path(config["output"]["tables_dir"]) / f"top30_{name}.csv")

    logger.info("Generating Table 1")
    mixed_top30 = top30_results.get("mixed")
    if mixed_top30 is None:
        mixed_top30 = analysis.generate_top30(
            ko_table,
            annotations,
            top_n=config["parameters"]["top30_n"],
            prevalence_threshold=config["parameters"]["prevalence_threshold"],
            mechanism_filter=None,
        )
    table1 = analysis.build_table1(mixed_top30)
    io.write_tidy(table1, Path(config["output"]["tables_dir"]) / "Table1_Top30_major_ABRGs.csv")
    table1_md = table1.to_markdown(index=False)
    docs_path = Path("docs")
    docs_path.mkdir(parents=True, exist_ok=True)
    table1_md_path = docs_path / "Table1_Top30_major_ABRGs.md"
    with table1_md_path.open("w", encoding="utf-8") as md_handle:
        md_handle.write(table1_md + "\n")

    logger.info("Generating figures")
    figures.plot_top_kos(top_kos, Path(config["output"]["figures_dir"]) / "top_kos_over_time.png")
    group_col = None
    for candidate in ["period", "season", "phase"]:
        if candidate in metadata.columns:
            group_col = candidate
            break
    annotate_samples = config["parameters"].get("plot_annotate_samples", False)
    figures.plot_pca(
        pca_coords,
        metadata,
        group_col=group_col,
        annotate_samples=annotate_samples,
        output_path=Path(config["output"]["figures_dir"]) / "pca.png",
    )
    figures.plot_pcoa(
        pcoa_coords,
        metadata,
        group_col=group_col,
        annotate_samples=annotate_samples,
        output_path=Path(config["output"]["figures_dir"]) / "pcoa.png",
    )

    if top30_results.get("mixed") is None:
        top30_results["mixed"] = analysis.generate_top30(
            ko_table,
            annotations,
            top_n=config["parameters"]["top30_n"],
            prevalence_threshold=config["parameters"]["prevalence_threshold"],
            mechanism_filter=None,
        )
    figures.plot_top30_heatmap(
        rel_abundance,
        metadata,
        top30_results["mixed"],
        Path(config["output"]["figures_dir"]) / "heatmap_top30_mixed.png",
        label_style=config["parameters"].get("heatmap_label_style", "KO_gene"),
    )
    if top30_results.get("efflux_only") is None:
        top30_results["efflux_only"] = analysis.generate_top30(
            ko_table,
            annotations,
            top_n=config["parameters"]["top30_n"],
            prevalence_threshold=config["parameters"]["prevalence_threshold"],
            mechanism_filter=["Efflux"],
        )
    figures.plot_top30_heatmap(
        rel_abundance,
        metadata,
        top30_results["efflux_only"],
        Path(config["output"]["figures_dir"]) / "heatmap_top30_efflux_only.png",
        label_style=config["parameters"].get("heatmap_label_style", "KO_gene"),
    )

    if "day" in metadata.columns:
        tidy_ts = analysis.top_kos_time_series_tidy(
            rel_abundance,
            metadata,
            top_n=config["parameters"].get("time_series_top_n", 10),
        )
        io.write_tidy(
            tidy_ts,
            Path(config["output"]["tables_dir"]) / "top_kos_time_series_tidy.csv",
        )
        figures.plot_top_kos_time_series(
            tidy_ts,
            Path(config["output"]["figures_dir"]) / "top_kos_time_series.png",
            rolling_window=config["parameters"].get("rolling_window_days", 7),
        )
        figures.plot_richness_over_time(
            richness_df,
            metadata,
            Path(config["output"]["figures_dir"]) / "richness_over_time.png",
        )

    logger.info("Writing manifest")
    manifest = {
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "version": __version__,
        "parameters": config["parameters"],
        "inputs": config["input"],
        "outputs": config["output"],
    }
    Path(config["output"]["manifest"]).parent.mkdir(parents=True, exist_ok=True)
    with Path(config["output"]["manifest"]).open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)

    logger.info("Pipeline complete")
