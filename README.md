# WWTP-ABRG-Resistome-Pipeline

Reproducible bioinformatics pipeline to analyze antibiotic resistance gene (ABRG) KEGG Orthologs (KOs) inferred by PICRUSt2 from 16S data across 236 days (SEM01–SEM14). The pipeline regenerates KO abundance tables, richness, top KOs over time, differential abundance (early vs late), PCA/PCoA, Bray–Curtis/PERMANOVA, co-occurrence networks, and publication-ready tables/figures.

## Repository layout

```
data/
  raw/           # PICRUSt2 KO tables and stratified outputs
  processed/     # Derived KO tables
metadata/        # Sample metadata and KO annotations
results/
  tables/
  figures/
  networks/
docs/            # Paper assets and templates
scripts/
notebooks/
wwtp_abrg/        # Pipeline package
config/
```

## Inputs

Provide the following inputs via `config.yaml`:

- **PICRUSt2 KO table**: KO abundance table (rows = KO IDs, columns = sample IDs).
- **Optional stratified KO table**: PICRUSt2 stratified output for taxon attribution.
- **Sample metadata**: includes `sample_id` and `period` (early/late).
- **KO annotations**: KO ID, gene name, mechanism, antibiotic class.
- **Optional**: taxonomy or 16S cluster tables (expandable in future).

Example files are available under `data/raw/` and `metadata/` for CI.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python -m wwtp_abrg.run --config config/example.yaml
```

Or via Makefile:

```bash
make setup
make all
```

## Outputs

The pipeline generates:

- **KO abundance tables**: `data/processed/ko_raw_abundance.csv`, `data/processed/ko_relative_abundance.csv`
- **KO richness**: `results/tables/ko_richness.csv`
- **Top KOs over time**: `results/tables/top_kos_over_time.csv`
- **Differential abundance**: `results/tables/differential_abundance.csv` and summary
- **PCA/PCoA**: coordinate tables and plots
- **Bray–Curtis + PERMANOVA**: distance matrix and PERMANOVA table
- **Co-occurrence network edges**: `results/networks/spearman_edges.csv`
- **Taxon attribution**: `results/tables/taxon_contributions.csv`
- **Top 30 major ABRGs**: `results/tables/top30_*.csv`
- **Manifest**: `results/run_manifest.json` with timestamp + parameters

## Figure/Table mapping (paper assets)

| Output | Description | File |
| --- | --- | --- |
| Table 1 | Top 30 major ABRGs with KO ID, gene name, mechanism, antibiotic class | `results/tables/top30_mixed.csv` (or scenario-specific) |
| Fig. PCA | PCA scatter | `results/figures/pca.png` |
| Fig. PCoA | PCoA scatter | `results/figures/pcoa.png` |
| Fig. Top KOs | Top KOs over time | `results/figures/top_kos_over_time.png` |

## Configuration

Edit `config/example.yaml` to set:

- `top_n`, `top30_n`, `prevalence_threshold`
- `correlation_r`, `p_value`
- input and output paths
- scenario settings for `top30` tables (efflux-only, mixed, prevalence-filtered, OTU-sum)

## Extending differential abundance

Differential abundance uses t-tests (default) with FDR correction. The module is structured for future drop-in DESeq2 integration.

## Run manifest

Each run writes `results/run_manifest.json` with a timestamp, package version, inputs, outputs, and parameters for reproducibility.

## License

MIT.
