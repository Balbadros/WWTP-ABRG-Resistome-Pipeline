# Figure legend templates

## Figure PCA
**Template:** Principal component analysis (PCA) of KO relative abundances across SEM01–SEM14 samples. Points are colored by `{group_col}` when available (period/season/phase). Samples cluster by temporal grouping, indicating shifts in the resistome profile.

**Short interpretation:** Samples from later time points show separation along PC1, suggesting compositional turnover in dominant ABRG KOs.

## Figure PCoA
**Template:** Principal coordinates analysis (PCoA) using Bray–Curtis dissimilarity on KO relative abundance profiles. Points are colored by `{group_col}`. Separation indicates compositional differences between early and late (or seasonal) phases.

**Short interpretation:** The Bray–Curtis structure captures temporal divergence in KO composition, supporting distinct resistome phases.

## Figure Top KOs Over Time
**Template:** Top 15 KOs by mean abundance plotted across samples. Bars summarize relative abundance, highlighting dominant ABRG signals through time.

**Short interpretation:** A small subset of KOs consistently contributes the majority of ABRG signal across days.

## Heatmap: Top-30 Major ABRGs (Mixed)
**Template:** Heatmap of the top 30 major ABRG KOs (mixed mechanisms) using log1p relative abundance. KO IDs (and gene names, when available) are shown on the y-axis; samples are ordered by day.

**Short interpretation:** Dominant ABRG KOs show stable prevalence with episodic spikes across the 236-day series.

## Heatmap: Top-30 Major ABRGs (Efflux-only)
**Template:** Heatmap of the top 30 efflux-related ABRG KOs using log1p relative abundance. KO IDs and gene names are shown; samples ordered by day.

**Short interpretation:** Efflux systems are persistent across the time series, indicating a baseline multi-drug resistance signal.

## Richness Over Time
**Template:** KO richness (count of detected KOs) plotted by day. 

**Short interpretation:** Richness varies across days, suggesting temporal fluctuations in detected ABRG diversity.

## Top KOs Time Series
**Template:** Rolling mean (window = `{rolling_window_days}` days) of relative abundance for the top `{time_series_top_n}` KOs.

**Short interpretation:** Dominant KOs exhibit coherent temporal trends, with peaks aligned to specific sampling windows.
