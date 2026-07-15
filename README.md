# Metabarcoding of sedimentary ancient DNA suggests centuries of species turnover and introductions

## Overview

This repository contains the datasets and R scripts used to reproduce the analyses presented in the manuscript:

> *Metabarcoding of sedimentary ancient DNA suggests centuries of species turnover and introductions*

The study investigates long-term biodiversity dynamics preserved in marine sedimentary ancient DNA (sedaDNA) using metabarcoding of mitochondrial COI and nuclear 18S markers. Analyses include diversity estimation, community composition, hierarchical modelling, temporal trends, and reconstruction of native and non-indigenous species through time.

---

## Repository structure

```
.
├── Data/
│   ├── metadata.xlsx
│   ├── Supplementary Data S1_XRF_Dating.xlsx
│   ├── Supplementary Data S2.1_COI.xlsx
│   ├── Supplementary Data S2.2_18S.xlsx
│   ├── Supplementary Data S3 Relative Abundance Both Genes.xlsx
│   └── Supplementary Data S4_SpStatus.xlsx
│
├── Scripts/
│   ├── XRF_Data.R
│   ├── Map_Study_Area.R
│   ├── Richness_Barplots_MDS.R
│   ├── Barplots_All_Comunity.R
│   ├── NIS_Barplots_GAMs_Timelines.R
│   ├── Alpha_Diversity_Models.R
│   ├── Beta_Diversity_PERMANOVA_dbRDA.R
│   ├── HMSC_Models.R
│   └── Supplementary_Tables.R
│
├── Figures/
│
├── Supplementary/
│
└── README.md
```

---

## Data

### metadata.xlsx

Metadata for every sediment sample, including

- Sampling site
- Sediment core
- Sediment depth
- Sediment age
- Technical replicate information
- Additional sample descriptors

---

### Supplementary Data S1

Chronology and XRF geochemical data used to characterise sediment archives.

---

### Supplementary Data S2

Processed metabarcoding datasets

* COI ASV table
* 18S ASV table

including taxonomic assignments.

---

### Supplementary Data S3

Merged relative abundance dataset for both genetic markers.

---

### Supplementary Data S4

Species classification table including native, cryptogenic and non-indigenous species.

---

## Scripts

### Map_Study_Area.R

Creates Figure 1 showing the study area and sampling locations.

---

### XRF_Data.R

Plots sediment chronology and XRF elemental profiles.

---

### Richness_Barplots_MDS.R

Calculates alpha diversity metrics and generates:

- Richness
- Shannon diversity
- NMDS ordinations
- Diversity figures

---

### Barplots_All_Comunity.R

Produces stacked barplots showing relative taxonomic composition for:

- Kingdoms
- Phyla
- Major eukaryotic groups

for both COI and 18S datasets.

---

### NIS_Barplots_GAMs_Timelines.R

Analyses non-indigenous species through time, including:

- Relative abundance
- Species richness
- GAM models
- Detection timelines
- Heatmaps

---

### Alpha_Diversity_Models.R

Mixed-effects models evaluating technical and biological replication.

---

### Beta_Diversity_PERMANOVA_dbRDA.R

Community analyses including

- Bray-Curtis dissimilarity
- PERMANOVA
- dbRDA
- variance partitioning

---

### HMSC_Models.R

Hierarchical Modelling of Species Communities (HMSC)

including

- variance partitioning
- model diagnostics
- beta coefficients
- model fit statistics

---

## Software

Analyses were performed in

- R ≥ 4.3

Main packages include

- vegan
- lme4
- emmeans
- Hmsc
- mgcv
- ggplot2
- tidyverse
- openxlsx
- pairwiseAdonis

---

## Reproducibility

All analyses were performed using processed metabarcoding datasets after bioinformatic filtering and quality control.

Input files correspond to those provided as Supplementary Data accompanying the manuscript.

---

## Citation

If using these data or scripts, please cite: https://www.biorxiv.org/content/10.64898/2025.12.02.691833v1

*Paper citation will be added upon publication.*



CC-BY 4.0 for data and scripts.
