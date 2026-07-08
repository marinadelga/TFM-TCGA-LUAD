# TFM-TCGA-LUAD
## Transcriptomic Risk Signature in Lung Adenocarcinoma (LUAD)
This repository contains the R scripts used in my Master's Thesis:
**"Development and validation of a transcriptomic risk signature in lung adenocarcinoma (LUAD)."**

## Analyses included:
- Differential gene expression analysis (limma-voom)
- Survival analysis (Kaplan–Meier and Cox regression)
- Construction of the Risk Score
- Time-dependent ROC analysis
- Decision Curve Analysis
- Functional analyses

## Dataset:
The analyses were performed using the **TCGA-LUAD** cohort.
Clinical, mutation and RNA-seq data can be downloaded from:
- **cBioPortal**: https://www.cbioportal.org/study/summary?id=luad_tcga_gdc

## Scripts:
- **ProyectoTFM.R**	Main workflow including survival analyses, Risk Score construction and downstream analyses
- **DEARaw.R**	Differential expression analysis using raw counts and limma-voom

## Author
**Marina Delgado Valero**
Master's Degree in Bioinformatics
Universidad Internacional de Valencia (VIU)

## Requirements
The scripts were developed in **R (version 4.0.2)** in **RStudio** and require the installation of the corresponding R packages used throughout the analyses.
