### TFM-TCGA-LUAD
# LUAD is heterogeneous; EGFR/KRAS mutations have limited prognostic value.
# This repository contains the R scripts used in my Master's Thesis for the development and validation of a four-gene transcriptomic risk signature in lung adenocarcinoma (LUAD).

# The analyses include:
# - Differential gene expression analysis (limma-voom)
# - Survival analysis (Kaplan–Meier and Cox regression)
# - Construction of the Risk Score
# - Time-dependent ROC analysis
# - Decision Curve Analysis
# - Functional analyses

### DATASET:
# The analyses were performed using the TCGA-LUAD cohort.
# Clinical, mutation and RNA-seq data can be downloaded from: cBioPortal: https://www.cbioportal.org/study/summary?id=luad_tcga_gdc

### Scripts:
# - ProyectoTFM.R	Main workflow including survival analyses, Risk Score construction and downstream analyses
# - DEARaw.R	Differential expression analysis using raw counts and limma-voom

### Author: Marina Delgado
# Master's Degree in Bioinformatics
# Universidad Internacional de Valencia (VIU)
