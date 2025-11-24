# Scripts for Co-occurrence analysis in ASD
Analysis Scripts for Co-occurrence of rare variants implicates gene pairs in cytoskeletal pathways and is associated with increased severity for ASD

This repository contains analysis scripts for identifying and characterizing co-occurring rare variants in autism spectrum disorder (ASD) using the RareComb methodology across Korean and European ancestry cohorts.

## Overview

The analysis pipeline consists of variant quality control, de novo variant (DNV) calling, rare variant filtering, RareComb-based gene pair identification, statistical validation, and downstream functional characterization.

---

## Data Processing & Quality Control

### Variants_filtering.py
Whole-genome variant quality control and rare variant annotation using gnomAD.

### modified_denovo.py
Modified de novo variant caller for DRAGEN IGG pVCF files.

### vep-configuration.json
VEP configuration for functional annotation with additional plugins.

### run_PRScs.sh
Polygenic score calculation using PRS-CS with European LD reference panels.

---

## Helper Functions

### function_definition.R
Collection of helper functions for RareComb analysis and gene annotation.

---

## Gene Pair Identification

### ident_gene_pair.R
RareComb pipeline for identifying co-occurring pathogenic variant pairs across multiple cohorts.

---

## Statistical Testing & Validation

### burden_test.R
Case-control burden tests for DNV and rare inherited variants in constrained genes.

### cross_evaluation.R
Cross-ancestry replication analysis of co-occurring gene pairs.

### Post_hoc_analysis.R
Post-hoc logistic regression testing gene-gene co-occurrence with covariate adjustment.

### Simulation_test.R
Permutation-based significance testing for disrupted gene pairs.

### Single_gene_burden.R
Single-gene rare variant burden test using binomial test.

---

## Phenotypic & Genetic Architecture Analysis

### gen_phe_association.R
Association testing between co-occurring rare variants and clinical phenotypes.

### pTDT_calculation.R
Polygenic transmission disequilibrium test for carriers of co-occurring rare variants.

---

## Functional Characterization

### CoExp_correlation.R
Cell-type-specific co-expression analysis in developing human brain.

### FET_geneset.R
Gene set enrichment analysis for disrupted gene pairs.

---

## Figure Generation

### Figure1.R
Main Figure 1 generation for the manuscript.

### Figure2.R
Main Figure 2 generation for the manuscript.

### Figure3.R
Main Figure 3 generation for the manuscript.

### Figure4.R
Main Figure 4 generation for the manuscript.

---

## Dependencies

### Python
- hail 0.2
- numpy, pandas

### R
- Core: data.table, dplyr, tidyr, tibble
- Statistics: broom, pwr, openxlsx
- Visualization: ggplot2, ggpubr, cowplot, ggrepel
- Genomics: biomaRt
- Parallel: parallel, doParallel, foreach
- Special: anndata, reticulate (for Python integration)

### External Tools
- VEP (Ensembl Variant Effect Predictor) v110+
- PRS-CS (https://github.com/getian107/PRScs)
- pyRareComb (https://pypi.org/project/pyrarecomb)

---

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
