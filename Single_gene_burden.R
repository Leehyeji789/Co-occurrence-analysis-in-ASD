################################################################################
# Single Gene Burden Test Analysis
# Binomial test comparing variant carriers in cases vs controls
################################################################################

library(tidyverse)
library(biomaRt)
library(openxlsx)

## Set working directory
setwd('path/to/project_root')
source('path/to/function_definition.R')

## Setup Biomart for gene annotation
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

################################################################################
# DEFINE HELPER FUNCTION
################################################################################

## Function to perform binomial test and extract results
get_p_value <- function(x, n, p) { result <- binom.test(x = x, n = n, p = p, alternative = 'greater'); return(list(p_value = result$p.value, lower_ci = result$conf.int[1], upper_ci = result$conf.int[2], estimate = result$estimate)) }

################################################################################
# 1. EUROPEAN ANCESTRY ANALYSIS
################################################################################

## Load EUR variant data (DNV + rare inherited variants)
eur <- readRDS('path/to/EUR_AllCohorts.PathogenicOnly.long_format.rds')
eur_info <- read.delim('path/to/EUR_sample_list.txt')
eur_info <- eur_info %>% filter(vcf_iid %in% eur$s)
ncase_eur <- 14433
ncntl_eur <- 36752

## Create gene count table per sample
count_gene_eur <- merge(eur %>% select(-Group), eur_info, by.x = 's', by.y = 'vcf_iid') %>% group_by(gene_id) %>% count(s)
df_eur <- pivot_wider(count_gene_eur, names_from = 'gene_id', values_from = 'n') %>% merge(., eur_info, by.x = 's', by.y = 'vcf_iid', all.y = T) %>% replace(is.na(.), 0) %>% pivot_longer(cols = 2:(ncol(.)-5), names_to = 'gene_id', values_to = 'n')

## Perform per-gene binomial test
pergene_test_eur <- df_eur %>% group_by(Group, gene_id) %>% summarise(sum = sum(n, na.rm = TRUE), .groups = 'drop') %>% pivot_wider(names_from = 'Group', values_from = 'sum', values_fill = list(sum = 0)) %>% mutate(total_n = Autism + Control, p = ncase_eur / (ncntl_eur + ncase_eur), binom_results = mapply(get_p_value, x = Autism, n = total_n, p = p, SIMPLIFY = FALSE), RR = (Autism / Control) / (ncase_eur / (ncntl_eur + ncase_eur))) %>% unnest_wider(binom_results) %>% mutate(padj = p.adjust(p_value, method = 'BH'))

## Annotate with gene symbols
bm_eur <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values = pergene_test_eur$gene_id, mart = mart)
eur_single <- merge(bm_eur, pergene_test_eur, by.x = 'ensembl_gene_id', by.y = 'gene_id', all.y = T) %>% select(ensembl_gene_id, hgnc_symbol, nCase = Autism, nControl = Control, nTotal = total_n, expected_p = p, estimate, p_value, lower_ci, upper_ci, padj, RR)

cat("EUR single gene burden test complete.\n")
cat("Significant genes (FDR < 0.2):\n")
print(eur_single %>% filter(padj < 0.2) %>% select(hgnc_symbol, nCase, nControl, p_value, padj, RR))

################################################################################
# 2. EAST ASIAN ANCESTRY ANALYSIS
################################################################################

## Load EAS variant data (DNV + rare inherited variants)
eas <- readRDS('path/to/EAS_AllCohorts.PathogenicOnly.long_format.rds')
eas_info <- read.delim('path/to/EAS_sample_list.txt')
eas_info <- eas_info %>% filter(vcf_iid %in% eas$s)
ncase_eas <- 1057
ncntl_eas <- 5053

## Create gene count table per sample
count_gene_eas <- merge(eas %>% select(-Group), eas_info, by.x = 's', by.y = 'vcf_iid') %>% group_by(gene_id) %>% count(s)
df_eas <- pivot_wider(count_gene_eas, names_from = 'gene_id', values_from = 'n') %>% merge(., eas_info, by.x = 's', by.y = 'vcf_iid', all.y = T) %>% replace(is.na(.), 0) %>% pivot_longer(cols = 2:(ncol(.)-5), names_to = 'gene_id', values_to = 'n')

## Perform per-gene binomial test
pergene_test_eas <- df_eas %>% group_by(Group, gene_id) %>% summarise(sum = sum(n, na.rm = TRUE), .groups = 'drop') %>% pivot_wider(names_from = 'Group', values_from = 'sum', values_fill = list(sum = 0)) %>% mutate(total_n = Autism + Control, p = ncase_eas / (ncntl_eas + ncase_eas), binom_results = mapply(get_p_value, x = Autism, n = total_n, p = p, SIMPLIFY = FALSE), RR = (Autism / Control) / (ncase_eas / (ncntl_eas + ncase_eas))) %>% unnest_wider(binom_results) %>% mutate(padj = p.adjust(p_value, method = 'BH'))

## Annotate with gene symbols
bm_eas <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values = pergene_test_eas$gene_id, mart = mart)
eas_single <- merge(bm_eas, pergene_test_eas, by.x = 'ensembl_gene_id', by.y = 'gene_id', all.y = T) %>% select(ensembl_gene_id, hgnc_symbol, nCase = Autism, nControl = Control, nTotal = total_n, expected_p = p, estimate, p_value, lower_ci, upper_ci, padj, RR)

cat("EAS single gene burden test complete.\n")
cat("Significant genes (FDR < 0.2):\n")
print(eas_single %>% filter(padj < 0.2) %>% select(hgnc_symbol, nCase, nControl, p_value, padj, RR))

################################################################################
# 3. CREATE SUPPLEMENTARY TABLE
################################################################################

TableSupple <- createWorkbook()
addWorksheet(TableSupple, "EAS_single_gene_burden")
addWorksheet(TableSupple, "EUR_single_gene_burden")

## Write EAS results
eas_single %>% select(ensembl_gene_id, hgnc_symbol, nVar_case = nCase, nVar_control = nControl, nVar_total = nTotal, expected_ratio = expected_p, estimate, lower_ci, upper_ci, p_value, FDR = padj, RR) %>% writeData(wb = TableSupple, sheet = "EAS_single_gene_burden", x = ., startCol = 1, startRow = 1, rowNames = FALSE)

## Write EUR results
eur_single %>% select(ensembl_gene_id, hgnc_symbol, nVar_case = nCase, nVar_control = nControl, nVar_total = nTotal, expected_ratio = expected_p, estimate, lower_ci, upper_ci, p_value, FDR = padj, RR) %>% writeData(wb = TableSupple, sheet = "EUR_single_gene_burden", x = ., startCol = 1, startRow = 1, rowNames = FALSE)

## Save Excel file
saveWorkbook(wb = TableSupple, file = "path/to/output/Supplementary_SingleGeneBurden.xlsx", overwrite = TRUE)