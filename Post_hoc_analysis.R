################################################################################
# Post-hoc Analysis: Testing Co-occurrence of Disrupted Gene Pairs in Cases
# Logistic regression with covariate adjustment
################################################################################

library(tidyverse)
library(data.table)
library(biomaRt)

## Set working directory
setwd('path/to/project_root')
source('path/to/function_definition.R')

################################################################################
# 1. LOAD DATA
################################################################################
## Load identified disrupted gene pairs
eur_candidates_file <- "path/to/eur_candidate_pairs.rds"
eas_candidates_file <- "path/to/eas_candidate_pairs.rds"
load(eur_candidates_file) 
load(eas_candidates_file) 

## Load input matrices (variant carrier status per sample)
eas_input <- fread('EAS_AllCohorts.PathogenicOnly.rarecomb_input_matrix.csv.gz')
eur_input <- fread("EUR_AllCohorts.PathogenicOnly.rarecomb_input_matrix.csv.gz")

## Load sample information and covariates
sm <- read.delim('path/to/somalier_ancestry.tsv') %>% dplyr::rename(sample_id = `X.sample_id`)
eur_info <- read.delim('path/to/EUR_sample_list.txt')
eas_info <- read.delim('path/to/EAS_sample_list.txt')
info <- bind_rows(eur_info, eas_info) %>% left_join(., sm, join_by(vcf_iid == sample_id))

## Load age data
sp <- fread('path/to/SPARK_WGS_mastertable.tsv') %>% select(vcf_iid = spid, age_m)
wes <- fread('path/to/SPARK_WES_mastertable.tsv') %>% select(vcf_iid = spid, age_m)
kor <- fread('path/to/Korean_supertable.txt') %>% select(vcf_iid, age_m)
ssc <- fread('path/to/SSC_age_individual.txt') %>% select(family, id = `id()`, age)
map <- fread('path/to/SSC_id_sample_map.csv')
ssc <- left_join(ssc, map, join_by(id == `SFARI ID`)) %>% select(vcf_iid = `Sample ID`, age_m = age) %>% filter(!is.na(vcf_iid))
age <- bind_rows(sp, wes, kor, ssc)
info <- left_join(info, age, join_by(vcf_iid)) %>% unique()

## Load sample QC metrics
sqc_sp_wgs <- fread('path/to/SPARK_WGS_sampleQC.txt')
sqc_sp_wes <- fread('path/to/SPARK_WES_sampleQC.txt')
sqc_ssc <- fread('path/to/SSC_sampleQC.txt')
sqc_kor_wgs <- fread('path/to/Korean_WGS_sampleQC.tsv.gz')
sqc_kor_wes <- fread('path/to/Korean_WES_sampleQC.tsv.gz')
sqc <- bind_rows(sqc_sp_wgs, sqc_kor_wes, sqc_ssc, sqc_kor_wgs, sqc_sp_wes) %>% select(-fam_id, -ROLE, -is_female)
info <- left_join(info, sqc, join_by(vcf_iid == s))

## Biomart for gene symbol annotation
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

################################################################################
# 2. DEFINE CAST FUNCTION FOR GENE-GENE INTERACTION TESTING
################################################################################

.get_p <- function(coefs, term) { if (!is.null(coefs) && term %in% rownames(coefs)) coefs[term, "Pr(>|z|)"] else NA_real_ }

CAST <- function(i, j, input, covar, covar_interest = c("sex","age_m", paste0("PC",1:10), "cohort","sequencing_type", "gq_mean_all", "dp_mean_all", "call_rate_all")) {
  
  ## Check if both genes exist in input
  if (!(i %in% colnames(input)) || !(j %in% colnames(input))) return(NULL)
  
  ## Convert to numeric carrier status (0/1)
  dom1 <- as.numeric(input[[i]] %in% c(1,"1",TRUE))
  dom2 <- as.numeric(input[[j]] %in% c(1,"1",TRUE))
  
  ## Count carriers
  dom12 <- dom1 * dom2
  m12 <- sum(dom12, na.rm=TRUE)
  m1  <- sum(dom1,  na.rm=TRUE); n1 <- sum(dom1 == 0, na.rm=TRUE)
  m2  <- sum(dom2,  na.rm=TRUE); n2 <- sum(dom2 == 0, na.rm=TRUE)
  
  ## Determine first/second gene (by carrier frequency)
  prop1 <- m1 / max(n1,1); prop2 <- m2 / max(n2,1)
  if (prop2 >= prop1) { first_vec <- dom2; second_vec <- dom1; first_name <- j; second_name <- i
  } else { first_vec <- dom1; second_vec <- dom2; first_name <- i; second_name <- j }
  
  ## Check variance in carrier status
  u_first  <- unique(first_vec[!is.na(first_vec)])
  u_second <- unique(second_vec[!is.na(second_vec)])
  base_out <- data.frame(genename1 = first_name, genename2 = second_name, N.double.carriers = m12, Ncarriers.gene1 = m1, Nnoncarriers.gene1 = n1, Ncarriers.gene2 = m2, Nnoncarriers.gene2 = n2, pvalueLRT = NA_real_, OR.dom = NA_real_, stringsAsFactors = FALSE)
  if (length(u_first) < 2 || length(u_second) < 2) return(base_out)
  
  ## Prepare covariates
  cov_use <- covar
  if ("cohort" %in% names(cov_use)) cov_use$cohort <- as.factor(cov_use$cohort)
  if ("sequencing_type" %in% names(cov_use)) cov_use$sequencing_type <- as.factor(cov_use$sequencing_type)
  if ("sex" %in% names(cov_use)) cov_use$sex <- as.numeric(cov_use$sex)
  if ("age_m" %in% names(cov_use)) cov_use$age_m <- as.numeric(cov_use$age_m)
  for (pc in intersect(paste0("PC",1:10), names(cov_use))) cov_use[[pc]] <- as.numeric(cov_use[[pc]])
  
  ## Build model dataframe
  model_df <- cbind(cov_use, first.gene = first_vec, second.gene = second_vec)
  model_df <- droplevels(model_df)
  
  ## Remove single-level factors and zero-variance numeric columns
  drop_cols <- logical(ncol(model_df))
  for (kk in seq_along(model_df)) {
    nm <- names(model_df)[kk]; x <- model_df[[kk]]
    if (nm %in% c("first.gene","second.gene")) next
    if (is.factor(x)) { if (nlevels(x) < 2) drop_cols[kk] <- TRUE
    } else if (is.numeric(x)) { ux <- unique(x); ux <- ux[!is.na(ux)]; if (length(ux) < 2) drop_cols[kk] <- TRUE }
  }
  if (any(drop_cols)) model_df <- model_df[, !drop_cols, drop=FALSE]
  
  ## Define explanatory variables
  expl_vars <- setdiff(names(model_df), c("first.gene","second.gene"))
  
  ## Fit logistic regression models (null and full)
  if (length(expl_vars) == 0) {
    mod_null <- glm(first.gene ~ 1, family=binomial, data=model_df)
    mod_full <- glm(first.gene ~ second.gene, family=binomial, data=model_df)
  } else {
    f_null <- as.formula(paste("first.gene ~", paste(expl_vars, collapse=" + ")))
    f_full <- as.formula(paste("first.gene ~", paste(c(expl_vars,"second.gene"), collapse=" + ")))
    mod_null <- glm(f_null, family=binomial, data=model_df)
    mod_full <- glm(f_full, family=binomial, data=model_df)
  }
  if (!isTRUE(mod_full$converged)) return(base_out)
  
  ## Likelihood ratio test for second.gene
  lrt <- mod_null$deviance - mod_full$deviance
  base_out$pvalueLRT <- pchisq(lrt, df=1, lower.tail=FALSE)
  cf <- coef(mod_full)
  if ("second.gene" %in% names(cf)) base_out$OR.dom <- unname(exp(cf["second.gene"]))
  
  ## Extract p-values for each covariate
  d1 <- try(drop1(mod_full, test="Chisq"), silent=TRUE)
  overall_p <- setNames(rep(NA_real_, length(covar_interest)), covar_interest)
  if (!inherits(d1,"try-error")) {
    rn <- rownames(d1)
    for (nm in covar_interest) { if (nm %in% rn) overall_p[nm] <- suppressWarnings(d1[nm, "Pr(>Chi)"]) }
  }
  
  ## Extract Wald p-values for numeric/binary variables
  coefs <- summary(mod_full)$coefficients
  for (nm in covar_interest) { if (is.na(overall_p[nm]) && nm %in% rownames(coefs)) overall_p[nm] <- .get_p(coefs, nm) }
  
  ## Append covariate p-values as columns
  pcols <- setNames(as.list(as.numeric(overall_p)), paste0("p_", names(overall_p)))
  base_out <- cbind(base_out, as.data.frame(pcols, check.names=FALSE))
  
  base_out
}

################################################################################
# 3. EUROPEAN ANCESTRY ANALYSIS
################################################################################

## Remove BOM and prepare data
names(eur_input) <- sub("^\ufeff", "", names(eur_input))
rawData_eur <- eur_input %>% filter(Output_group == 1)
input_eur <- column_to_rownames(rawData_eur, var = 'Sample_Name')
input_case_eur <- input_eur[, 1:(ncol(input_eur)-1)]

## Prepare covariates for EUR cases
info_df <- info %>% as.data.frame()
rownames(info_df) <- NULL
phenotype_eur <- column_to_rownames(info_df, var = 'vcf_iid')
covar_case_eur <- phenotype_eur[rownames(input_case_eur), ] %>% filter(Group == 'Autism') %>% select('sex','cohort', PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, age_m, dp_mean_all, gq_mean_all, call_rate_all) %>% mutate(sex = case_when(sex == 'Male' ~ 1, sex == 'Female' ~ 0), sequencing_type = ifelse(cohort == 'SSC', 'WGS', ifelse(grepl('WGS', cohort), 'WGS', 'WES'))) %>% separate(col = 'cohort', into = c('cohort', NA), sep = '_')

## Prepare gene pair list for EUR
pair_list_eur <- eur_candidate[, c(5,2)] %>% mutate(Item_1 = paste0('Input_', Item_1), Item_2 = paste0('Input_', Item_2))

## Run CAST for each EUR gene pair
results_eur <- lapply(1:nrow(pair_list_eur), function(idx) { CAST(pair_list_eur[idx, 1], pair_list_eur[idx, 2], input_case_eur, covar_case_eur) })
results_df_eur <- do.call(rbind, results_eur) %>% as.data.frame(stringsAsFactors = FALSE) %>% mutate(FDR = p.adjust(pvalueLRT, method = 'BH'), Item_1 = gsub('Input_', '', genename1), Item_2 = gsub('Input_', '', genename2))

## Annotate with gene symbols
bm1_eur <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), filters='ensembl_gene_id', values=results_df_eur$Item_1, mart=mart)
bm2_eur <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), filters='ensembl_gene_id', values=results_df_eur$Item_2, mart=mart)
results_df_eur <- left_join(results_df_eur, bm1_eur, join_by(Item_1==ensembl_gene_id)) %>% rename(gene_symbol_1=hgnc_symbol) %>% left_join(., bm2_eur, join_by(Item_2==ensembl_gene_id)) %>% rename(gene_symbol_2=hgnc_symbol) %>% select(Item_1, gene_symbol_1, Item_2, gene_symbol_2, 3:27)

## Save EUR results
write.table(results_df_eur, 'path/to/output/EUR_PostHoc_CoOccurrence_Results.txt', sep = '\t', row.names = F, col.names = T, quote=F)

################################################################################
# 4. EAST ASIAN ANCESTRY ANALYSIS
################################################################################

## Prepare data
rawData_eas <- eas_input %>% filter(Output_group == 1)
input_eas <- column_to_rownames(rawData_eas, var = 'Sample_Name')
input_case_eas <- input_eas[, 1:(ncol(input_eas)-1)]

## Prepare covariates for EAS cases
phenotype_eas <- column_to_rownames(info_df, var = 'vcf_iid')
covar_case_eas <- phenotype_eas[rownames(input_case_eas), ] %>% filter(Group == 'Autism') %>% select('sex','cohort', PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, age_m, dp_mean_all, gq_mean_all, call_rate_all) %>% mutate(sex = case_when(sex == 'Male' ~ 1, sex == 'Female' ~ 0), sequencing_type = ifelse(cohort == 'SSC', 'WGS', ifelse(grepl('WGS', cohort), 'WGS', 'WES'))) %>% separate(col = 'cohort', into = c('cohort', NA), sep = '_')

## Prepare gene pair list for EAS
pair_list_eas <- eas_candidate[, c(5,2)] %>% mutate(Item_1 = paste0('Input_', Item_1), Item_2 = paste0('Input_', Item_2))

## Run CAST for each EAS gene pair
results_eas <- lapply(1:nrow(pair_list_eas), function(idx) { CAST(pair_list_eas[idx, 1], pair_list_eas[idx, 2], input_case_eas, covar_case_eas) })
results_df_eas <- do.call(rbind, results_eas) %>% as.data.frame(stringsAsFactors = FALSE)
num_cols <- names(results_df_eas)[3:23]
results_df_eas[num_cols] <- lapply(results_df_eas[num_cols], as.numeric)
results_df_eas <- results_df_eas %>% mutate(FDR = p.adjust(pvalueLRT, method = 'BH'), Item_1 = gsub('Input_', '', genename1), Item_2 = gsub('Input_', '', genename2))

## Annotate with gene symbols
bm1_eas <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), filters='ensembl_gene_id', values=results_df_eas$Item_1, mart=mart)
bm2_eas <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), filters='ensembl_gene_id', values=results_df_eas$Item_2, mart=mart)
results_df_eas <- left_join(results_df_eas, bm1_eas, join_by(Item_1==ensembl_gene_id)) %>% rename(gene_symbol_1=hgnc_symbol) %>% left_join(., bm2_eas, join_by(Item_2==ensembl_gene_id)) %>% rename(gene_symbol_2=hgnc_symbol) %>% select(Item_1, gene_symbol_1, Item_2, gene_symbol_2, 3:27)

## Save EAS results
write.table(results_df_eas, 'path/to/output/EAS_PostHoc_CoOccurrence_Results.txt', sep = '\t', row.names = F, col.names = T, quote=F)