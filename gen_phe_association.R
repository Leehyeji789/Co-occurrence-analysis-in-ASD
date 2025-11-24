################################################################################
# Phenotype Association Analysis for Co-occurring Rare Variants
# Analysis script for Figure 4
################################################################################

## Load helper functions
source("path/to/function_definition.R") 

## Packages
library(tidyverse)
library(broom)
library(epitools)


############################################################
## Global settings
############################################################

project_dir <- "path/to/project_root"
setwd(project_dir)

################################################################################
# 1. DEFINE FUNCTIONS
################################################################################

# Function to prepare sample table with genetic factors and phenotypes
pheno_assoc_tb_case <- function(sample_tb, combo_car, dnv_car, clin, ps){
  all <- sample_tb %>%
    mutate(isDNV = ifelse(vcf_iid %in% dnv_car, 1, 0), isComb = ifelse(vcf_iid %in% combo_car, 1, 0)) %>%
    merge(., clin, by= 'ind_id', all.x=T) %>%
    merge(., ps, by = 'vcf_iid', all.x=T) %>%
    mutate(highPS = ifelse(PS >= 1.282, 1, 0), N_comb = case_when(isComb == 0 ~ 'Non-carrier', isComb ==1 ~ 'Carrier')) %>%
    filter(!is.na(isComb))
  
  case <- all %>% filter(Group == 'Autism')
  cntl <- all %>% filter(Group == 'Control')
  
  return(list(all, case, cntl))
}

# Function for genetic-phenotype association testing in cases
gene_phe_association <- function(phe_, gen_, case, linear, excludeDNV, adjustPS, adjustComb){
  # Build regression formula
  equation <- paste0(gen_, ' ~ ', phe_)
  if(excludeDNV == T) equation <- paste0(equation, ' + isDNV')
  if(adjustPS == T) equation <- paste0(equation, '+ PS')
  if(adjustComb == T) equation <- paste0(equation, '+ isComb')
  
  # Adjust for sex if multiple sexes present
  nSex <- case %>% filter(!is.na(eval(parse(text=phe_))) & !is.na(eval(parse(text=gen_)))) %>% count(imputed_sex) %>% filter(!is.na(imputed_sex)) %>% nrow()
  if(nSex > 1) equation <- paste0(equation, ' + imputed_sex')
  
  # Adjust for cohort if multiple cohorts present
  ncohort <- case %>% filter(!is.na(eval(parse(text=phe_))) & !is.na(eval(parse(text=gen_)))) %>% count(cohort) %>% filter(!is.na(cohort)) %>% nrow()
  if(ncohort > 1) equation <- paste0(equation, '+ cohort')
  
  # Adjust for ancestry if multiple ancestries present
  nAncestry <- case %>% filter(!is.na(eval(parse(text=phe_))) & !is.na(eval(parse(text=gen_)))) %>% count(ancestry) %>% filter(!is.na(ancestry)) %>% nrow()
  if(nAncestry > 1) equation <- paste0(equation, '+ ancestry')
  
  formula <- as.formula(equation)
  nY <- case %>% filter(!is.na(eval(parse(text=phe_))) & !is.na(eval(parse(text=gen_)))) %>% count(eval(parse(text=gen_))) %>% nrow()
  
  if(nY == 1){
    return(NULL)
  } else {
    if(linear==T){
      n <- case %>% filter(!is.na(eval(parse(text=phe_))) & !is.na(eval(parse(text=gen_)))) %>% nrow()
      linreg <- lm(formula, data = case)
      CI <- confint.default(linreg)[2,]
      tidy(linreg) %>% .[2,] %>% mutate(beta = estimate, l_ci = CI[1], u_ci = CI[2], n = n, phe = phe_, gen = gen_) %>% dplyr::select(-term)
    } else {
      n <- case %>% filter(!is.na(eval(parse(text=phe_))) & !is.na(eval(parse(text=gen_)))) %>% count(eval(parse(text=gen_)))
      logit <- glm(data = case, formula, family = 'binomial')
      CI <- exp(confint.default(logit))[2,]
      tidy(logit) %>% .[2,] %>% mutate(or = exp(estimate), l_ci = CI[1], u_ci = CI[2], n_1 = n[n$`eval(parse(text = gen_))`==1,]$n, n_0 = n[n$`eval(parse(text = gen_))`==0,]$n, n = n_1 + n_0, phe = phe_, gen = gen_) %>% dplyr::select(-term)
    }
  }
}

# Function for genetic-phenotype association testing in controls/parents
gene_phe_association_cntl <- function(phe_, gen_, case, linear, excludeDNV, adjustPS, adjustComb){
  # Build regression formula with similar logic as gene_phe_association
  equation <- paste0(gen_, ' ~ ', phe_)
  if(excludeDNV == T) equation <- paste0(equation, ' + isDNV')
  if(adjustPS == T) equation <- paste0(equation, '+ PS')
  if(adjustComb == T) equation <- paste0(equation, '+ isComb')
  
  nSex <- case %>% filter(!is.na(eval(parse(text=phe_))) & !is.na(eval(parse(text=gen_)))) %>% count(imputed_sex) %>% filter(!is.na(imputed_sex)) %>% nrow()
  if(nSex > 1) equation <- paste0(equation, ' + imputed_sex')
  
  ncohort <- case %>% filter(!is.na(eval(parse(text=phe_))) & !is.na(eval(parse(text=gen_)))) %>% count(cohort) %>% filter(!is.na(cohort)) %>% nrow()
  nPScohort <- case %>% filter(!is.na(eval(parse(text=phe_))) & !is.na(eval(parse(text=gen_)))) %>% count(cohort, eval(parse(text=gen_))) %>% filter(`eval(parse(text = gen_))` == 1) %>% nrow()
  if((ncohort > 1) & !(adjustPS == 1 & nPScohort == 1)) equation <- paste0(equation, '+ cohort')
  
  nAncestry <- case %>% filter(!is.na(eval(parse(text=phe_))) & !is.na(eval(parse(text=gen_)))) %>% count(ancestry) %>% filter(!is.na(ancestry)) %>% nrow()
  if(nAncestry > 1) equation <- paste0(equation, '+ ancestry')
  
  formula <- as.formula(equation)
  nY <- case %>% filter(!is.na(eval(parse(text=phe_))) & !is.na(eval(parse(text=gen_)))) %>% count(eval(parse(text=gen_))) %>% nrow()
  
  if(nY == 1){
    return(NULL)
  } else {
    if(linear==T){
      n <- case %>% filter(!is.na(eval(parse(text=phe_))) & !is.na(eval(parse(text=gen_)))) %>% nrow()
      linreg <- lm(formula, data = case)
      CI <- confint.default(linreg)[2,]
      tidy(linreg) %>% .[2,] %>% mutate(beta = estimate, l_ci = CI[1], u_ci = CI[2], n = n, phe = phe_, gen = gen_) %>% dplyr::select(-term)
    } else {
      n <- case %>% filter(!is.na(eval(parse(text=phe_))) & !is.na(eval(parse(text=gen_)))) %>% count(eval(parse(text=gen_)))
      logit <- glm(data = case, formula, family = 'binomial')
      CI <- exp(confint.default(logit))[2,]
      tidy(logit) %>% .[2,] %>% mutate(or = exp(estimate), l_ci = CI[1], u_ci = CI[2], n_1 = n[n$`eval(parse(text = gen_))`==1,]$n, n_0 = n[n$`eval(parse(text = gen_))`==0,]$n, n = n_1 + n_0, phe = phe_, gen = gen_) %>% dplyr::select(-term)
    }
  }
}

################################################################################
# 2. LOAD DATA
################################################################################

# Load polygenic scores
ps_eur <- rbind(read.delim('path/to/SSC_EUR.PS.txt') %>% select(vcf_iid, PS), read.delim('path/to/SPARK_EUR.PS.txt') %>% select(vcf_iid, PS))
ps_eas <- rbind(read.delim('path/to/SSC_EAS.PS.txt') %>% select(vcf_iid, PS), read.delim('path/to/SPARK_EAS.PS.txt') %>% select(vcf_iid, PS), read.delim('path/to/Korean_ASD.PS.txt') %>% select(vcf_iid, PS))

# Load co-occurring rare variant results
eur_candidates_file <- "path/to/eur_candidate_pairs.rds"
eas_candidates_file <- "path/to/eas_candidate_pairs.rds"
load(eur_candidates_file) 
load(eas_candidates_file) 

# Load clinical phenotype data
ssc_clin <- read.delim('path/to/SSC.pheno_table.txt') %>% dplyr::rename('ind_id'='sample_id') %>% merge(., read.csv('path/to/SSC_offspring.csv') %>% select(ind_id = SSC_ID, first_word, first_walk, first_phrase), by = 'ind_id')
sp_clin <- read.delim('path/to/SPARK.pheno_table.txt') %>% dplyr::rename('ind_id'='sample_id') %>% merge(., read.csv('path/to/SPARK.core.csv') %>% select(ind_id = SPARK_ID, first_word, first_walk), by = 'ind_id')
kor_clin <- read.delim('path/to/Korean.pheno_table.txt') %>% select(-Group) %>% merge(., read.delim('path/to/Korean_ASD_supertable.txt') %>% select(vcf_iid, Group), by = 'vcf_iid')
dev <- read.delim('path/to/Korean_ASD_clinical.txt') %>% select(vcf_iid, ADIR_walking, ADIR_word)
kor_clin <- merge(kor_clin, dev, by = 'vcf_iid')

clin <- bind_rows(ssc_clin %>% select(ind_id, age_m, ADOS_total, SCQ_current, SCQ_lifetime, SRS, ADOS_SA, ADIR_A, ADIR_Bv, ADOS_RRB, ADIR_C, RBSR, FSIQ, Non_verbal_IQ, VABS, VABS_communication, VABS_daily_living, VABS_socialization, VABS_motor_skill, DCDQ, ID, first_word, first_walk), sp_clin %>% select(ind_id, age_m, ADOS_total, SCQ_current, SCQ_lifetime, SRS, ADOS_SA, ADIR_A, ADIR_Bv, ADOS_RRB, ADIR_C, RBSR, FSIQ, Non_verbal_IQ, VABS, VABS_communication, VABS_daily_living, VABS_socialization, VABS_motor_skill, DCDQ, ID, first_word, first_walk), kor_clin %>% select(ind_id = vcf_iid, age_m, ADOS_total, SCQ_current, SCQ_lifetime, SRS, ADOS_SA, ADIR_A, ADIR_Bv, ADOS_RRB, ADIR_C, RBSR, FSIQ, Non_verbal_IQ, VABS, VABS_communication, VABS_daily_living, VABS_socialization, VABS_motor_skill, DCDQ, ID, first_word = ADIR_word, first_walk = ADIR_walking)) %>% unique()

# Load parent clinical phenotype data
ssc_par_clin <- read.csv('path/to/SSC_parents.csv') %>% mutate(vcf_fid = as.character(family)) %>% select(ind_id = SSC_ID, vcf_fid, BAPQ_PL = BAPQ_pragmatic, BAPQ_rigidity = BAPQ_rigid, BAPQ_aloof = BAPQ_aloof, BAPQ_overall, SRS = SRS_T)
kor_par_clin <- read.delim('path/to/Korean_ASD_parents_clinical.txt') %>% filter(ROLE == 'parent' & Group == 'Control') %>% select(1, 2, 5, 7, 11, 14, 16:16, 59:61, 66:75, ADOS_RRB = ADOS_RBB) %>% mutate(ind_id = vcf_iid)
par_clin <- bind_rows(ssc_par_clin, kor_par_clin)

# Define DNV carriers (using constraint genes)
constraint_genes <- cons %>% filter(oe_lof_upper < 0.2) %>% pull(transcript_id) %>% unique()
dnv_car_eur <- eur_input %>% filter(transmission_pattern == 'DNV') %>% filter((VariantType == 'protein_truncating'&transcript_id %in% constraint_genes) | (VariantType == 'missense'&mis_mean_score>=0.7&transcript_id %in% constraint_genes)) %>% pull(s) %>% unique()
dnv_car_eas <- eas_input %>% filter(transmission_pattern == 'DNV') %>% filter((VariantType == 'protein_truncating'&transcript_id %in% constraint_genes) | (VariantType == 'missense'&mis_mean_score>=0.7&transcript_id %in% constraint_genes)) %>% pull(s) %>% unique()
dnv_car <- c(dnv_car_eur, dnv_car_eas)

# Define carriers of co-occurring rare variants
combo_car <- c(eur_cand_inh %>% pull(s) %>% unique(), eas_cand_inh %>% pull(s) %>% unique(), eur_candidate %>% pull(Control_Samples) %>% paste(., collapse='|') %>% str_split(., '\\|') %>% unlist() %>% unique(), eas_candidate %>% pull(Control_Samples) %>% paste(., collapse='|') %>% str_split(., '\\|') %>% unlist() %>% unique())
combo_car_both <- c(eur_cand_inh %>% filter(Inherit_pattern %in% c('Maternal + Paternal')) %>% pull(s) %>% unique(), eas_cand_inh %>% filter(Inherit_pattern %in% c('Maternal + Paternal')) %>% pull(s) %>% unique(), eur_cand_cont_inh %>% filter(Inherit_pattern %in% c('Maternal + Paternal')) %>% pull(s) %>% unique(), eas_cand_cont_inh %>% filter(Inherit_pattern %in% c('Maternal + Paternal')) %>% pull(s) %>% unique())
combo_car_mat <- c(eur_cand_inh %>% filter((Inherit_pattern %in% c('Both Maternal'))) %>% pull(s) %>% unique(), eas_cand_inh %>% filter((Inherit_pattern %in% c('Both Maternal'))) %>% pull(s) %>% unique(), eur_cand_cont_inh %>% filter((Inherit_pattern %in% c('Both Maternal'))) %>% pull(s) %>% unique(), eas_cand_cont_inh %>% filter((Inherit_pattern %in% c('Both Maternal'))) %>% pull(s) %>% unique())
combo_car_pat <- c(eur_cand_inh %>% filter((Inherit_pattern %in% c('Both Paternal'))) %>% pull(s) %>% unique(), eas_cand_inh %>% filter((Inherit_pattern %in% c('Both Paternal'))) %>% pull(s) %>% unique(), eur_cand_cont_inh %>% filter((Inherit_pattern %in% c('Both Paternal'))) %>% pull(s) %>% unique(), eas_cand_cont_inh %>% filter((Inherit_pattern %in% c('Both Paternal'))) %>% pull(s) %>% unique())
combo_car_uni <- c(combo_car_pat, combo_car_mat)

'path/to/output/carrier_definitions.rds'

# Load sample information
all_samples_eur <- eur_input %>% pull(s) %>% unique()
eur <- read.delim('path/to/SFARI_sample_info.txt') %>% filter(overlap_wgs == ''|is.na(overlap_wgs)) %>% select(-overlap_wgs) %>% filter(vcf_iid %in% all_samples_eur)

all_samples_eas <- eas_input %>% pull(s) %>% unique()
sf_eas <- read.delim('path/to/SFARI_sample_info.txt') %>% filter(overlap_wgs == ''|is.na(overlap_wgs)) %>% filter(vcf_iid %in% all_samples_eas) %>% select(-overlap_wgs)
kor <- read.delim('path/to/Korean_ASD_supertable.txt') %>% filter(vcf_iid %in% all_samples_eas) %>% mutate(ancestry = 'EAS', familytype = ifelse(vcf_fid %in% c(read.delim('path/to/Korean_ASD_supertable.txt')%>% group_by(vcf_fid) %>% summarise(nASD = sum(Group == 'Autism')) %>% filter(nASD >= 2) %>% pull(vcf_fid)), 'Multiplex', 'Simplex'), sex = ifelse(sex == 'M', 'Male', 'Female')) %>% select(vcf_iid, vcf_fid, Group, ROLE, imputed_sex = sex, ancestry, familytype, cohort = data) %>% mutate(ind_id=vcf_iid)
eas <- bind_rows(sf_eas, kor)

################################################################################
# 3. RUN ANALYSES
################################################################################

# Define phenotype domains
domains <- c('ADOS_total', 'ADOS_SA', 'ADOS_RRB', 'SCQ_current', 'SCQ_lifetime', 'SRS', 'RBSR', 'FSIQ', 'Non_verbal_IQ', 'VABS', 'DCDQ', 'first_word', 'first_walk')

# === ANALYSIS 1: All cases (EUR + EAS combined) ===
eur_phe_tb <- pheno_assoc_tb_case(eur, combo_car, dnv_car, clin, ps_eur)
eas_phe_tb <- pheno_assoc_tb_case(eas, combo_car, dnv_car, clin, ps_eas)
tb_case <- bind_rows(eur_phe_tb[[2]], eas_phe_tb[[2]])

res_dnv <- lapply(domains, gene_phe_association, 'isDNV', tb_case, F, F, F, F) %>% bind_rows %>% dplyr::rename(pval = p.value)
res_ps <- lapply(domains, gene_phe_association, 'PS', tb_case, T, T, F, F) %>% bind_rows %>% dplyr::rename(pval = p.value)
res_comb <- lapply(domains, gene_phe_association, 'isComb', tb_case, F, T, F, F) %>% bind_rows %>% dplyr::rename(pval = p.value)

# Analysis by inheritance pattern
case <- pheno_assoc_tb_case(eur, combo_car_both, dnv_car, clin, ps_eur)[[2]]
case_eas <- pheno_assoc_tb_case(eas, combo_car_both, dnv_car, clin, ps_eas)[[2]]
tb_case_both <- bind_rows(case, case_eas)
res_comb_both <- lapply(domains, gene_phe_association, 'isComb', tb_case_both, F, T, F, F) %>% bind_rows %>% dplyr::rename(pval = p.value) %>% mutate(gen = 'Co-occurring rare variants\n(Maternal + Paternal)')

case <- pheno_assoc_tb_case(eur, combo_car_uni, dnv_car, clin, ps_eur)[[2]]
case_eas <- pheno_assoc_tb_case(eas, combo_car_uni, dnv_car, clin, ps_eas)[[2]]
tb_case_uni <- bind_rows(case, case_eas)
res_comb_uni <- lapply(domains, gene_phe_association, 'isComb', tb_case_uni, F, T, F, F) %>% bind_rows %>% dplyr::rename(pval = p.value) %>% mutate(gen = 'Co-occurring rare variants\n(Uni-parental)')

res_all <- bind_rows(res_dnv, res_ps, res_comb, res_comb_both, res_comb_uni)

# === ANALYSIS 2: Male cases ===
tb_case_male <- bind_rows(eur_phe_tb[[2]], eas_phe_tb[[2]]) %>% filter(imputed_sex == 'Male')
res_comb_male <- lapply(domains, gene_phe_association, 'isComb', tb_case_male, F, T, F, F) %>% bind_rows %>% dplyr::rename(pval = p.value)

case <- pheno_assoc_tb_case(eur, combo_car_both, dnv_car, clin, ps_eur)[[2]]
case_eas <- pheno_assoc_tb_case(eas, combo_car_both, dnv_car, clin, ps_eas)[[2]]
tb_case_both_male <- bind_rows(case, case_eas) %>% filter(imputed_sex == 'Male')
res_comb_both_male <- lapply(domains, gene_phe_association, 'isComb', tb_case_both_male, F, T, F, F) %>% bind_rows %>% dplyr::rename(pval = p.value) %>% mutate(gen = 'Co-occurring rare variants\n(Maternal + Paternal)')

case <- pheno_assoc_tb_case(eur, combo_car_uni, dnv_car, clin, ps_eur)[[2]]
case_eas <- pheno_assoc_tb_case(eas, combo_car_uni, dnv_car, clin, ps_eas)[[2]]
tb_case_uni_male <- bind_rows(case, case_eas) %>% filter(imputed_sex == 'Male')
res_comb_uni_male <- lapply(domains, gene_phe_association, 'isComb', tb_case_uni_male, F, T, F, F) %>% bind_rows %>% dplyr::rename(pval = p.value) %>% mutate(gen = 'Co-occurring rare variants\n(Uni-parental)')

res_dnv_male <- lapply(domains, gene_phe_association, 'isDNV', tb_case_both_male, F, F, F, F) %>% bind_rows %>% dplyr::rename(pval = p.value)
res_ps_male <- lapply(domains, gene_phe_association, 'PS', tb_case_both_male, T, T, F, F) %>% bind_rows %>% dplyr::rename(pval = p.value)

res_male <- bind_rows(res_dnv_male, res_ps_male, res_comb_male, res_comb_both_male, res_comb_uni_male)

# === ANALYSIS 3: Parents ===
par_eur <- eur_phe_tb[[1]] %>% as.data.frame() %>% select(-c(12:39)) %>% filter(ROLE == 'parent') %>% merge(., par_clin %>% select(1, 3:7, 10, 12:23), by.x = 'ind_id', by.y = 'ind_id', all.x=T)
par_eas <- eas_phe_tb[[1]] %>% as.data.frame() %>% select(-c(12:39)) %>% filter(ROLE == 'parent') %>% merge(., par_clin %>% select(1, 3:7, 10, 12:23), by.x = 'ind_id', by.y = 'ind_id', all.x=T)

combo_car_fid_eur <- eur_cand_inh %>% pull(vcf_fid) %>% unique()
dnv_car_fid_eur <- eur %>% filter(vcf_iid %in% dnv_car) %>% pull(vcf_fid) %>% unique()
par_eur <- par_eur %>% mutate(isChild_combo_carrier = ifelse(vcf_fid %in% combo_car_fid_eur, 1, 0), isChild_dnv_carrier = ifelse(vcf_fid %in% dnv_car_fid_eur, 1, 0))

combo_car_fid_eas <- eas_cand_inh %>% pull(vcf_fid) %>% unique()
dnv_car_fid_eas <- eas %>% filter(vcf_iid %in% dnv_car) %>% pull(vcf_fid) %>% unique()
par_eas <- par_eas %>% mutate(isChild_combo_carrier = ifelse(vcf_fid %in% combo_car_fid_eas, 1, 0), isChild_dnv_carrier = ifelse(vcf_fid %in% dnv_car_fid_eas, 1, 0))

tb_par <- bind_rows(par_eur, par_eas)
domains_par <- c('BAPQ_overall', 'BAPQ_aloof', 'BAPQ_PL', 'BAPQ_rigidity', 'SRS')

res_par_child_comb <- lapply(domains_par, gene_phe_association_cntl, 'isChild_combo_carrier', tb_par, F, F, T, F) %>% bind_rows %>% dplyr::rename(pval = p.value)
res_par_comb <- lapply(domains_par, gene_phe_association_cntl, 'isComb', tb_par, F, F, T, F) %>% bind_rows %>% dplyr::rename(pval = p.value)
res_par_dnv <- lapply(domains_par, gene_phe_association_cntl, 'isChild_dnv_carrier', tb_par, F, F, T, F) %>% bind_rows %>% dplyr::rename(pval = p.value)
res_par_ps <- lapply(domains_par, gene_phe_association_cntl, 'PS', tb_par, T, F, F, F) %>% bind_rows %>% dplyr::rename(pval = p.value)

res_par <- bind_rows(res_par_child_comb, res_par_comb, res_par_dnv, res_par_ps)

################################################################################
# 4. SAVE RESULTS
################################################################################

save(res_all, file = 'path/to/output/Phenotypic_association_results_AllCases.rds')
save(res_male, file = 'path/to/output/Phenotypic_association_results_MaleCases.rds')
save(res_par, file = 'path/to/output/Phenotypic_association_results_Parents.rds')
save(tb_case_both_male, file = 'path/to/output/Case_data_male_biparental.rds')

# Save intermediate data needed for pTDT analysis
save(eur_phe_tb, eas_phe_tb, file = 'path/to/output/Phenotype_association_data.rds')
save(combo_car, combo_car_both, combo_car_mat, combo_car_pat, combo_car_uni, dnv_car, 
     file = 'path/to/output/Carrier_definitions.rds')

cat("Analysis complete. Results saved.\n")