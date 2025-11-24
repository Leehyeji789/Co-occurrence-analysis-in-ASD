################################################################################
# pTDT (polygenic Transmission Disequilibrium Test) Analysis
# Calculate pTDT deviation and test significance with ANCOVA adjustment
################################################################################

library(tidyverse)
library(broom)

# Set working directory
setwd('/path/to/data/')

################################################################################
# 1. LOAD DATA
################################################################################

# Load phenotype association data (contains eur_phe_tb and eas_phe_tb)
load('path/to/output/Phenotype_association_data.rds')

# Load carrier definitions (contains combo_car_both, combo_car_mat, combo_car_pat, dnv_car)
load('path/to/output/Carrier_definitions.rds')

# Combine EUR and EAS samples (includes parents and children with genetic factors)
tb_all <- bind_rows(eur_phe_tb[[1]] %>% as.data.frame(), eas_phe_tb[[1]] %>% as.data.frame())

################################################################################
# 2. CALCULATE pTDT
################################################################################

# Extract parent polygenic scores
mom <- tb_all %>% filter(ROLE == 'parent' & imputed_sex == 'Female') %>% select(vcf_iid, vcf_fid, PS) %>% filter(!is.na(PS))
dad <- tb_all %>% filter(ROLE == 'parent' & imputed_sex == 'Male') %>% select(vcf_iid, vcf_fid, PS) %>% filter(!is.na(PS))

# Extract children and merge with parental PS
child <- tb_all %>% filter(ROLE == 'child') %>% filter(!is.na(PS)) %>%
  merge(., mom %>% select(vcf_fid, PS_maternal = PS), by = 'vcf_fid') %>%
  merge(., dad %>% select(vcf_fid, PS_paternal = PS), by = 'vcf_fid') %>%
  rowwise() %>%
  mutate(PS_parental = (PS_maternal + PS_paternal) / 2, pTDT = PS - PS_parental, maternal_pTDT = PS - PS_maternal, paternal_pTDT = PS - PS_paternal)

################################################################################
# 3. DEFINE GENETIC SUBGROUPS
################################################################################

child <- child %>% mutate(group = factor(case_when(isDNV == 1 ~ 'DNV', isDNV == 0 & isComb == 1 & (vcf_iid %in% combo_car_both) ~ 'Co-occurring rare variants\n(Maternal + Paternal)', isDNV == 0 & isComb == 1 & (vcf_iid %in% c(combo_car_mat, combo_car_pat)) ~ 'Co-occurring rare variants\n(Uni-parental)', TRUE ~ 'Others'), levels = c('Co-occurring rare variants\n(Maternal + Paternal)', 'Co-occurring rare variants\n(Uni-parental)', 'DNV', 'Others')))

# Print sample sizes per group
cat("Sample sizes per genetic subgroup:\n")
print(child %>% count(group, Group, imputed_sex) %>% arrange(imputed_sex) %>% pivot_wider(names_from = 'Group', values_from = 'n'))

################################################################################
# 4. ANCOVA-ADJUSTED pTDT TEST FOR AUTISM CASES
################################################################################

perform_adjusted_pTDT_test <- function(child_data, test_group, diagnosis_group) {
  set.seed(123)
  data <- child_data %>% filter(group == test_group & Group == diagnosis_group)
  n <- nrow(data)
  
  if (n < 10) {
    warning(paste("Insufficient sample size for group:", test_group, diagnosis_group, "(n =", n, ")"))
    return(NULL)
  }
  
  # Fit ANCOVA model: adjust pTDT for ancestry, cohort, and sex
  model_corrected <- lm(pTDT ~ ancestry + cohort + imputed_sex, data = data)
  
  # Extract residuals and add back the mean fitted value (covariate-adjusted response)
  adjusted_response <- residuals(model_corrected) + mean(model_corrected$fitted.values)
  
  # Test whether adjusted pTDT significantly differs from 0
  t_test <- t.test(adjusted_response, mu = 0, alternative = 'two.sided')
  
  tibble(group = test_group, diagnosis = diagnosis_group, pval = t_test$p.value, mean = mean(adjusted_response), sd = sd(adjusted_response), ci_lower = t_test$conf.int[1], ci_upper = t_test$conf.int[2], n = n) %>% mutate(se = sd / sqrt(n))
}

# Define groups to test
autism_groups <- c('DNV', 'Co-occurring rare variants\n(Maternal + Paternal)', 'Co-occurring rare variants\n(Uni-parental)', 'Others')
control_groups <- c('DNV', 'Co-occurring rare variants\n(Maternal + Paternal)', 'Others')

# Run tests for autism cases
stats_autism <- bind_rows(lapply(autism_groups, function(grp) perform_adjusted_pTDT_test(child, grp, 'Autism')))
stats_autism <- stats_autism %>% mutate(padj = p.adjust(pval, method = 'BH'), sig = ifelse(pval < 0.05, 'Yes', 'No'), sig_adj = ifelse(padj < 0.05, 'Yes', 'No'))

# Run tests for controls
stats_control <- bind_rows(lapply(control_groups, function(grp) perform_adjusted_pTDT_test(child, grp, 'Control')))
stats_control <- stats_control %>% mutate(padj = p.adjust(pval, method = 'BH'), sig = ifelse(pval < 0.05, 'Yes', 'No'), sig_adj = ifelse(padj < 0.05, 'Yes', 'No'))

################################################################################
# 5. COMBINE RESULTS
################################################################################

stats_combined <- bind_rows(stats_autism, stats_control)

cat("\n=== pTDT Results for Autism Cases ===\n")
print(stats_autism %>% select(group, mean, se, ci_lower, ci_upper, pval, padj, sig_adj, n))

cat("\n=== pTDT Results for Controls ===\n")
print(stats_control %>% select(group, mean, se, ci_lower, ci_upper, pval, padj, sig_adj, n))

################################################################################
# 6. SAVE RESULTS
################################################################################

save(child, stats_autism, stats_control, stats_combined, file = 'path/to/output/pTDT_analysis_results.rds')

# Export to Excel-friendly format
write.csv(stats_autism, 'path/to/output/pTDT_autism_cases.csv', row.names = FALSE)
write.csv(stats_control, 'path/to/output/pTDT_controls.csv', row.names = FALSE)
write.csv(stats_combined, 'path/to/output/pTDT_combined_results.csv', row.names = FALSE)