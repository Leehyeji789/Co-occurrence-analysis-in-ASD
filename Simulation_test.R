################################################################################
# Permutation Test Analysis for Disrupted Gene Pairs
# Step 1: Run 10,000 permutations per gene pair
#        - Randomly shuffle case/control labels 10,000 times
#        - Recalculate co-occurrence for each permutation
# Step 2: Calculate empirical p-values from permutation distribution
################################################################################

library(tidyverse)
library(parallel)
library(data.table)
library(biomaRt)
library(openxlsx)

## Set working directory
setwd('path/to/project_root')

## Setup Biomart for gene annotation
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

################################################################################
# 1. DEFINE PERMUTATION FUNCTION
################################################################################

run_permutation_analysis_fork_single <- function(df_path, comb, out_file, num_iterations, num_cores) {
  cat("\n========== [Permutation Analysis:", out_file, "] START ==========\n")
  
  ## Load input matrix (variant carrier status)
  if (grepl("\\.rds$|\\.RDS$", df_path, ignore.case = TRUE)) { cat("* Loading df from RDS:", df_path, "\n"); df <- readRDS(df_path)
  } else { cat("* Loading df with fread:", df_path, "\n"); df <- fread(df_path) }
  
  ## Prepare gene pair list
  setDT(comb)
  setnames(comb, old = c("Item_1","Item_2"), new = c("geneA","geneB"), skip_absent = TRUE)
  n_pairs <- nrow(comb)
  cat("* Number of gene pairs:", n_pairs, "\n")
  
  ## Define binomial test function
  run_binomial_test <- function(g1_vec, g2_vec) {
    n_all <- length(g1_vec); observed_count <- sum(g1_vec == 1 & g2_vec == 1)
    if (n_all == 0) return(c(NA, NA, NA))
    gene1_count <- sum(g1_vec == 1); gene2_count <- sum(g2_vec == 1)
    expected_count <- (gene1_count / n_all) * (gene2_count / n_all) * n_all
    expected_prop <- ifelse(n_all > 0, expected_count / n_all, NA_real_)
    test_res <- binom.test(x = observed_count, n = n_all, p = expected_prop, alternative = "greater")
    c(test_res$p.value, expected_prop, observed_count / n_all)
  }
  
  ## Run permutations in parallel using fork (mclapply)
  cat(sprintf("* Running %d permutations with %d cores...\n", num_iterations, num_cores))
  all_out <- mclapply(X = seq_len(num_iterations), mc.cores = num_cores, FUN = function(iter_i) {
    set.seed(iter_i)
    perm_group <- sample(df$Output_group)
    idx_case <- which(perm_group == 1); idx_cont <- which(perm_group == 0)
    out_mat <- matrix(NA_real_, nrow = n_pairs, ncol = 6)
    for (r in seq_len(n_pairs)) {
      g1 <- paste0("Input_", comb$geneA[r]); g2 <- paste0("Input_", comb$geneB[r])
      if (length(idx_case) > 0) { tmp_c <- run_binomial_test(df[[g1]][idx_case], df[[g2]][idx_case]); out_mat[r, 1] <- tmp_c[1]; out_mat[r, 2] <- tmp_c[2]; out_mat[r, 3] <- tmp_c[3] }
      if (length(idx_cont) > 0) { tmp_d <- run_binomial_test(df[[g1]][idx_cont], df[[g2]][idx_cont]); out_mat[r, 4] <- tmp_d[1]; out_mat[r, 5] <- tmp_d[2]; out_mat[r, 6] <- tmp_d[3] }
    }
    out_mat
  })
  
  ## Combine results from all iterations
  cat("* Binding all iteration results...\n")
  big_mat <- do.call(rbind, all_out)
  iter_vec <- rep(seq_len(num_iterations), each = n_pairs); geneA_vec <- rep(comb$geneA, times = num_iterations); geneB_vec <- rep(comb$geneB, times = num_iterations)
  final_dt <- data.table(iteration = iter_vec, Item_1 = geneA_vec, Item_2 = geneB_vec, Case_p_value = big_mat[,1], Case_exp_prop = big_mat[,2], Case_obs_prop = big_mat[,3], Cont_p_value = big_mat[,4], Cont_exp_prop = big_mat[,5], Cont_obs_prop = big_mat[,6])
  
  ## Save results
  cat("* Saving final result to:", out_file, "\n")
  saveRDS(final_dt, file = out_file)
  cat("\n========== [", out_file, "] DONE ==========\n")
  invisible(final_dt)
}

################################################################################
# 2. LOAD CANDIDATE GENE PAIRS
################################################################################

eur_candidates_file <- "path/to/eur_candidate_pairs.rds"
eas_candidates_file <- "path/to/eas_candidate_pairs.rds"
load(eur_candidates_file) 
load(eas_candidates_file)

## Extract gene pair combinations
eas_comb <- eas_candidate %>% select(Item_1, Item_2)
eur_comb <- eur_candidate %>% select(Item_1, Item_2)

################################################################################
# 3. RUN PERMUTATION TESTS
################################################################################

## EAS permutation (10,000 iterations)
run_permutation_analysis_fork_single(df_path = "path/to/EAS_AllCohorts.PathogenicOnly.rarecomb_input_matrix.csv.gz", comb = eas_comb, out_file = "path/to/output/EAS_permutation_results.rds", num_iterations = 10000, num_cores = 30)

## EUR permutation (10,000 iterations)
run_permutation_analysis_fork_single(df_path = "path/to/EUR_AllCohorts.PathogenicOnly.rarecomb_input_matrix.csv.gz", comb = eur_comb, out_file = "path/to/output/EUR_permutation_results.rds", num_iterations = 10000, num_cores = 30)

################################################################################
# 4. ANALYZE PERMUTATION RESULTS - EAS
################################################################################

## Load permutation results
eas_perm <- readRDS('path/to/output/EAS_permutation_results.rds')

## Add gene pair ID and merge with observed results
eas_perm <- eas_perm %>% mutate(id = paste0(Item_1, '/', Item_2)) %>% left_join(., eas_candidate %>% select(Item_1, Item_2, Case_pvalue_more, Cont_pvalue_more, Case_Exp_Prob_Combo, Case_Obs_Prob_Combo, Cont_Obs_Prob_Combo, Cont_Exp_Prob_Combo))

## Calculate delta (observed - expected proportion) for each iteration
eas_perm <- eas_perm %>% group_by(id) %>% mutate(simulated_delta = Case_obs_prop - Case_exp_prop, original_delta = Case_Obs_Prob_Combo - Case_Exp_Prob_Combo) %>% ungroup()

## Summarize permutation distribution per gene pair
eas_summary <- eas_perm %>% group_by(id) %>% summarise(Item_1 = first(Item_1), Item_2 = first(Item_2), original_delta = first(original_delta), mean_sim = mean(simulated_delta), sd_sim = sd(simulated_delta), z_score = (first(original_delta) - mean(simulated_delta)) / sd(simulated_delta), empirical_P = 2 * (1 - pnorm(abs(z_score))), quantile_original = mean(simulated_delta >= first(original_delta))) %>% ungroup() %>% unique()

## Annotate with gene symbols
bm1_eas <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values = eas_summary$Item_1, mart = mart) %>% unique()
bm2_eas <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values = eas_summary$Item_2, mart = mart) %>% unique()
eas_table <- left_join(eas_summary, bm1_eas, join_by(Item_1 == ensembl_gene_id)) %>% rename(gene_symbol_1 = hgnc_symbol) %>% left_join(., bm2_eas, join_by(Item_2 == ensembl_gene_id)) %>% rename(gene_symbol_2 = hgnc_symbol) %>% select(Item_1, gene_symbol_1, Item_2, gene_symbol_2, original_delta, mean_simulated_delta = mean_sim, sd_simulated_delta = sd_sim, empirical_P) %>% unique()

cat("EAS permutation analysis complete: ", nrow(eas_table), " gene pairs\n")

################################################################################
# 5. ANALYZE PERMUTATION RESULTS - EUR
################################################################################

## Load permutation results
eur_perm <- readRDS('path/to/output/EUR_permutation_results.rds')

## Add gene pair ID and merge with observed results
eur_perm <- eur_perm %>% mutate(id = paste0(Item_1, '/', Item_2)) %>% left_join(., eur_candidate %>% select(Item_1, Item_2, Case_pvalue_more, Cont_pvalue_more, Case_Exp_Prob_Combo, Case_Obs_Prob_Combo, Cont_Obs_Prob_Combo, Cont_Exp_Prob_Combo))

## Calculate delta (observed - expected proportion) for each iteration
eur_perm <- eur_perm %>% group_by(id) %>% mutate(simulated_delta = Case_obs_prop - Case_exp_prop, original_delta = Case_Obs_Prob_Combo - Case_Exp_Prob_Combo) %>% ungroup()

## Summarize permutation distribution per gene pair
eur_summary <- eur_perm %>% group_by(id) %>% summarise(Item_1 = first(Item_1), Item_2 = first(Item_2), original_delta = first(original_delta), mean_sim = mean(simulated_delta), sd_sim = sd(simulated_delta), z_score = (first(original_delta) - mean(simulated_delta)) / sd(simulated_delta), empirical_P = 2 * (1 - pnorm(abs(z_score))), quantile_original = mean(simulated_delta >= first(original_delta))) %>% ungroup() %>% unique()

## Annotate with gene symbols
bm1_eur <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values = eur_summary$Item_1, mart = mart) %>% unique()
bm2_eur <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values = eur_summary$Item_2, mart = mart) %>% unique()
eur_table <- left_join(eur_summary, bm1_eur, join_by(Item_1 == ensembl_gene_id)) %>% rename(gene_symbol_1 = hgnc_symbol) %>% left_join(., bm2_eur, join_by(Item_2 == ensembl_gene_id)) %>% rename(gene_symbol_2 = hgnc_symbol) %>% select(Item_1, gene_symbol_1, Item_2, gene_symbol_2, original_delta, mean_simulated_delta = mean_sim, sd_simulated_delta = sd_sim, empirical_P) %>% unique()

cat("EUR permutation analysis complete: ", nrow(eur_table), " gene pairs\n")

################################################################################
# 6. CREATE SUPPLEMENTARY TABLE (EXCEL)
################################################################################

TableSupple <- createWorkbook()
addWorksheet(TableSupple, "EAS_pairs")
addWorksheet(TableSupple, "EUR_pairs")
writeData(wb = TableSupple, sheet = "EAS_pairs", x = eas_table, startCol = 1, startRow = 1, rowNames = FALSE)
writeData(wb = TableSupple, sheet = "EUR_pairs", x = eur_table, startCol = 1, startRow = 1, rowNames = FALSE)
saveWorkbook(wb = TableSupple, file = "path/to/output/Supplementary_PermutationTest.xlsx", overwrite = TRUE)