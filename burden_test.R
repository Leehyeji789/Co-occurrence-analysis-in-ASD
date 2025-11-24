############################################################
## Burden tests for DNV and rare variants in EUR / EAS
##   1) DNV
##   2) DNV in constraint genes
##   3) Rare inherited variants
##   4) Rare inherited variants in constraint genes
##   5) Ultra-rare inherited variants
##   6) Ultra-rare inherited variants in constraint genes
############################################################

library(tidyverse)
library(data.table)
library(broom)

############################################################
## Global settings
############################################################

project_dir <- "path/to/project_root"
setwd(project_dir)

############################################################
## Constraint matrix
############################################################

constraint_file <- "resources/gnomad.v4.0.constraint_metrics.tsv"

cons <- read.delim(constraint_file) %>%
  dplyr::select(
    gene_symbol   = gene,
    transcript_id = transcript,
    pLI           = lof.pLI,
    oe_lof_upper  = lof.oe_ci.upper,
    mis_z         = mis.z_score
  )

## Example definition of constraint genes
cons_genes <- cons %>%
  dplyr::filter(oe_lof_upper < 0.2) %>%
  dplyr::pull(transcript_id)

############################################################
## Burden test function (adjusted counts + logistic regression)
############################################################

burden_test <- function(table) {
  table_child <- table %>%
    dplyr::filter(ROLE == "child") %>%
    dplyr::mutate(
      VariantType = ifelse(
        VariantType == "missense",
        "Damaging missense",
        "Protein truncating"
      ),
      VariantType = factor(
        VariantType,
        levels = c("Protein truncating", "Damaging missense")
      )
    )

  adj_model <- lm(n ~ cohort + sex, data = table_child)

  table_child <- table_child %>%
    dplyr::mutate(
      n_adj_raw = coef(adj_model)["(Intercept)"] + residuals(adj_model),
      mean_raw  = mean(n, na.rm = TRUE),
      mean_adj  = mean(n_adj_raw, na.rm = TRUE),
      shift     = mean_raw - mean_adj,
      n_adj     = n_adj_raw + shift
    )

  glm_mis <- glm(
    data = table_child %>%
      dplyr::mutate(Group = ifelse(Group == "Autism", 1, 0)) %>%
      dplyr::filter(VariantType == "Damaging missense"),
    Group ~ n_adj + cohort + sex,
    family = "binomial"
  )

  glm_ptv <- glm(
    data = table_child %>%
      dplyr::mutate(Group = ifelse(Group == "Autism", 1, 0)) %>%
      dplyr::filter(VariantType == "Protein truncating"),
    Group ~ n_adj + cohort + sex,
    family = "binomial"
  )

  logit_results <- function(glm_model) {
    ci <- exp(confint.default(glm_model))[2, ]
    broom::tidy(glm_model) %>%
      dplyr::mutate(
        odds_ratio = exp(estimate),
        l_ci       = ci[1],
        u_ci       = ci[2]
      ) %>%
      dplyr::select(
        term, estimate, std.error, statistic,
        p.value, odds_ratio, l_ci, u_ci
      ) %>%
      dplyr::slice(2)
  }

  mis_results <- logit_results(glm_mis)
  ptv_results <- logit_results(glm_ptv)

  burden_test_results <- dplyr::bind_rows(
    dplyr::mutate(mis_results, VariantType = "Damaging missense"),
    dplyr::mutate(ptv_results, VariantType = "Protein truncating")
  )

  group_stats <- table_child %>%
    dplyr::group_by(Group, VariantType) %>%
    dplyr::summarise(
      mean_value  = mean(n_adj, na.rm = TRUE),
      sd          = sd(n_adj, na.rm = TRUE),
      sample_size = dplyr::n_distinct(s),
      .groups     = "drop"
    )

  burden_test_results_combined <- burden_test_results %>%
    dplyr::left_join(group_stats, by = "VariantType") %>%
    dplyr::mutate(
      VariantType = factor(
        VariantType,
        levels = c("Protein truncating", "Damaging missense")
      )
    ) %>%
    dplyr::select(
      VariantType, Group,
      mean_value, sd, sample_size,
      estimate, std.error, statistic,
      odds_ratio, l_ci, u_ci,
      p.value
    )

  list(results = burden_test_results_combined)
}

############################################################
## Helper: build sample-level count table for burden_test
############################################################

make_burden_table <- function(variants,
                              sample_info,
                              adj_table = NULL) {
  # variants: s, VariantType, ...
  # sample_info: vcf_iid, Group, ROLE, sex, cohort
  # adj_table: optional (s, AdjustFactor_age)

  counts <- variants %>%
    dplyr::group_by(VariantType) %>%
    dplyr::count(s, name = "n")

  wide <- counts %>%
    tidyr::pivot_wider(
      names_from  = VariantType,
      values_from = n
    )

  merged <- wide %>%
    dplyr::right_join(
      sample_info %>% dplyr::filter(ROLE == "child"),
      by = c("s" = "vcf_iid")
    )

  merged[is.na(merged)] <- 0

  long <- merged %>%
    tidyr::pivot_longer(
      cols      = c("missense", "protein_truncating"),
      names_to  = "VariantType",
      values_to = "n"
    )

  if (!is.null(adj_table)) {
    long <- long %>%
      dplyr::left_join(adj_table, by = "s") %>%
      dplyr::mutate(
        AdjustFactor_age = dplyr::if_else(
          is.na(AdjustFactor_age) | AdjustFactor_age == 0,
          1,
          AdjustFactor_age
        ),
        n = n * AdjustFactor_age
      ) %>%
      dplyr::select(-AdjustFactor_age)
  }

  long
}

############################################################
## Data input (paths are anonymized)
############################################################

eur_info <- read.delim("path/to/EUR_sample_list.txt")
eas_info <- read.delim("path/to/EAS_sample_list.txt")

paternal_adj <- fread("path/to/paternal_age_adjustment.txt") %>%
  dplyr::select(s = vcf_iid, AdjustFactor_age)

paternal_adj$AdjustFactor_age[paternal_adj$AdjustFactor_age == 0] <- 1

eur_raw <- readRDS("path/to/EUR_AllCohorts.PathogenicOnly.long_format")
eas_raw <- readRDS("path/to/EAS_AllCohorts.PathogenicOnly.long_format")

############################################################
## Define variant subsets for each ancestry
## (user should fill in rare / ultra-rare filters)
############################################################

#### EUR ####

eur_dnv <- eur_raw %>%
  dplyr::filter(transmission_pattern == "DNV")

eur_dnv_cons <- eur_dnv %>%
  dplyr::filter(transcript_id %in% cons_genes)

eur_inherited <- eur_raw %>%
  dplyr::filter(transmission_pattern != "DNV")

## Rare inherited variants
eur_rare <- eur_inherited

eur_rare_cons <- eur_rare %>%
  dplyr::filter(transcript_id %in% cons_genes)

## Ultra-rare inherited variants
eur_ultra <- eur_inherited

eur_ultra_cons <- eur_ultra %>%
  dplyr::filter(transcript_id %in% cons_genes)

#### EAS ####

eas_dnv <- eas_raw %>%
  dplyr::filter(transmission_pattern == "DNV")

eas_dnv_cons <- eas_dnv %>%
  dplyr::filter(transcript_id %in% cons_genes)

eas_inherited <- eas_raw %>%
  dplyr::filter(transmission_pattern != "DNV" | is.na(transmission_pattern))

## Rare inherited variants
eas_rare <- eas_inherited

eas_rare_cons <- eas_rare %>%
  dplyr::filter(transcript_id %in% cons_genes)

## Ultra-rare inherited variants
eas_ultra <- eas_inherited

eas_ultra_cons <- eas_ultra %>%
  dplyr::filter(transcript_id %in% cons_genes)

############################################################
## Run burden tests – EUR
############################################################

eur_dnv_tab <- make_burden_table(
  variants    = eur_dnv,
  sample_info = eur_info,
  adj_table   = paternal_adj
)
eur_dnv_bt <- burden_test(eur_dnv_tab)

eur_dnv_cons_tab <- make_burden_table(
  variants    = eur_dnv_cons,
  sample_info = eur_info,
  adj_table   = paternal_adj
)
eur_dnv_cons_bt <- burden_test(eur_dnv_cons_tab)

eur_rare_tab <- make_burden_table(
  variants    = eur_rare,
  sample_info = eur_info
)
eur_rare_bt <- burden_test(eur_rare_tab)

eur_rare_cons_tab <- make_burden_table(
  variants    = eur_rare_cons,
  sample_info = eur_info
)
eur_rare_cons_bt <- burden_test(eur_rare_cons_tab)

eur_ultra_tab <- make_burden_table(
  variants    = eur_ultra,
  sample_info = eur_info
)
eur_ultra_bt <- burden_test(eur_ultra_tab)

eur_ultra_cons_tab <- make_burden_table(
  variants    = eur_ultra_cons,
  sample_info = eur_info
)
eur_ultra_cons_bt <- burden_test(eur_ultra_cons_tab)

############################################################
## Run burden tests – EAS
############################################################

eas_dnv_tab <- make_burden_table(
  variants    = eas_dnv,
  sample_info = eas_info,
  adj_table   = paternal_adj
)
eas_dnv_bt <- burden_test(eas_dnv_tab)

eas_dnv_cons_tab <- make_burden_table(
  variants    = eas_dnv_cons,
  sample_info = eas_info,
  adj_table   = paternal_adj
)
eas_dnv_cons_bt <- burden_test(eas_dnv_cons_tab)

eas_rare_tab <- make_burden_table(
  variants    = eas_rare,
  sample_info = eas_info
)
eas_rare_bt <- burden_test(eas_rare_tab)

eas_rare_cons_tab <- make_burden_table(
  variants    = eas_rare_cons,
  sample_info = eas_info
)
eas_rare_cons_bt <- burden_test(eas_rare_cons_tab)

eas_ultra_tab <- make_burden_table(
  variants    = eas_ultra,
  sample_info = eas_info
)
eas_ultra_bt <- burden_test(eas_ultra_tab)

eas_ultra_cons_tab <- make_burden_table(
  variants    = eas_ultra_cons,
  sample_info = eas_info
)
eas_ultra_cons_bt <- burden_test(eas_ultra_cons_tab)

############################################################
## Collect and export results
############################################################

add_meta <- function(res_tbl, ancestry, test_label) {
  res_tbl %>%
    dplyr::mutate(
      Ancestry = ancestry,
      Test     = test_label
    ) %>%
    dplyr::select(
      Ancestry, Test, VariantType, Group,
      mean_value, sd, sample_size,
      estimate, std.error, statistic,
      odds_ratio, l_ci, u_ci,
      p.value
    )
}

eur_results <- dplyr::bind_rows(
  add_meta(eur_dnv_bt$results,       "EUR", "DNV"),
  add_meta(eur_dnv_cons_bt$results,  "EUR", "DNV_in_constraint_genes"),
  add_meta(eur_rare_bt$results,      "EUR", "Rare_inherited"),
  add_meta(eur_rare_cons_bt$results, "EUR", "Rare_inherited_in_constraint_genes"),
  add_meta(eur_ultra_bt$results,     "EUR", "Ultra_rare_inherited"),
  add_meta(eur_ultra_cons_bt$results,"EUR", "Ultra_rare_inherited_in_constraint_genes")
)

eas_results <- dplyr::bind_rows(
  add_meta(eas_dnv_bt$results,       "EAS", "DNV"),
  add_meta(eas_dnv_cons_bt$results,  "EAS", "DNV_in_constraint_genes"),
  add_meta(eas_rare_bt$results,      "EAS", "Rare_inherited"),
  add_meta(eas_rare_cons_bt$results, "EAS", "Rare_inherited_in_constraint_genes"),
  add_meta(eas_ultra_bt$results,     "EAS", "Ultra_rare_inherited"),
  add_meta(eas_ultra_cons_bt$results,"EAS", "Ultra_rare_inherited_in_constraint_genes")
)

all_results <- dplyr::bind_rows(eur_results, eas_results) %>%
  dplyr::group_by(Ancestry) %>%
  dplyr::mutate(
    padj = p.adjust(p.value, method = "BH")
  ) %>%
  dplyr::ungroup()

dir.create("results", showWarnings = FALSE)

write.table(
  all_results,
  file      = "results/burden_results_DNV_and_rare_variants.txt",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)