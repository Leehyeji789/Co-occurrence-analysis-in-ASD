############################################################
## Cross-ancestry replication analysis
############################################################

## Packages
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(broom)
library(pwr)
library(openxlsx)

source("path/to/function_definition.R")


############################################################
## Helper functions
############################################################

build_cross_replication_summary <- function(candidate_dt, gene_by_sample_dt) {
  # candidate_dt       : data.table with gene_symbol_1, gene_symbol_2, id
  # gene_by_sample_dt  : data.table with s, gene_symbol, Group
  
  n_total_dt <- gene_by_sample_dt[, .(n_total = uniqueN(s)), by = Group]
  
  pairs_dt <- unique(candidate_dt[, .(id, g1 = gene_symbol_1, g2 = gene_symbol_2)])
  
  g1_car <- merge(
    pairs_dt[, .(id, gene_symbol = g1)],
    gene_by_sample_dt,
    by = "gene_symbol",
    allow.cartesian = TRUE
  )[, .(id, s, Group)]
  g1_car <- unique(g1_car)
  
  g2_car <- merge(
    pairs_dt[, .(id, gene_symbol = g2)],
    gene_by_sample_dt,
    by = "gene_symbol",
    allow.cartesian = TRUE
  )[, .(id, s, Group)]
  g2_car <- unique(g2_car)
  
  n_g1_dt   <- g1_car[, .(n_g1   = uniqueN(s)), by = .(id, Group)]
  n_g2_dt   <- g2_car[, .(n_g2   = uniqueN(s)), by = .(id, Group)]
  n_both_dt <- merge(g1_car, g2_car, by = c("id", "s", "Group"))[
    , .(n_both = uniqueN(s)), by = .(id, Group)
  ]
  
  base_dt <- CJ(
    id    = unique(pairs_dt$id),
    Group = unique(n_total_dt$Group)
  )
  
  cross_replication_summary_dt <-
    merge(base_dt,  n_total_dt, by = "Group",          all.x = TRUE, sort = FALSE) |>
    merge(          n_g1_dt,    by = c("id", "Group"), all.x = TRUE, sort = FALSE) |>
    merge(          n_g2_dt,    by = c("id", "Group"), all.x = TRUE, sort = FALSE) |>
    merge(          n_both_dt,  by = c("id", "Group"), all.x = TRUE, sort = FALSE)
  
  cross_replication_summary_dt[, `:=`(
    n_total = fifelse(is.na(n_total), 0L, n_total),
    n_g1    = fifelse(is.na(n_g1),    0L, n_g1),
    n_g2    = fifelse(is.na(n_g2),    0L, n_g2),
    n_both  = fifelse(is.na(n_both),  0L, n_both)
  )]
  
  as.data.frame(cross_replication_summary_dt)
}

safe_binom_greater <- function(x, n, p) {
  vals <- c(x, n, p)
  if (any(is.na(vals))) return(NA_real_)
  if (n <= 0) return(NA_real_)
  if (p <= 0)  return(ifelse(x > 0, 0, 1))
  if (p >= 1)  return(1)
  suppressWarnings(stats::binom.test(x, n, p = p, alternative = "greater")$p.value)
}

safe_fisher_greater <- function(case_yes, case_no, cont_yes, cont_no) {
  vals <- c(case_yes, case_no, cont_yes, cont_no)
  if (any(is.na(vals))) return(NA_real_)
  m <- matrix(vals, nrow = 2, byrow = TRUE)
  if (sum(m) == 0) return(NA_real_)
  suppressWarnings(stats::fisher.test(m, alternative = "greater")$p.value)
}

safe_prop_twosided <- function(x1, n1, x2, n2) {
  vals <- c(x1, n1, x2, n2)
  if (any(is.na(vals))) return(NA_real_)
  if (n1 <= 0 || n2 <= 0) return(NA_real_)
  x1 <- max(min(x1, n1), 0)
  x2 <- max(min(x2, n2), 0)
  suppressWarnings(stats::prop.test(c(x1, x2), c(n1, n2), correct = TRUE)$p.value)
}

cohen_h <- function(p1, p2) {
  2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))
}

calculate_power_safe <- function(p1, p2, n_case, n_control, alpha = 0.05) {
  if (any(is.na(c(p1, p2, n_case, n_control))) || n_case < 2 || n_control < 2) {
    return(list(effect_size = NA_real_, power = NA_real_))
  }
  p1c <- min(max(p1, 0), 1)
  p2c <- min(max(p2, 0), 1)
  h <- cohen_h(p1c, p2c)
  pw <- tryCatch(
    pwr.2p2n.test(h = h, n1 = n_case, n2 = n_control, sig.level = alpha)$power,
    error = function(e) NA_real_
  )
  list(effect_size = h, power = pw)
}


############################################################
## 1) EAS pairs in EUR
############################################################

## EUR input for gene-by-sample table
load("rarecomb/eur_candidates.rds")  # contains: eur_input, eur_candidate
gene_by_sample_eur <- eur_input %>%
  dplyr::filter(!is.na(gene_symbol), Group %in% c("Case", "Control")) %>%
  dplyr::distinct(s, gene_symbol, Group)

## EAS candidate pairs
load("rarecomb/eas_candidates.rds")  # contains: eas_input, eas_candidate
eas_candidate <- eas_candidate %>%
  dplyr::mutate(
    id = paste0(
      pmin(gene_symbol_1, gene_symbol_2), "-",
      pmax(gene_symbol_1, gene_symbol_2)
    )
  )

setDT(eas_candidate)
setDT(gene_by_sample_eur)

cross_replication_summary_eas_in_eur <-
  build_cross_replication_summary(eas_candidate, gene_by_sample_eur)

eas_pairs_in_eur <- cross_replication_summary_eas_in_eur %>%
  dplyr::group_by(id, Group) %>%
  dplyr::summarise(
    N_total = n_total[1L],
    N_g1    = n_g1[1L],
    N_g2    = n_g2[1L],
    N_both  = n_both[1L],
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    N_total = tidyr::replace_na(N_total, 0L),
    N_g1    = tidyr::replace_na(N_g1,    0L),
    N_g2    = tidyr::replace_na(N_g2,    0L),
    N_both  = tidyr::replace_na(N_both,  0L),
    p_exp   = dplyr::if_else(N_total > 0, (N_g1 / N_total) * (N_g2 / N_total), NA_real_),
    p_obs   = dplyr::if_else(N_total > 0, N_both / N_total, NA_real_),
    Exp_Count_Combo = dplyr::if_else(is.na(p_exp) | N_total == 0L, NA_real_, N_total * p_exp),
    Obs_Count_Combo = N_both,
    pvalue_more     = purrr::map2_dbl(
      N_both, N_total,
      ~ safe_binom_greater(.x, .y, p = p_exp[cur_group_id()])
    ),
    fold_enrichment = dplyr::if_else(!is.na(p_exp) & p_exp > 0, p_obs / p_exp, NA_real_),
    log2_enrichment = dplyr::if_else(is.finite(fold_enrichment), log2(fold_enrichment), NA_real_)
  ) %>%
  dplyr::mutate(Group = factor(Group, levels = c("Case", "Control"))) %>%
  tidyr::pivot_wider(
    id_cols    = id,
    names_from = Group,
    values_from = c(
      N_total, N_g1, N_g2,
      Obs_Count_Combo, Exp_Count_Combo, pvalue_more,
      fold_enrichment, log2_enrichment
    ),
    names_glue = "{Group}_{.value}"
  ) %>%
  dplyr::mutate(
    dplyr::across(
      c(
        Case_N_total, Control_N_total,
        Case_N_g1, Control_N_g1, Case_N_g2, Control_N_g2,
        Case_Obs_Count_Combo, Control_Obs_Count_Combo,
        Case_Exp_Count_Combo, Control_Exp_Count_Combo
      ),
      ~ tidyr::replace_na(., 0)
    )
  ) %>%
  dplyr::mutate(
    Case_Obs_count_g1 = Case_N_g1,
    Case_Obs_count_g2 = Case_N_g2,
    Cont_Obs_count_g1 = Control_N_g1,
    Cont_Obs_count_g2 = Control_N_g2,
    Case_Prob_count_g1 = dplyr::if_else(
      Case_N_total > 0, Case_N_g1 / Case_N_total, NA_real_
    ),
    Case_Prob_count_g2 = dplyr::if_else(
      Case_N_total > 0, Case_N_g2 / Case_N_total, NA_real_
    ),
    Cont_Prob_count_g1 = dplyr::if_else(
      Control_N_total > 0, Control_N_g1 / Control_N_total, NA_real_
    ),
    Cont_Prob_count_g2 = dplyr::if_else(
      Control_N_total > 0, Control_N_g2 / Control_N_total, NA_real_
    )
  ) %>%
  dplyr::mutate(
    Case_Exp_Prob_Combo = dplyr::if_else(
      Case_N_total > 0, Case_Exp_Count_Combo / Case_N_total, NA_real_
    ),
    Cont_Exp_Prob_Combo = dplyr::if_else(
      Control_N_total > 0, Control_Exp_Count_Combo / Control_N_total, NA_real_
    ),
    Case_Obs_Prob_Combo = dplyr::if_else(
      Case_N_total > 0, Case_Obs_Count_Combo / Case_N_total, NA_real_
    ),
    Cont_Obs_Prob_Combo = dplyr::if_else(
      Control_N_total > 0, Control_Obs_Count_Combo / Control_N_total, NA_real_
    )
  ) %>%
  dplyr::mutate(
    Case_NonCombo    = pmax(Case_N_total    - Case_Obs_Count_Combo,    0),
    Control_NonCombo = pmax(Control_N_total - Control_Obs_Count_Combo, 0),
    CC_pvalue_fisher_greater = purrr::pmap_dbl(
      list(Case_Obs_Count_Combo, Case_NonCombo,
           Control_Obs_Count_Combo, Control_NonCombo),
      ~ safe_fisher_greater(..1, ..2, ..3, ..4)
    ),
    CC_pvalue_prop_twosided = purrr::pmap_dbl(
      list(Case_Obs_Count_Combo, Case_N_total,
           Control_Obs_Count_Combo, Control_N_total),
      ~ safe_prop_twosided(..1, ..2, ..3, ..4)
    )
  ) %>%
  dplyr::mutate(
    Case_Obs_Prob_Combo = tidyr::replace_na(Case_Obs_Prob_Combo, 0),
    Cont_Obs_Prob_Combo = tidyr::replace_na(Cont_Obs_Prob_Combo, 0),
    Case_Prob_safe      = pmin(pmax(Case_Obs_Prob_Combo, 1e-8), 1 - 1e-8),
    Cont_Prob_safe      = pmin(pmax(Cont_Obs_Prob_Combo, 1e-8), 1 - 1e-8)
  ) %>%
  dplyr::mutate(
    res = purrr::pmap(
      list(Case_Prob_safe, Cont_Prob_safe,
           Case_N_total, Control_N_total),
      ~ calculate_power_safe(..1, ..2, ..3, ..4, alpha = 0.05)
    ),
    Effect_Size = purrr::map_dbl(res, "effect_size"),
    Power       = purrr::map_dbl(res, "power")
  ) %>%
  dplyr::select(-res) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    CC_FDR_fisher_greater = dplyr::if_else(
      is.na(CC_pvalue_fisher_greater),
      NA_real_,
      p.adjust(CC_pvalue_fisher_greater, method = "BH")
    ),
    CC_FDR_prop_twosided = dplyr::if_else(
      is.na(CC_pvalue_prop_twosided),
      NA_real_,
      p.adjust(CC_pvalue_prop_twosided, method = "BH")
    ),
    Effect_Size = dplyr::if_else(
      is.na(Effect_Size),
      2 * asin(sqrt(Case_Prob_safe)) -
        2 * asin(sqrt(Cont_Prob_safe)),
      Effect_Size
    )
  )


############################################################
## 2) EUR pairs in EAS
############################################################

gene_by_sample_eas <- eas_input %>%
  dplyr::filter(!is.na(gene_symbol), Group %in% c("Case", "Control")) %>%
  dplyr::distinct(s, gene_symbol, Group)

eur_candidate <- eur_candidate %>%
  dplyr::mutate(
    id = paste0(
      pmin(gene_symbol_1, gene_symbol_2), "-",
      pmax(gene_symbol_1, gene_symbol_2)
    )
  )

setDT(eur_candidate)
setDT(gene_by_sample_eas)

cross_replication_summary_eur_in_eas <-
  build_cross_replication_summary(eur_candidate, gene_by_sample_eas)

eur_pairs_in_eas <- cross_replication_summary_eur_in_eas %>%
  dplyr::group_by(id, Group) %>%
  dplyr::summarise(
    N_total = n_total[1L],
    N_g1    = n_g1[1L],
    N_g2    = n_g2[1L],
    N_both  = n_both[1L],
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    N_total = tidyr::replace_na(N_total, 0L),
    N_g1    = tidyr::replace_na(N_g1,    0L),
    N_g2    = tidyr::replace_na(N_g2,    0L),
    N_both  = tidyr::replace_na(N_both,  0L),
    p_exp   = dplyr::if_else(N_total > 0, (N_g1 / N_total) * (N_g2 / N_total), NA_real_),
    p_obs   = dplyr::if_else(N_total > 0, N_both / N_total, NA_real_),
    Exp_Count_Combo = dplyr::if_else(is.na(p_exp) | N_total == 0L, NA_real_, N_total * p_exp),
    Obs_Count_Combo = N_both,
    pvalue_more     = purrr::map2_dbl(
      N_both, N_total,
      ~ safe_binom_greater(.x, .y, p = p_exp[cur_group_id()])
    ),
    fold_enrichment = dplyr::if_else(!is.na(p_exp) & p_exp > 0, p_obs / p_exp, NA_real_),
    log2_enrichment = dplyr::if_else(is.finite(fold_enrichment), log2(fold_enrichment), NA_real_)
  ) %>%
  dplyr::mutate(Group = factor(Group, levels = c("Case", "Control"))) %>%
  tidyr::pivot_wider(
    id_cols    = id,
    names_from = Group,
    values_from = c(
      N_total, N_g1, N_g2,
      Obs_Count_Combo, Exp_Count_Combo, pvalue_more,
      fold_enrichment, log2_enrichment
    ),
    names_glue = "{Group}_{.value}"
  ) %>%
  dplyr::mutate(
    dplyr::across(
      c(
        Case_N_total, Control_N_total,
        Case_N_g1, Control_N_g1, Case_N_g2, Control_N_g2,
        Case_Obs_Count_Combo, Control_Obs_Count_Combo,
        Case_Exp_Count_Combo, Control_Exp_Count_Combo
      ),
      ~ tidyr::replace_na(., 0)
    )
  ) %>%
  dplyr::mutate(
    Case_Obs_count_g1 = Case_N_g1,
    Case_Obs_count_g2 = Case_N_g2,
    Cont_Obs_count_g1 = Control_N_g1,
    Cont_Obs_count_g2 = Control_N_g2,
    Case_Prob_count_g1 = dplyr::if_else(
      Case_N_total > 0, Case_N_g1 / Case_N_total, NA_real_
    ),
    Case_Prob_count_g2 = dplyr::if_else(
      Case_N_total > 0, Case_N_g2 / Case_N_total, NA_real_
    ),
    Cont_Prob_count_g1 = dplyr::if_else(
      Control_N_total > 0, Control_N_g1 / Control_N_total, NA_real_
    ),
    Cont_Prob_count_g2 = dplyr::if_else(
      Control_N_total > 0, Control_N_g2 / Control_N_total, NA_real_
    )
  ) %>%
  dplyr::mutate(
    Case_Exp_Prob_Combo = dplyr::if_else(
      Case_N_total > 0, Case_Exp_Count_Combo / Case_N_total, NA_real_
    ),
    Cont_Exp_Prob_Combo = dplyr::if_else(
      Control_N_total > 0, Control_Exp_Count_Combo / Control_N_total, NA_real_
    ),
    Case_Obs_Prob_Combo = dplyr::if_else(
      Case_N_total > 0, Case_Obs_Count_Combo / Case_N_total, NA_real_
    ),
    Cont_Obs_Prob_Combo = dplyr::if_else(
      Control_N_total > 0, Control_Obs_Count_Combo / Control_N_total, NA_real_
    )
  ) %>%
  dplyr::mutate(
    Case_NonCombo    = pmax(Case_N_total    - Case_Obs_Count_Combo,    0),
    Control_NonCombo = pmax(Control_N_total - Control_Obs_Count_Combo, 0),
    CC_pvalue_fisher_greater = purrr::pmap_dbl(
      list(Case_Obs_Count_Combo, Case_NonCombo,
           Control_Obs_Count_Combo, Control_NonCombo),
      ~ safe_fisher_greater(..1, ..2, ..3, ..4)
    ),
    CC_pvalue_prop_twosided = purrr::pmap_dbl(
      list(Case_Obs_Count_Combo, Case_N_total,
           Control_Obs_Count_Combo, Control_N_total),
      ~ safe_prop_twosided(..1, ..2, ..3, ..4)
    )
  ) %>%
  dplyr::mutate(
    Case_Obs_Prob_Combo = tidyr::replace_na(Case_Obs_Prob_Combo, 0),
    Cont_Obs_Prob_Combo = tidyr::replace_na(Cont_Obs_Prob_Combo, 0),
    Case_Prob_safe      = pmin(pmax(Case_Obs_Prob_Combo, 1e-8), 1 - 1e-8),
    Cont_Prob_safe      = pmin(pmax(Cont_Obs_Prob_Combo, 1e-8), 1 - 1e-8)
  ) %>%
  dplyr::mutate(
    res = purrr::pmap(
      list(Case_Prob_safe, Cont_Prob_safe,
           Case_N_total, Control_N_total),
      ~ calculate_power_safe(..1, ..2, ..3, ..4, alpha = 0.05)
    ),
    Effect_Size = purrr::map_dbl(res, "effect_size"),
    Power       = purrr::map_dbl(res, "power")
  ) %>%
  dplyr::select(-res) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    CC_FDR_fisher_greater = dplyr::if_else(
      is.na(CC_pvalue_fisher_greater),
      NA_real_,
      p.adjust(CC_pvalue_fisher_greater, method = "BH")
    ),
    CC_FDR_prop_twosided = dplyr::if_else(
      is.na(CC_pvalue_prop_twosided),
      NA_real_,
      p.adjust(CC_pvalue_prop_twosided, method = "BH")
    ),
    Effect_Size = dplyr::if_else(
      is.na(Effect_Size),
      2 * asin(sqrt(Case_Prob_safe)) -
        2 * asin(sqrt(Cont_Prob_safe)),
      Effect_Size
    )
  )


############################################################
## 3) Effect sizes and plots
############################################################

ori_res <- dplyr::bind_rows(
  eur_candidate %>% dplyr::mutate(ancestry = "EUR"),
  eas_candidate %>% dplyr::mutate(ancestry = "EAS")
)

upd_res <- dplyr::bind_rows(
  eur_pairs_in_eas %>% dplyr::mutate(ancestry = "EAS"),
  eas_pairs_in_eur %>% dplyr::mutate(ancestry = "EUR")
)

## EUR pairs in EAS
df_eur_pairs <- dplyr::full_join(
  ori_res,
  upd_res %>% dplyr::filter(ancestry == "EAS"),
  by = "id"
) %>%
  dplyr::mutate(
    case_rate      = Case_Obs_Count_Combo.x / Case_Exp_Count_Combo.x,
    ctrl_rate      = Cont_Obs_Count_Combo    / Cont_Exp_Count_Combo,
    upd_case_rate  = Case_Obs_Count_Combo.y  / Case_Exp_Count_Combo.y,
    upd_ctrl_rate  = Control_Obs_Count_Combo / Control_Exp_Count_Combo,
    ES             = (case_rate - ctrl_rate) /
      sqrt(1 / (Case_Exp_Count_Combo.x + 1e-8) +
             1 / (Cont_Exp_Count_Combo + 1e-8)),
    upd_ES         = (upd_case_rate - upd_ctrl_rate) /
      sqrt(1 / (Case_Exp_Count_Combo.y + 1e-8) +
             1 / (Control_Exp_Count_Combo + 1e-8))
  ) %>%
  dplyr::filter(is.finite(ES), is.finite(upd_ES)) %>%
  dplyr::filter(Case_Obs_Count_Combo.y != 0 | Control_Obs_Count_Combo != 0)

fit_eur <- lm(upd_ES ~ ES, data = df_eur_pairs)
summary_fit_eur <- broom::tidy(fit_eur)
b0_eur <- coef(fit_eur)[1]
b1_eur <- coef(fit_eur)[2]
pval_eur <- summary_fit_eur$p.value[2]
p_label_eur <- ifelse(
  pval_eur < 1e-5,
  formatC(pval_eur, format = "e", digits = 2),
  round(pval_eur, 5)
)

p_eur_pairs_in_eas <- ggplot(df_eur_pairs, aes(x = ES, y = upd_ES)) +
  geom_point(col = "gray20", alpha = 0.6, size = 0.5) +
  geom_smooth(method = "lm", color = "steelblue", fill = "lightblue", se = TRUE) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = sprintf("y = %.2f x + %.2f\np = %s", b1_eur, b0_eur, p_label_eur),
    hjust = -0.1, vjust = 1.1, size = 3.5, color = "black"
  ) +
  labs(
    x = "East-Asian ancestry",
    y = "European ancestry",
    title = "Effect size of EUR pairs"
  ) +
  theme_v1

## EAS pairs in EUR
df_eas_pairs <- dplyr::full_join(
  ori_res,
  upd_res %>% dplyr::filter(ancestry == "EUR"),
  by = "id"
) %>%
  dplyr::mutate(
    case_rate      = Case_Obs_Count_Combo.x / Case_Exp_Count_Combo.x,
    ctrl_rate      = Cont_Obs_Count_Combo    / Cont_Exp_Count_Combo,
    upd_case_rate  = Case_Obs_Count_Combo.y  / Case_Exp_Count_Combo.y,
    upd_ctrl_rate  = Control_Obs_Count_Combo / Control_Exp_Count_Combo,
    ES             = (case_rate - ctrl_rate) /
      sqrt(1 / (Case_Exp_Count_Combo.x + 1e-8) +
             1 / (Cont_Exp_Count_Combo + 1e-8)),
    upd_ES         = (upd_case_rate - upd_ctrl_rate) /
      sqrt(1 / (Case_Exp_Count_Combo.y + 1e-8) +
             1 / (Control_Exp_Count_Combo + 1e-8))
  ) %>%
  dplyr::filter(is.finite(ES), is.finite(upd_ES)) %>%
  dplyr::filter(Case_Obs_Count_Combo.y != 0 | Control_Obs_Count_Combo != 0)

fit_eas <- lm(upd_ES ~ ES, data = df_eas_pairs)
summary_fit_eas <- broom::tidy(fit_eas)
b0_eas <- coef(fit_eas)[1]
b1_eas <- coef(fit_eas)[2]
pval_eas <- summary_fit_eas$p.value[2]
p_label_eas <- ifelse(
  pval_eas < 1e-5,
  formatC(pval_eas, format = "e", digits = 2),
  round(pval_eas, 5)
)

p_eas_pairs_in_eur <- ggplot(df_eas_pairs, aes(x = ES, y = upd_ES)) +
  geom_point(col = "gray20", alpha = 0.6, size = 0.5) +
  geom_smooth(method = "lm", color = "steelblue", fill = "lightblue", se = TRUE) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = sprintf("y = %.2f x + %.2f\np = %s", b1_eas, b0_eas, p_label_eas),
    hjust = -0.1, vjust = 1.1, size = 3.5, color = "black"
  ) +
  labs(
    x = "East-Asian ancestry",
    y = "European ancestry",
    title = "Effect size of EAS pairs"
  ) +
  theme_v1

p_combined <- ggarrange(
  p_eur_pairs_in_eas,
  p_eas_pairs_in_eur,
  labels = c("A", "B"),
  font.label = list(size = 18, face = "bold"),
  vjust = 1.1
)

ggsave(
  "figures/cross_ancestry_replication.pdf",
  p_combined,
  width = 6.75,
  height = 3
)


############################################################
## 4) Export tables
############################################################

eur_pairs_in_eas_export <- df_eur_pairs %>%
  dplyr::mutate(FDR = p.adjust(Case_pvalue_more.y, method = "BH")) %>%
  dplyr::select(
    gene_id_1   = Item_1,
    gene_symbol_1,
    gene_id_2   = Item_2,
    gene_symbol_2,
    Case_Obs_Count_gene1 = Case_Obs_count_g1,
    Case_Obs_Count_gene2 = Case_Obs_count_g2,
    Case_Exp_Count_Combo = Case_Exp_Count_Combo.y,
    Case_Obs_Count_Combo = Case_Obs_Count_Combo.y,
    Case_pvalue_more     = Case_pvalue_more.y,
    Case_Adj_Pval_BH     = FDR,
    Cont_Obs_Count_I1    = Cont_Obs_count_g1,
    Cont_Obs_Count_I2    = Cont_Obs_count_g2,
    Cont_Exp_Count_Combo = Control_Exp_Count_Combo,
    Cont_Obs_Count_Combo = Control_Obs_Count_Combo,
    Cont_pvalue_more     = Control_pvalue_more
  )

eas_pairs_in_eur_export <- df_eas_pairs %>%
  dplyr::mutate(FDR = p.adjust(Case_pvalue_more.y, method = "BH")) %>%
  dplyr::select(
    gene_id_1   = Item_1,
    gene_symbol_1,
    gene_id_2   = Item_2,
    gene_symbol_2,
    Case_Obs_Count_gene1 = Case_Obs_count_g1,
    Case_Obs_Count_gene2 = Case_Obs_count_g2,
    Case_Exp_Count_Combo = Case_Exp_Count_Combo.y,
    Case_Obs_Count_Combo = Case_Obs_Count_Combo.y,
    Case_pvalue_more     = Case_pvalue_more.y,
    Case_Adj_Pval_BH     = FDR,
    Cont_Obs_Count_I1    = Cont_Obs_count_g1,
    Cont_Obs_Count_I2    = Cont_Obs_count_g2,
    Cont_Exp_Count_Combo = Control_Exp_Count_Combo,
    Cont_Obs_Count_Combo = Control_Obs_Count_Combo,
    Cont_pvalue_more     = Control_pvalue_more
  )

wb <- createWorkbook()
addWorksheet(wb, "EUR_pairs_in_EAS")
writeData(wb, "EUR_pairs_in_EAS", eur_pairs_in_eas_export)

addWorksheet(wb, "EAS_pairs_in_EUR")
writeData(wb, "EAS_pairs_in_EUR", eas_pairs_in_eur_export)

saveWorkbook(
  wb,
  "tables/cross_ancestry_replication.xlsx",
  overwrite = TRUE
)