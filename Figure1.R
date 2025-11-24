############################################################
## Main Figure 1
############################################################

## Load helper functions
source("path/to/function_definition.R") 

## Packages
library(tidyverse)
library(data.table)
library(ggpubr)
library(cowplot)
library(ggrepel)
library(broom)
library(openxlsx)
library(readxl)
library(biomaRt)

############################################################
## Global settings
############################################################

project_dir <- "path/to/project_root"
setwd(project_dir)

clr_case <- "black"
clr_cntl <- "grey"

############################################################
## Load precomputed burden-test results
############################################################

burden_results_file <- "results/burden_results_DNV_and_rare_variants.txt"
burden_results <- read.delim(burden_results_file, stringsAsFactors = FALSE)

## Candidate gene pairs (EUR / EAS)
eur_candidates_file <- "path/to/eur_candidate_pairs.rds"
eas_candidates_file <- "path/to/eas_candidate_pairs.rds"

load(eur_candidates_file) 
load(eas_candidates_file) 

############################################################
## Helper: plot burden barplot from precomputed results
############################################################
plot_burden_from_results <- function(results, ancestry, test_label, title) {
  dat <- results %>% dplyr::filter(Ancestry == ancestry, Test == test_label)

  plot_dat <- dat %>%
    dplyr::group_by(Group, VariantType) %>%
    dplyr::summarise(mean = mean_value[1], sd = sd[1], nSp = sample_size[1], .groups = "drop") %>%
    dplyr::mutate(
      VariantType = dplyr::if_else(VariantType %in% c("missense", "Damaging missense"), "Damaging missense", "Protein truncating"),
      VariantType = factor(VariantType, levels = c("Protein truncating", "Damaging missense")),
      se = sd / sqrt(nSp)
    )

  annot_dat <- dat %>%
    dplyr::group_by(VariantType) %>%
    dplyr::summarise(
      ymax       = max(mean_value + sd / sqrt(sample_size), na.rm = TRUE),
      odds_ratio = odds_ratio[1],
      padj       = padj[1],
      .groups    = "drop"
    ) %>%
    dplyr::mutate(VariantType = factor(VariantType, levels = c("Protein truncating", "Damaging missense")))

  ymax_global <- max(plot_dat$mean + plot_dat$se, na.rm = TRUE)

  ggplot(plot_dat, aes(Group, mean)) +
    geom_bar(stat = "identity", aes(fill = Group), show.legend = FALSE) +
    geom_segment(aes(x = Group, xend = Group, y = pmax(mean - se, 0), yend = mean + se), col = "black", linewidth = 0.35) +
    geom_segment(data = annot_dat, aes(x = 1, xend = 2, y = ymax + ymax / 15, yend = ymax + ymax / 15), col = "black", linewidth = 0.35, inherit.aes = FALSE) +
    geom_text(
      data = annot_dat,
      aes(
        x = 1.5,
        y = ymax + ymax / 8,
        label = paste0("OR = ", round(odds_ratio, 2), "\n", "P = ", formatC(padj, digits = 1, format = "e"))
      ),
      col = "black",
      size = 10 * 0.3514321,
      inherit.aes = FALSE
    ) +
    facet_grid(~VariantType) +
    scale_fill_manual(values = c(Autism = clr_case, Control = clr_cntl)) +
    labs(x = "", y = "Variants per person (adjusted)", title = title) +
    coord_cartesian(ylim = c(0, ymax_global * 1.3)) +
    theme_v1
}

############################################################
## 1. Burden of rare inherited variants per ancestry
############################################################

## EUR – rare inherited variants
p1a <- plot_burden_from_results(
  results   = burden_results,
  ancestry  = "EUR",
  test_label = "Rare_inherited",
  title     = "Burden of rare inherited variants in EUR"
)

## EAS – rare inherited variants
p1b <- plot_burden_from_results(
  results   = burden_results,
  ancestry  = "EAS",
  test_label = "Rare_inherited",
  title     = "Burden of rare inherited variants in EAS"
)



############################################################
## 2. RareComb co-occurrence and inheritance patterns
############################################################
## Co-occurrence plots: EUR
p_dat_eur <- eur_candidate %>%
  dplyr::select(
    itemset_id,
    gene_symbol_1,
    gene_symbol_2,
    starts_with("Case_Obs_Count_Combo"),
    starts_with("Cont_Obs_Count_Combo"),
    starts_with("Case_Exp_Count_Combo"),
    starts_with("Cont_Exp_Count_Combo")
  ) %>%
  tidyr::pivot_longer(
    cols = -c(itemset_id, gene_symbol_1, gene_symbol_2),
    names_pattern = "(.*)(Obs_Count_Combo|Exp_Count_Combo)$",
    names_to  = c("Group", ".value")
  ) %>%
  dplyr::mutate(
    Group = ifelse(Group == "Case_", "Autism", "Control")
  )

max_eur <- max(p_dat_eur$Obs_Count_Combo, p_dat_eur$Exp_Count_Combo, na.rm = TRUE)

p1c_core <- p_dat_eur %>%
  ggplot(aes(Exp_Count_Combo, Obs_Count_Combo, colour = Group)) +
  geom_point(show.legend = TRUE, size = 2) +
  geom_line(aes(group = itemset_id), alpha = 0.5, show.legend = FALSE, linewidth = 0.05) +
  ggrepel::geom_text_repel(
    data = p_dat_eur %>%
      dplyr::filter(Group == "Autism", Obs_Count_Combo >= 8),
    aes(
      Exp_Count_Combo,
      Obs_Count_Combo,
      label = paste0(gene_symbol_1, "-", gene_symbol_2)
    ),
    box.padding = 0.5,
    max.overlaps = 10,
    show.legend = FALSE,
    size = 10 * 0.3514321
  ) +
  coord_cartesian(xlim = c(0, max_eur), ylim = c(0, max_eur)) +
  labs(
    x = "Expected # of samples with co-occurrence",
    y = "Observed # of samples with co-occurrence",
    title = paste0("Co-occurrence events in ", nrow(eur_candidate), " gene pairs (EUR)")
  ) +
  scale_color_manual(values = c(Autism = clr_case, Control = clr_cntl)) +
  geom_abline(size = 0.35, linetype = "dashed") +
  theme_v1

eur_cand_inh_plot <- eur_cand_inh %>%
  dplyr::mutate(
    Inherit_pattern = ifelse(
      Inherit_pattern == "unkown phase",
      "Unknown phase",
      Inherit_pattern
    )
  ) %>%
  plot_inherit_pattern() +
  labs(x = "Inheritance pattern of co-occurrence events") +
  ylim(c(0, 80)) +
  theme(
    legend.position = c(0.6, 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p1c <- cowplot::ggdraw() +
  cowplot::draw_plot(p1c_core) +
  cowplot::draw_plot(eur_cand_inh_plot, x = 0.5, y = 0.08, width = 0.4, height = 0.35)

## Co-occurrence plots: EAS
p_dat_eas <- eas_candidate %>%
  dplyr::select(
    itemset_id,
    gene_symbol_1,
    gene_symbol_2,
    starts_with("Case_Obs_Count_Combo"),
    starts_with("Cont_Obs_Count_Combo"),
    starts_with("Case_Exp_Count_Combo"),
    starts_with("Cont_Exp_Count_Combo")
  ) %>%
  tidyr::pivot_longer(
    cols = -c(itemset_id, gene_symbol_1, gene_symbol_2),
    names_pattern = "(.*)(Obs_Count_Combo|Exp_Count_Combo)$",
    names_to  = c("Group", ".value")
  ) %>%
  dplyr::mutate(
    Group = ifelse(Group == "Case_", "Autism", "Control")
  )

max_eas <- max(p_dat_eas$Obs_Count_Combo, p_dat_eas$Exp_Count_Combo, na.rm = TRUE)

p1d_core <- p_dat_eas %>%
  ggplot(aes(Exp_Count_Combo, Obs_Count_Combo, colour = Group)) +
  geom_point(show.legend = TRUE, size = 2) +
  geom_line(aes(group = itemset_id), alpha = 0.5, show.legend = FALSE, linewidth = 0.05) +
  ggrepel::geom_text_repel(
    data = p_dat_eas %>%
      dplyr::filter(Group == "Autism", Obs_Count_Combo >= 3),
    aes(
      Exp_Count_Combo,
      Obs_Count_Combo,
      label = paste0(gene_symbol_1, "-", gene_symbol_2)
    ),
    box.padding = 0.5,
    max.overlaps = 30,
    show.legend = FALSE,
    size = 10 * 0.3514321
  ) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 10)) +
  labs(
    x = "Expected # of samples with co-occurrence",
    y = "Observed # of samples with co-occurrence",
    title = paste0("Co-occurrence events in ", nrow(eas_candidate), " gene pairs (EAS)")
  ) +
  scale_color_manual(values = c(Autism = clr_case, Control = clr_cntl)) +
  geom_abline(size = 0.35, linetype = "dashed") +
  theme_v1

eas_cand_inh_plot <- eas_cand_inh %>%
  dplyr::mutate(
    Inherit_pattern = ifelse(
      Inherit_pattern == "unkown phase",
      "Unknown phase",
      Inherit_pattern
    )
  ) %>%
  plot_inherit_pattern() +
  labs(x = "Inheritance pattern of co-occurrence events") +
  theme(
    legend.position = c(0.9, 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p1d <- cowplot::ggdraw() +
  cowplot::draw_plot(p1d_core) +
  cowplot::draw_plot(eas_cand_inh_plot, x = 0.5, y = 0.08, width = 0.27, height = 0.35)



############################################################
## 3. Gene constraint (LOEUF) and gene length
############################################################
constraint_file  <- "path/to/constraint_metrics.txt"
sfari_file       <- "path/to/sfari_geneset.csv"
fu2022_tada_file <- "path/to/Fu2022_TADA_geneset.xlsx"

## Constraint + SFARI + TADA72
constraint_dat <- read.delim(constraint_file) %>% dplyr::select(gene_id, oe_lof_upper)

sfari_dat <- read.csv(sfari_file) %>%
  dplyr::filter(gene.score %in% c(1, 2))
sfari_dat$ensembl.id[sfari_dat$gene.symbol == "ADA"]  <- "ENSG00000105976"
sfari_dat$ensembl.id[sfari_dat$gene.symbol == "MET"]  <- "ENSG00000196839"
sfari_dat$ensembl.id[sfari_dat$gene.symbol == "PHB1"] <- "ENSG00000167085"
sfari_dat$ensembl.id[sfari_dat$gene.symbol == "AR"]   <- "ENSG00000169083"
sfari_ids <- sfari_dat$ensembl.id

fu72_id <- readxl::read_excel(fu2022_tada_file) %>%
  dplyr::filter(ASD72) %>%
  dplyr::pull(gene_id)

## Gene sets matrix
mt <- tibble(gene_id = cdsGenes) %>%
  dplyr::mutate(
    DNV       = as.integer(gene_id %in% fu72_id),
    SFARI     = as.integer(gene_id %in% sfari_ids),
    EUR_pairs = as.integer(gene_id %in% (eur_candidate %>% output_gene_id())),
    EAS_pairs = as.integer(gene_id %in% (eas_candidate %>% output_gene_id())),
    Allgenes  = 1L
  ) %>%
  dplyr::left_join(constraint_dat, by = "gene_id") %>%
  tidyr::pivot_longer(cols = 2:6, names_to = "set", values_to = "value") %>%
  dplyr::filter(value == 1) %>%
  dplyr::select(-value) %>%
  dplyr::mutate(set = factor(set, levels = c("DNV", "SFARI", "EAS_pairs", "EUR_pairs", "Allgenes")))

## Gene length & symbol
gene_length <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "cds_length"),
  filters    = "ensembl_gene_id",
  values     = mt$gene_id,
  mart       = mart
)

gene_symbol <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = mt$gene_id,
  mart       = mart
)

max_cds_length <- gene_length %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::summarise(max_CDS_length = max(cds_length, na.rm = TRUE), .groups = "drop")

mt <- mt %>%
  dplyr::left_join(max_cds_length, by = c("gene_id" = "ensembl_gene_id")) %>%
  dplyr::left_join(gene_symbol,   by = c("gene_id" = "ensembl_gene_id")) %>%
  dplyr::rename(gene = hgnc_symbol)

## LOEUF distribution plot
median_loeuf <- mt %>%
  dplyr::group_by(set) %>%
  dplyr::filter(!is.na(oe_lof_upper)) %>%
  dplyr::summarise(median = stats::median(oe_lof_upper, na.rm = TRUE),
                   n      = dplyr::n(),
                   .groups = "drop")

tmp_loeuf <- mt %>%
  dplyr::left_join(median_loeuf, by = "set") %>%
  dplyr::mutate(label = paste0(set, "\n(n=", n, ")"))

labels_loeuf <- tmp_loeuf %>%
  dplyr::distinct(set, label) %>%
  { stats::setNames(.$label, .$set) }

wtest_p <- function(a, b) stats::wilcox.test(a, b)$p.value

p_values_loeuf <- c(
  wtest_p(mt$oe_lof_upper[mt$set == "DNV"],
          mt$oe_lof_upper[mt$set == "Allgenes" & !(mt$gene_id %in% mt$gene_id[mt$set == "DNV"])]),
  wtest_p(mt$oe_lof_upper[mt$set == "SFARI"],
          mt$oe_lof_upper[mt$set == "Allgenes" & !(mt$gene_id %in% mt$gene_id[mt$set == "SFARI"])]),
  wtest_p(mt$oe_lof_upper[mt$set == "EAS_pairs"],
          mt$oe_lof_upper[mt$set == "Allgenes" & !(mt$gene_id %in% mt$gene_id[mt$set == "EAS_pairs"])]),
  wtest_p(mt$oe_lof_upper[mt$set == "EUR_pairs"],
          mt$oe_lof_upper[mt$set == "Allgenes" & !(mt$gene_id %in% mt$gene_id[mt$set == "EUR_pairs"])]),
  wtest_p(mt$oe_lof_upper[mt$set == "EAS_pairs"],
          mt$oe_lof_upper[mt$set == "DNV"]),
  wtest_p(mt$oe_lof_upper[mt$set == "EUR_pairs"],
          mt$oe_lof_upper[mt$set == "DNV"])
) %>% p.adjust(method = "BH")

p1e <- tmp_loeuf %>%
  dplyr::filter(set %in% c("DNV", "SFARI", "EAS_pairs", "EUR_pairs", "Allgenes")) %>%
  ggplot() +
  geom_violin(aes(set, oe_lof_upper, fill = set), show.legend = FALSE, size = 0.1, width = 1) +
  scale_x_discrete(labels = labels_loeuf) +
  scale_fill_manual(values = c("DNV" = "#67322e",
                               "SFARI" = "#c38f16",
                               "EAS_pairs" = "#4e6d58",
                               "EUR_pairs" = "#41507b",
                               "Allgenes" = "grey")) +
  labs(x = "", y = "LOEUF score", title = "Distribution of gene constraint across gene sets") +
  stat_summary(aes(set, oe_lof_upper, fill = set), fun = stats::median,
               geom = "crossbar", size = 0.1, width = 0.5, color = "black", show.legend = FALSE) +
  coord_cartesian(ylim = c(NA, 2.85)) +
  theme_v1 +
  geom_segment(data = tibble(a = 1), x = 1, xend = 5, y = 2.1, yend = 2.1, size = 0.25) +
  geom_text(data = tibble(a = 1), x = 3,   y = 2.15, label = paste0("P = ", formatC(p_values_loeuf[1], format = "e", digits = 1))) +
  geom_segment(data = tibble(a = 1), x = 2, xend = 5, y = 2.3, yend = 2.3, size = 0.25) +
  geom_text(data = tibble(a = 1), x = 3.5, y = 2.35, label = paste0("P = ", formatC(p_values_loeuf[2], format = "e", digits = 1))) +
  geom_segment(data = tibble(a = 1), x = 3, xend = 5, y = 2.5, yend = 2.5, size = 0.25) +
  geom_text(data = tibble(a = 1), x = 4,   y = 2.55, label = paste0("P = ", formatC(p_values_loeuf[3], format = "e", digits = 1))) +
  geom_segment(data = tibble(a = 1), x = 4, xend = 5, y = 2.7, yend = 2.7, size = 0.25) +
  geom_text(data = tibble(a = 1), x = 4.5, y = 2.75, label = paste0("P = ", formatC(p_values_loeuf[4], format = "e", digits = 1))) +
  geom_segment(data = tibble(a = 1), x = 1, xend = 3, y = 2.6, yend = 2.6, size = 0.25) +
  geom_text(data = tibble(a = 1), x = 2,   y = 2.65, label = paste0("P = ", formatC(p_values_loeuf[5], format = "e", digits = 1))) +
  geom_segment(data = tibble(a = 1), x = 1, xend = 4, y = 2.8, yend = 2.8, size = 0.25) +
  geom_text(data = tibble(a = 1), x = 2.5, y = 2.85, label = paste0("P = ", formatC(p_values_loeuf[6], format = "e", digits = 1)))


############################################################
## 4. Final multi-panel Figure 1
############################################################
## Left panel (A–F)
figure1_left <- ggarrange(
  ggarrange(p1b, p1a, nrow = 1, labels = c("A", "B"), font.label = list(size = 18)),
  p1e,
  ncol = 1,
  labels = c(NA, "C"),
  font.label = list(size = 18)
)

## Right panel (D–E)
figure1_right <- ggarrange(
  p1d, p1c,
  ncol = 1,
  labels = c("D", "E"),
  font.label = list(size = 18)
)

## Combine (Left | Right)
figure1 <- ggarrange(
  figure1_left, figure1_right,
  nrow = 1,
  widths = c(2, 1)
)

ggsave("Figure1_main.pdf", figure1, width = 13, height = 11)