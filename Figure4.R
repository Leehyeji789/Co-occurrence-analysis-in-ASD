############################################################
## Main Figure 4
############################################################

## Load helper functions
source("path/to/function_definition.R")

## Packages
library(tidyverse)
library(ggpubr)
library(egg)

############################################################
## Global settings
############################################################

project_dir <- "path/to/project_root"
setwd(project_dir)

################################################################################
# LOAD ANALYSIS RESULTS
################################################################################

load('path/to/output/Phenotypic_association_results_AllCases.rds')
load('path/to/output/Phenotypic_association_results_MaleCases.rds')
load('path/to/output/Phenotypic_association_results_Parents.rds')
load('path/to/output/Case_data_male_biparental.rds')

domains <- c('ADOS_total', 'ADOS_SA', 'ADOS_RRB', 'SCQ_current', 'SCQ_lifetime', 'SRS', 'RBSR', 'FSIQ', 'Non_verbal_IQ', 'VABS', 'DCDQ', 'first_word', 'first_walk')

################################################################################
# PANEL A: HEATMAP OF ASSOCIATIONS IN ALL CASES
################################################################################

# Combine all cases and male/female stratified results
p3a <- bind_rows(res_all %>% mutate(set = 'All'), res_male %>% mutate(set = 'Male')) %>%
  merge(., cross_join(tibble(phe = domains), tibble(gen = c('isDNV', 'PS', 'isComb', 'Co-occurring rare variants\n(Maternal + Paternal)', 'Co-occurring rare variants\n(Uni-parental)'))), by = c('phe', 'gen'), all.y = T) %>%
  mutate(pval = ifelse((n_1 < 5 & !is.na(n_1)) | n < 30, 1, pval)) %>%
  filter(phe %in% c('ADOS_total', 'SCQ_current', 'SCQ_lifetime', 'SRS', 'ADOS_SA', 'ADOS_RRB', 'RBSR', 'FSIQ', 'Non_verbal_IQ', 'VABS', 'DCDQ', 'first_word', 'first_walk')) %>%
  mutate(gen = ifelse(gen == 'isDNV', 'DNV', ifelse(gen == 'isComb', 'Co-occurring rare variants', gen)), phe = factor(phe, level = c('ADOS_total', 'ADOS_SA', 'ADOS_RRB', 'SCQ_current', 'SCQ_lifetime', 'SRS', 'RBSR', 'FSIQ', 'Non_verbal_IQ', 'VABS', 'DCDQ', 'first_word', 'first_walk')), estimate = ifelse(estimate > 0.5, 0.5, ifelse(estimate < -0.5, -0.5, estimate)), pval = ifelse(pval < 0.0001, 0.0001, pval), sig = ifelse(pval < 0.05, 'Y', 'N')) %>%
  filter(gen %in% c('DNV', 'PS', 'Co-occurring rare variants')) %>%
  filter(ifelse(set %in% c('Male'), gen %in% c('Co-occurring rare variants'), TRUE)) %>%
  mutate(gen = factor(ifelse(set %in% c('Male'), paste0(gen, ' (', set, ')'), gen), levels = rev(c('DNV', 'PS', 'Co-occurring rare variants', 'Co-occurring rare variants (Male)')))) %>%
  ggplot(aes()) +
  geom_tile(aes(phe, gen), fill = 'white', size = 0.25, color = 'black', show.legend = F) +
  geom_point(aes(phe, gen, size = -log10(pval), fill = estimate), shape = 22) +
  geom_point(aes(phe, gen, shape = sig), size = 8, show.legend = F) +
  scale_shape_manual(values = c('Y' = 8, 'N' = NA)) +
  scale_size_continuous(range = c(0, 25), breaks = c(0.5, 1, 2, 4)) +
  scale_fill_gradient2(high = "#6c1d0e", low = "#235070", midpoint = 0, limits = c(-0.55, 0.55)) +
  theme_minimal(base_size = 10) +
  theme(legend.background = element_blank(), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x.bottom = element_blank(), axis.line.y = element_blank(), axis.ticks = element_blank(), panel.border = element_blank(), axis.title = element_text(size = 10), axis.text = element_text(size = 10), plot.title = element_text(size = 12)) +
  labs(x = '', y = '', title = 'Association between genetic factors and phenotypes in the autism case group') +
  coord_fixed(ratio = 1)

################################################################################
# PANEL B: ODDS RATIO PLOT FOR MALE CASES
################################################################################

tmp1 <- res_male %>%
  merge(., cross_join(tibble(phe = domains), tibble(gen = c('isDNV', 'PS', 'isComb', 'Co-occurring rare variants\n(Maternal + Paternal)', 'Co-occurring rare variants\n(Uni-parental)'))), by = c('phe', 'gen'), all.y = T) %>%
  mutate(pval = ifelse((n_1 < 5 & !is.na(n_1)) | n < 30, 1, pval)) %>%
  filter(phe %in% c('ADOS_total', 'SCQ_current', 'SCQ_lifetime', 'SRS', 'ADOS_SA', 'ADOS_RRB', 'RBSR', 'FSIQ', 'Non_verbal_IQ', 'VABS', 'DCDQ', 'first_word', 'first_walk')) %>%
  mutate(gen = factor(ifelse(gen == 'isDNV', 'DNV', ifelse(gen == 'isComb', 'Co-occurring rare variants', gen)), levels = rev(c('DNV', 'PS', 'Co-occurring rare variants', 'Co-occurring rare variants\n(Maternal + Paternal)', 'Co-occurring rare variants\n(Uni-parental)'))), phe = factor(phe, level = c('ADOS_total', 'ADOS_SA', 'ADOS_RRB', 'SCQ_current', 'SCQ_lifetime', 'SRS', 'RBSR', 'FSIQ', 'Non_verbal_IQ', 'VABS', 'DCDQ', 'first_word', 'first_walk')), estimate = ifelse(estimate > 0.5, 0.5, ifelse(estimate < -0.5, -0.5, estimate)), pval = ifelse(pval < 0.0001, 0.0001, pval), sig = ifelse(pval < 0.05, 'Y', 'N')) %>%
  filter(gen %in% c('Co-occurring rare variants', 'Co-occurring rare variants\n(Maternal + Paternal)', 'Co-occurring rare variants\n(Uni-parental)') & phe %in% c('ADOS_total')) %>%
  mutate(gen = factor(gen, levels = c('Co-occurring rare variants', 'Co-occurring rare variants\n(Maternal + Paternal)', 'Co-occurring rare variants\n(Uni-parental)')))

label <- paste0(tmp1$gen, '\n(n = ', tmp1$n_1, ')') %>% unique()

p3b <- tmp1 %>%
  ggplot() +
  geom_point(aes(gen, or), size = 2.5) +
  geom_segment(aes(y = l_ci, yend = u_ci, x = gen, xend = gen)) +
  geom_point(aes(gen, or, fill = sig), size = 3, shape = 21, show.legend = F) +
  geom_hline(yintercept = 1, linetype = 'dashed', size = .35) +
  scale_fill_manual(values = c('Y' = 'black', 'N' = 'white')) +
  scale_x_discrete(labels = label) +
  theme_v1 +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(y = 'Odds ratio', x = '')

################################################################################
# PANEL C: HEATMAP OF ASSOCIATIONS IN PARENTS
################################################################################

p4a <- res_par %>%
  merge(., cross_join(tibble(phe = c('BAPQ_overall', 'BAPQ_aloof', 'BAPQ_PL', 'BAPQ_rigidity', 'SRS')), tibble(gen = c('isChild_dnv_carrier', 'PS', 'isComb', 'isChild_combo_carrier'))), by = c('phe', 'gen'), all.y = T) %>%
  mutate(pval = ifelse((n_1 < 5 & !is.na(n_1)) | n < 30, 1, pval)) %>%
  filter(phe %in% c('BAPQ_overall', 'BAPQ_aloof', 'BAPQ_PL', 'BAPQ_rigidity', 'SRS')) %>%
  mutate(gen = factor(case_when(gen == 'isDNV' ~ 'DNV', gen == 'isComb' ~ 'Co-occurring rare variants', gen == 'isChild_dnv_carrier' ~ 'Parents of cases\nwith DNV', gen == 'isChild_combo_carrier' ~ 'Parents of cases\nwith co-occurring rare variants', gen == 'PS' ~ 'PS'), levels = rev(c('PS', 'Co-occurring rare variants', 'Parents of cases\nwith DNV', 'Parents of cases\nwith co-occurring rare variants'))), phe = factor(phe, level = c('BAPQ_overall', 'BAPQ_aloof', 'BAPQ_PL', 'BAPQ_rigidity', 'SRS')), pval = ifelse(pval < 0.0001, 0.0001, pval), sig = ifelse(pval < 0.05, 'Y', 'N')) %>%
  ggplot(aes()) +
  geom_tile(aes(phe, gen), fill = 'white', size = 0.25, color = 'black', show.legend = F) +
  geom_point(aes(phe, gen, size = -log10(pval), fill = estimate), shape = 22) +
  geom_point(aes(phe, gen, shape = sig), size = 8, show.legend = F) +
  scale_shape_manual(values = c('Y' = 8, 'N' = NA)) +
  scale_size_continuous(range = c(0, 25), breaks = c(0.5, 1, 2, 4)) +
  scale_fill_gradient2(high = "#6c1d0e", low = "#235070", midpoint = 0, limits = c(-0.55, 0.55)) +
  theme_minimal(base_size = 10, base_family = 'Arial') +
  theme(legend.background = element_blank(), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x.bottom = element_blank(), axis.line.y = element_blank(), axis.ticks = element_blank(), panel.border = element_blank(), axis.title = element_text(size = 10), axis.text = element_text(size = 10), plot.title = element_text(size = 12)) +
  labs(x = '', y = '', title = 'Association between genetic factors and phenotypes in the control parents group') +
  coord_fixed(ratio = 1)

################################################################################
# COMBINE PANELS INTO FIGURE 4
################################################################################

p4 <- ggarrange(
  p3a, 
  ggarrange(p3b, p4a, labels = c('B', 'C'), nrow = 1, font.label = list(size = 18), widths = c(.75, 1)),
  labels = c('A', NA), 
  ncol = 1, 
  font.label = list(size = 18), 
  heights = c(1.2, 1)
)

################################################################################
# SAVE FIGURE
################################################################################

ggsave('Figure4_main.pdf', p4, width = 14, height = 9, dpi = 300, device = cairo_pdf)
cat("Figure 4 generated and saved successfully.\n")