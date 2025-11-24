############################################################
## Main Figure 3
############################################################

## Load helper functions
source("path/to/function_definition.R")

## Packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(anndata)
library(reticulate)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggpubr)
library(patchwork)

############################################################
## Global settings
############################################################

project_dir <- "path/to/project_root"
setwd(project_dir)

############################################################
## Load data
############################################################

## Candidate gene pairs (if needed elsewhere)
eur_candidates_file <- "path/to/eur_candidate_pairs.rds"
eas_candidates_file <- "path/to/eas_candidate_pairs.rds"
load(eur_candidates_file)
load(eas_candidates_file)

## Correlation results (per disrupted gene pair)
o_dir <- "path/to/Cor_Exp_genepair/"   

## Gene-pair list (two columns: Gene1, Gene2)
genepair <- read.delim(
  "path/to/gene_pair_list.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

## Cell-type metadata
celltype_file <- "path/to/Metadata.csv"

## Pseudobulk expression atlas 
pb_h5ad <- "path/to/Pseudobulk_mincells_30_lognorm_mean_0327.h5ad"

############################################################
## 1. Load cell-type annotation from metadata
############################################################
celltype_df <- read.csv(celltype_file, stringsAsFactors = FALSE)

## Common Leiden order
leiden_order <- c(
  # RG
  "C4","C18","C37",
  # Neuroblast
  "C11","C33","C3","C24","C26",
  # Excitatory
  "C12","C29","C25","C16","C17","C14","C36","C22","C5","C10","C1","C15","C34","C35",
  # Inhibitory
  "C20","C23","C21","C28","C8","C31","C6","C7","C19",
  # OPC
  "C9","C32",
  # Oligo
  "C0","C38",
  # Astro
  "C2","C27",
  # Micro
  "C39","C13",
  # Endo
  "C30"
)

celltype_df <- celltype_df %>%
  mutate(Leiden = factor(Leiden, levels = leiden_order))

############################################################
## 2. Load correlation results
############################################################

load_corr_results <- function(genepair, o_dir) {
  res_list <- lapply(seq_len(nrow(genepair)), function(i) {
    g1 <- genepair$Gene1[i]
    g2 <- genepair$Gene2[i]
    folder <- paste0(g1, "-", g2)
    fpath  <- file.path(
      o_dir,
      folder,
      paste0("Table_BTS_psbulk_leiden_corr_", g1, "_", g2, ".txt")
    )
    
    if (!file.exists(fpath)) {
      message("Missing file: ", fpath)
      return(NULL)
    }
    
    tb <- data.table::fread(fpath)
    
    tb %>%
      mutate(
        ID    = paste0(g1, "-", g2),
        Gene1 = g1,
        Gene2 = g2
      )
  })
  
  bind_rows(res_list) %>%
    mutate(fdr = p.adjust(pvalue, method = "BH"))
}

## Correlation results for disrupted gene pairs
disrupted_res <- load_corr_results(genepair, o_dir)

filtered_data_disrupted <- disrupted_res %>%
  mutate(Leiden = factor(paste0("C", Leiden), levels = leiden_order))

############################################################
## 3. Panel A: barplot
##    Cell-type counts of significant correlations
############################################################

cluster_cols <- c(
  "Astro"      = "#1f77b4ff",
  "Excitatory" = "#ff7f0eff",
  "Inhibitory" = "#d62728ff",
  "Micro"      = "#9467bdff",
  "Oligo"      = "#bdbd22ff",
  "OPC"        = "#17becfff",
  "RG"         = "#ffbb78ff",
  "Neuroblast" = "#e377c2ff",
  "Endo"       = "grey",
  "Others"     = "grey50"
)

sig_counts <- filtered_data_disrupted %>%
  filter(!is.na(correlation), fdr < 0.05, correlation > 0.5) %>%
  count(Leiden)

all_leiden <- tibble(Leiden = factor(leiden_order, levels = leiden_order))

p3a <- all_leiden %>%
  left_join(sig_counts, by = "Leiden") %>%
  mutate(n = ifelse(is.na(n), 0L, n)) %>%
  left_join(
    celltype_df %>% select(Leiden, cluster) %>% distinct(),
    by = "Leiden"
  ) %>%
  ggplot(aes(x = Leiden, y = n, fill = cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cluster_cols, drop = FALSE, na.value = "grey80") +
  labs(
    x     = "Cell type",
    y     = "Number of significant correlations",
    title = "Cell-type specific gene expression correlation of disrupted gene pairs (R > 0.5)"
  ) +
  theme_v1 +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x.bottom = element_text(angle = 90, vjust = 0.5)
  )

############################################################
## 4. Panel B: distribution of #cell-types per disrupted pair
############################################################

p3b <- filtered_data_disrupted %>%
  filter(!is.na(correlation), fdr < 0.05, correlation > 0.5) %>%
  count(ID) %>%
  count(n) %>%
  rename(n_cluster = n, n_pairs = nn) %>%
  mutate(
    n_cluster = ifelse(n_cluster >= 10, "≥10", as.character(n_cluster)),
    n_cluster = factor(
      n_cluster,
      levels = c(as.character(1:9), "≥10")
    )
  ) %>%
  ggplot(aes(x = n_cluster, y = n_pairs)) +
  geom_bar(stat = "identity", fill = "navy") +
  labs(
    x     = "Number of cell types",
    y     = "Number of gene pairs",
    title = "Distribution of disrupted gene pairs by number of cell types\nwith significant positive correlation"
  ) +
  theme_v1

############################################################
## 5. Cluster-based GO for disrupted gene pairs (Panel D)
############################################################

## 5-1. build heatmap matrix for clustering
hm_ids <- filtered_data_disrupted %>%
  filter(fdr < 0.05, correlation > 0.5) %>%
  count(ID) %>%
  filter(n >= 2) %>%
  pull(ID)

heatmap_data_disrupted <- filtered_data_disrupted %>%
  filter(ID %in% hm_ids, fdr < 0.05, correlation > 0.5) %>%
  select(Leiden, ID, correlation) %>%
  mutate(Leiden = factor(Leiden, levels = leiden_order)) %>%
  pivot_wider(names_from = Leiden, values_from = correlation) %>%
  tibble::column_to_rownames("ID") %>%
  as.matrix()

heatmap_data_disrupted[is.na(heatmap_data_disrupted)] <- 0

## 5-2. hierarchical clustering on disrupted pairs (rows)
d  <- dist(heatmap_data_disrupted, method = "manhattan")
hc <- hclust(d, method = "ward.D2")
row_cluster <- cutree(hc, k = 3)

row_anno <- data.frame(cluster = as.character(row_cluster))
rownames(row_anno) <- rownames(heatmap_data_disrupted)

## 5-3. GO enrichment per disrupted-pair cluster
nCluster <- 3
tb_cluster <- lapply(1:nCluster, function(i) {
  row_anno %>%
    mutate(ID = rownames(.)) %>%
    filter(cluster == as.character(i)) %>%
    pull(ID) %>%
    strsplit("-", fixed = TRUE) %>%
    unlist() %>%
    unique()
})
names(tb_cluster) <- paste0("cluster", 1:nCluster)

ego_corr <- compareCluster(
  geneClusters  = tb_cluster,
  OrgDb         = "org.Hs.eg.db",
  fun           = "enrichGO",
  ont           = "ALL",
  keyType       = "SYMBOL",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  minGSSize     = 5,
  maxGSSize     = 2000,
  readable      = FALSE
)

go_percluster <- pairwise_termsim(ego_corr, method = "JC")

## 5-4. Dot plot for a specific cluster (e.g. cluster2) – Panel D
p3d <- go_percluster@compareClusterResult %>%
  filter(Cluster == "cluster2") %>%
  arrange(p.adjust) %>%
  slice_head(n = 15) %>%
  ## Convert GeneRatio from "a/b" to numeric a/b
  mutate(
    GR_num = as.numeric(sub("/.*", "", GeneRatio)),
    GR_den = as.numeric(sub(".*/", "", GeneRatio)),
    GeneRatioVal = GR_num / GR_den
  ) %>%
  ggplot(aes(
    x     = GeneRatioVal,
    y     = reorder(Description, GeneRatioVal),
    size  = Count,
    color = -log10(p.adjust)
  )) +
  geom_point() +
  scale_size_continuous(range = c(3, 8)) +
  labs(
    x     = "Gene ratio",
    y     = "",
    title = "Enriched GO terms in genes from disrupted-pair cluster 2"
  ) +
  theme_v1

############################################################
## 6. Pseudobulk expression & smooth curves (Panels E/F/G)
############################################################

use_python("path/to/miniforge3/bin/python", required = FALSE)

ad      <- anndata::read_h5ad(pb_h5ad)
ad.raw  <- as.data.frame(ad$X)
ad.meta <- as.data.frame(ad[["obs"]])

## Subset clusters used for dynamic plot
ad.meta2 <- ad.meta[ad.meta$leiden_0.6 %in% c(0:28, 30:39), ]

plot_corr_plot <- function(gene1, gene2, leiden_code, corr_df) {
  cl_num <- as.numeric(leiden_code)
  
  meta_sub <- ad.meta2[ad.meta2$leiden_0.6 == cl_num, ]
  meta_sub <- meta_sub %>% arrange(Age)
  
  expr_sub <- ad.raw[rownames(meta_sub), c(gene1, gene2), drop = FALSE]
  
  expr_long <- expr_sub %>%
    mutate(Age = as.numeric(as.character(meta_sub$Age))) %>%
    pivot_longer(
      cols      = c(all_of(gene1), all_of(gene2)),
      names_to  = "Gene",
      values_to = "Expression"
    )
  
  cor_row <- corr_df %>%
    filter(
      Leiden == cl_num,
      (Gene1 == gene1 & Gene2 == gene2) |
        (Gene1 == gene2 & Gene2 == gene1)
    ) %>%
    slice_head(n = 1)
  
  R_est <- cor_row$correlation[1]
  p_adj <- cor_row$fdr[1]
  
  ggplot(expr_long, aes(x = Age, y = Expression, color = Gene)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_smooth(se = TRUE, method = "loess", span = 0.5, linewidth = 1) +
    scale_color_manual(values = c(gene1 = "blue", gene2 = "red")) +
    labs(
      title = paste0("Cluster ", leiden_code),
      x     = "Gestational age (days)",
      y     = "Expression level"
    ) +
    annotate(
      "text",
      x      = Inf,
      y      = Inf,
      label  = paste0("R = ", round(R_est, 2)),
      hjust  = 1.3,
      vjust  = 2.0,
      size   = 3,
      color  = "black"
    ) +
    annotate(
      "text",
      x      = Inf,
      y      = Inf,
      label  = paste0("padj = ", formatC(p_adj, digits = 1, format = "e")),
      hjust  = 1.3,
      vjust  = 3.5,
      size   = 3,
      color  = "black"
    ) +
    theme_classic() +
    theme_v1
}

## Example disrupted gene pairs and clusters for panels E/F/G
## (replace these gene/cluster choices as needed)
p3e <- plot_corr_plot("IFT122", "CPLANE1", "4", filtered_data_disrupted)
p3f <- plot_corr_plot("TIAM2", "WDR19",  "4", filtered_data_disrupted)
p3g <- plot_corr_plot("ANK2",  "BBX",    "4", filtered_data_disrupted)


############################################################
## Merge
############################################################

panel_top <- ggarrange(
  p3a, p3b, p3d,
  nrow       = 1,
  labels     = c("A", "B", "D"),
  vjust      = 1.1,
  font.label = list(size = 18),
  widths     = c(1.4, 0.5, 1.1)
)

panel_bottom <- ggarrange(
  NULL,
  ggarrange(
    p3e, p3f, p3g,
    ncol       = 1,
    labels     = c("E", "F", "G"),
    vjust      = 1.1,
    font.label = list(size = 18)
  ),
  nrow       = 1,
  labels     = c("C", NA),
  vjust      = 1.1,
  font.label = list(size = 18)
)

p3 <- panel_top / panel_bottom + plot_layout(heights = c(1, 2.6))

## Save figure (example)
ggsave(
  "Figure3_main.pdf",
  p3,
  dpi   = 300,
  width = 16,
  height = 11.4
)