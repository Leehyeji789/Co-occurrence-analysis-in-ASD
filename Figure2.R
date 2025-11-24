############################################################
## Main Figure 2
############################################################

## Load helper functions
source("path/to/function_definition.R") 

## Packages
library(tidyr)
library(ggrepel)
library(writexl)
library(foreach)
library(doParallel)
library(cowplot)
library(readxl)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)

############################################################
## Global settings
############################################################

project_dir <- "path/to/project_root"
setwd(project_dir)

############################################################
## Load data
############################################################
eur_candidates_file <- "path/to/eur_candidate_pairs.rds"
eas_candidates_file <- "path/to/eas_candidate_pairs.rds"
load(eur_candidates_file) 
load(eas_candidates_file) 

# FET results
fet_results_df = read_delim('path/to/FET_res_Disorder_Geneset.tsv')
fet_results_df_deg = read_delim('path/to/FET_res_CellType_DEGs.tsv')

############################################################
## A. GO enrichment
############################################################
tb_go <- list(
  DNV        = fu72,
  Candidate = c(eur_candidate %>% output_genes(),
                 eas_candidate %>% output_genes())
)

ego_mf <- compareCluster(
  geneClusters  = tb_go,
  OrgDb         = "org.Hs.eg.db",
  fun           = "enrichGO",
  keyType       = "SYMBOL",
  pAdjustMethod = "BH",
  ont           = "MF",
  qvalueCutoff  = 0.05,
  minGSSize     = 5,
  maxGSSize     = 2000,
  readable      = FALSE
)
edox2 <- pairwise_termsim(ego_mf, method = "JC")

## Figure 2A (wide layout)
p2a <- emapplot(
  edox2,
  showCategory       = 9,
  color              = "p.adjust",
  layout             = "nicely",
  cex_line           = 0.35,
  node_scale         = 2,
  cex_label_category = 1.2,
  pie                = "count",
  label_style        = "ggforce",
  label_format       = 10
) +
  scale_fill_manual(values = c("DNV" = "#8C3030", "Candidate" = "#806BA5")) +
  theme_v1 +
  coord_fixed(ratio = 0.5)

############################################################
## B-C. Gene set enrichment test
############################################################
plotFetResultsHeatmap <- function(df, orderA, orderB, title, labsA = "", labsB = "") {
  df$fisher_padj <- p.adjust(df$P.value, method = "BH", n = nrow(df))
  df1 <- df
  df1$Cluster <- factor(df1$Cluster, levels = orderA)
  df1$Disease <- factor(df1$Disease, levels = rev(orderB))
  
  df1$OR    <- ifelse(is.infinite(df1$OddsRatio), 0, df1$OddsRatio)
  df1$log2OR <- log2(df1$OR)
  df1$log2OR <- ifelse(is.infinite(df1$log2OR), 0, df1$log2OR)
  
  lim <- max(abs(df1$OR), na.rm = TRUE) %>% ceiling()
  
  ggplot(df1, aes(Cluster, Disease)) +
    geom_tile(fill = "white", linewidth = 0.25, color = "black", show.legend = FALSE) +
    geom_point(aes(fill = log2OR), shape = 22, size = 1.5, color = "black", stroke = 0.5) +
    geom_tile(
      data = df1[df1$P.value < 0.05, ],
      aes(fill = log2OR),
      height = 1, width = 1, linewidth = 0.25, color = "black", show.legend = FALSE
    ) +
    geom_point(
      data = df1[df1$fisher_padj < 0.05 & df1$OR != 0, ],
      shape = 8, size = 2, color = "black", show.legend = FALSE
    ) +
    scale_fill_gradient2(low = "#6A90CA", mid = "white", high = "#CD2836",
                         limits = c(-log2(lim), log2(lim))) +
    labs(title = title, x = labsA, y = labsB) +
    theme_minimal(base_size = 10) +
    theme(
      axis.ticks        = element_blank(),
      plot.title        = element_text(size = 12),
      panel.border      = element_blank(),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.background  = element_blank(),
      legend.text       = element_text(size = 10),
      axis.text.x       = element_text(size = 10, angle = 90, vjust = 0, colour = "black"),
      axis.text.y       = element_text(size = 10, colour = "black")
    )
}


developmental     <- c("ADHD", "CHD", "DDD", "EP")
neurodegenerative <- c("AD", "ALS", "MS", "PD")
neuropsychiatric  <- c("ANX", "BP", "MDD", "SCZ(G)", "SCZ(R)")
other_disorders   <- c("ABM", "EA", "RTE", "VBI")

orders <- c(developmental, neurodegenerative, neuropsychiatric, other_disorders)

## Panel B: Disorder-associated gene sets
tmp <- fet_results_df
sig_res <- tmp %>% dplyr::filter(padj < 0.05) %>% dplyr::pull(Cluster)
tmp <- tmp %>% dplyr::filter(Cluster %in% c(sig_res, fet_results_df$Cluster))

p2b <- plotFetResultsHeatmap(tmp, rev(orders), rev(unique(tmp$Disease)), title = "Enrichment") +
  coord_flip()


## Panel C: Developing human brain cell-type DEGs
desired_order <- c(
  # RG
  "C4", "C18", "C37",
  # Neuroblast
  "C11", "C33", "C3", "C24", "C26",
  # Excitatory
  "C12", "C29", "C25", "C16", "C17", "C14", "C36", "C22", "C5", "C10", "C1", "C15", "C34", "C35",
  # Inhibitory
  "C20", "C23", "C21", "C28", "C8", "C31", "C6", "C7", "C19",
  # Astro
  "C2", "C27",
  # Micro
  "C39", "C13",
  # Oligo
  "C0", "C38",
  # OPC
  "C9", "C32",
  # Endo
  "C30"
)

p2c <- plotFetResultsHeatmap(
  fet_results_df_deg %>% dplyr::select(-padj),
  rev(desired_order),
  rev(unique(fet_results_df_deg$Disease)),
  title = "Cell-type specific DEGs in developing human brain"
) + coord_flip()


############################################################
## Merge
############################################################
## Final Figure 2
p2 <- ggarrange(
  p2a,
  p2b,
  p2c,
  nrow       = 1,
  labels     = c("A", "B", "C"),
  vjust      = 1.1,
  font.label = list(size = 18),
  widths     = c(1.8, 0.54, 0.5)
)

p2

ggsave("Figure2_main.pdf", p2, width = 13, height = 7)