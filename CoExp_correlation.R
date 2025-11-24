stringsAsFactors = FALSE

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(ggrepel))
suppressMessages(library(anndata))
library(reticulate)
library(tidyr)
library(writexl)
library(foreach)
library(doParallel)
library(cowplot)
library(readxl)
library(data.table)

##------------------------------------------------------------
## Settings
##------------------------------------------------------------

# Number of cores
num_cores <- 36
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Python path (modify as appropriate)
use_python("path/to/miniforge3/bin/python")

# Output directory for correlation results
o_dir <- "path/to/Cor_Exp_genepair/"

# Gene-pair list
genepair <- read.delim("path/to/gene_pair_list.txt", header = TRUE, stringsAsFactors = FALSE)

##------------------------------------------------------------
## Load pseudobulk data
##------------------------------------------------------------

ad      <- anndata::read_h5ad("path/to/Pseudobulk_mincells_30_lognorm_mean_0327.h5ad")
ad.raw  <- as.data.frame(ad$X)        # pseudobulk expression
ad.meta <- as.data.frame(ad[["obs"]]) # metadata

# Exclude dataset with unspecific age
ad.meta  <- ad.meta[ad.meta$Dataset != "Nagy", ]

# Subset clusters used in dynamic plot
ad.meta2 <- ad.meta[ad.meta$leiden_0.6 %in% c(0:28, 30:39), ]

##------------------------------------------------------------
## Function to process one gene pair
##------------------------------------------------------------

process_gene_pair <- function(gene1, gene2) {
  message(paste0("Start ", gene1, "-", gene2))
  
  # Create folder for this gene pair
  output_folder <- paste0(gene1, "-", gene2)
  out_dir_pair  <- file.path(o_dir, output_folder)
  if (!dir.exists(out_dir_pair)) dir.create(out_dir_pair, recursive = TRUE)
  
  res <- data.frame()
  
  # Iterate over clusters
  clusters_sorted <- sort(unique(ad.meta2$leiden_0.6))
  
  for (cl_id in clusters_sorted) {
    # Subset data for current cluster
    ad.meta_sub <- ad.meta2[ad.meta2$leiden_0.6 == cl_id, ]
    ad.raw_sub  <- ad.raw[rownames(ad.meta_sub), ]
    
    # Subset two genes of interest
    ad.raw_sub_gene <- ad.raw_sub[, c(gene1, gene2)]
    
    # Extract cell-type annotations
    cell_type        <- unique(ad.meta_sub$`Cell Type`)
    leiden_label     <- unique(ad.meta_sub$Leiden)
    cluster_annot    <- unique(ad.meta_sub$cluster_annotated)
    
    # Correlation test
    cor_res <- cor.test(ad.raw_sub_gene[, gene1], ad.raw_sub_gene[, gene2])
    
    # Collect results
    res <- dplyr::bind_rows(
      res,
      data.frame(
        Gene1        = gene1,
        Gene2        = gene2,
        Leiden       = unique(ad.meta_sub$leiden_0.6),
        Cell_Type    = cell_type,
        Cell_Subtype = cluster_annot,
        correlation  = cor_res$estimate,
        pvalue       = cor_res$p.value
      )
    )
  }
  
  # Save result table for this gene pair
  out_file <- file.path(
    out_dir_pair,
    paste0("Table_BTS_psbulk_leiden_corr_", gene1, "_", gene2, ".txt")
  )
  write.table(
    res,
    out_file,
    sep       = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote     = FALSE
  )
  
  # Clean up objects to free memory
  rm(output_folder, out_dir_pair, ad.raw_sub, ad.meta_sub,
     ad.raw_sub_gene, cor_res, res, cell_type, leiden_label, cluster_annot)
  gc()
}

##------------------------------------------------------------
## Run in parallel for each gene pair
##------------------------------------------------------------

foreach(i = 1:nrow(genepair),
        .packages = c("anndata", "ggplot2", "dplyr", "tidyr",
                      "openxlsx", "ggrepel", "reticulate",
                      "writexl", "cowplot")) %dopar% {
  gene1 <- genepair[i, 1]
  gene2 <- genepair[i, 2]
  process_gene_pair(gene1, gene2)
}

# Stop cluster
stopCluster(cl)