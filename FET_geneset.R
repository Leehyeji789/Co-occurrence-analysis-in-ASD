library(dplyr)
library(readr)
library(openxlsx)
## Load helper functions
source("path/to/function_definition.R") 

############################################
## 1. Generic FET function
############################################

run_fet_sets <- function(gene_set_list, disease_gene_dict, bg, padj_method = "BH") {
  res <- lapply(names(gene_set_list), function(set_name) {
    set_genes <- intersect(gene_set_list[[set_name]], bg)
    
    lapply(names(disease_gene_dict), function(ds_name) {
      disease_genes <- intersect(disease_gene_dict[[ds_name]], bg)
      
      overlap <- length(intersect(set_genes, disease_genes))
      a <- overlap
      b <- length(set_genes)      - overlap
      c <- length(disease_genes)  - overlap
      d <- length(bg) - (a + b + c)
      
      ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
      
      data.frame(
        Set              = set_name,
        DiseaseSet       = ds_name,
        P.value          = ft$p.value,
        OddsRatio        = unname(ft$estimate),
        CI_low           = ft$conf.int[1],
        CI_high          = ft$conf.int[2],
        overlap          = a,
        set_only         = b,
        disease_only     = c,
        bg_only          = d,
        overlapped_genes = paste0(intersect(set_genes, disease_genes), collapse = ","),
        stringsAsFactors = FALSE
      )
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  res %>% mutate(padj = p.adjust(P.value, method = padj_method))
}

############################################
## 2. Disease gene sets: DNV vs Candidate
############################################

disease_gene_dict <- list(
  DNV       = fu72_id,
  Candidate = c(eur_candidate %>% output_gene_id(),
                eas_candidate %>% output_gene_id())
)

## background gene set
bg <- cdsGenes

############################################
## 3. FET for disorder gene lists
############################################
disorder_tbl   <- read.delim("path/to/disorder_gene_list.txt")
risk_gene_list <- split(disorder_tbl$gene_id, disorder_tbl$Disorder)

fet_disorder_results <- run_fet_sets(
  gene_set_list      = risk_gene_list,
  disease_gene_dict  = disease_gene_dict,
  bg                 = bg
)

head(fet_disorder_results)

## Save disorder FET results
write.table(
  fet_disorder_results,
  file      = "path/to/FET_res_Disorder_Geneset.tsv",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

############################################
## 4. FET for cell-type DEGs
############################################
cell_deg_tbl  <- read.delim("path/to/celltype_DEG_gene_list.txt")
cell_deg_list <- split(cell_deg_tbl$gene_id, cell_deg_tbl$Cluster)

fet_celltype_results <- run_fet_sets(
  gene_set_list      = cell_deg_list,
  disease_gene_dict  = disease_gene_dict,
  bg                 = bg
)

head(fet_celltype_results)

## Save cell-type DEG FET results
write.table(
  fet_celltype_results,
  file      = "path/to/FET_res_CellType_DEGs.tsv",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)