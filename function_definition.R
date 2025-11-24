## RareComb helpers (minimal set used in the pipeline)

library(dplyr)
library(tidyr)
library(tibble)
library(biomaRt)
library(parallel)
library(doParallel)
library(foreach)

## number of cores for parallel sections
if (!exists("numCores")) {
  nc <- parallel::detectCores()
  numCores <- max(1L, nc - 1L)
}

## Ensembl mart (used for gene annotation)
if (!exists("mart")) {
  mart <- biomaRt::useMart(
    "ENSEMBL_MART_ENSEMBL",
    dataset = "hsapiens_gene_ensembl"
  )
}

## ------------------------------------------------------------------
## 1. Boolean matrix for RareComb (samples x genes, 0/1)
## ------------------------------------------------------------------

## d: long-format variant table with at least:
##    - s        : sample ID
##    - gene_id  : Ensembl gene ID
##    - Group    : "Case" / "Control"
make_boolean_mt <- function(d) {
  gene_mat <- d %>%
    distinct(s, gene_id, .keep_all = TRUE) %>%
    mutate(n = 1L) %>%
    pivot_wider(
      names_from  = gene_id,
      values_from = n,
      values_fill = list(n = 0L),
      names_prefix = "Input_"
    )
  
  grp <- d %>%
    distinct(s, Group) %>%
    mutate(Output_group = as.integer(Group == "Case"))
  
  gene_mat %>%
    left_join(grp, by = "s") %>%
    rename(Sample_Name = s)
}

## ------------------------------------------------------------------
## 2. Small utilities
## ------------------------------------------------------------------

output_gene_id <- function(output) {
  c(output$Item_1, output$Item_2, output$Item_3, output$Item_4, output$Item_5) %>% unique()
}

output_genes <- function(output){
  c(output$gene_symbol_1, output$gene_symbol_2, output$gene_symbol_3, output$gene_symbol_4, output$gene_symbol_5) %>% unique()
}


## ------------------------------------------------------------------
## 3. Gene symbol annotation for RareComb itemsets
## ------------------------------------------------------------------

## output: RareComb result table
## k     : item length (2, 3, 4, 5)
annot_symbol <- function(output, k = 2L) {
  if (k == 2L) {
    output <- output %>%
      separate(Item_1, c(NA, "Item_1")) %>%
      separate(Item_2, c(NA, "Item_2"))
  } else if (k == 3L) {
    output <- output %>%
      separate(Item_1, c(NA, "Item_1")) %>%
      separate(Item_2, c(NA, "Item_2")) %>%
      separate(Item_3, c(NA, "Item_3"))
  } else if (k == 4L) {
    output <- output %>%
      separate(Item_1, c(NA, "Item_1")) %>%
      separate(Item_2, c(NA, "Item_2")) %>%
      separate(Item_3, c(NA, "Item_3")) %>%
      separate(Item_4, c(NA, "Item_4"))
  } else if (k == 5L) {
    output <- output %>%
      separate(Item_1, c(NA, "Item_1")) %>%
      separate(Item_2, c(NA, "Item_2")) %>%
      separate(Item_3, c(NA, "Item_3")) %>%
      separate(Item_4, c(NA, "Item_4")) %>%
      separate(Item_5, c(NA, "Item_5"))
  }
  
  bm <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"),
    filters    = "ensembl_gene_id",
    values     = output_gene_id(output),
    mart       = mart
  )
  
  if (k >= 1L) {
    output <- merge(
      bm, output,
      by.x = "ensembl_gene_id", by.y = "Item_1", all.y = TRUE
    ) %>%
      rename(
        Item_1        = ensembl_gene_id,
        gene_symbol_1 = hgnc_symbol,
        Item_1_chr    = chromosome_name
      )
  }
  if (k >= 2L) {
    output <- merge(
      bm, output,
      by.x = "ensembl_gene_id", by.y = "Item_2", all.y = TRUE
    ) %>%
      rename(
        Item_2        = ensembl_gene_id,
        gene_symbol_2 = hgnc_symbol,
        Item_2_chr    = chromosome_name
      )
  }
  if (k >= 3L) {
    output <- merge(
      bm, output,
      by.x = "ensembl_gene_id", by.y = "Item_3", all.y = TRUE
    ) %>%
      rename(
        Item_3        = ensembl_gene_id,
        gene_symbol_3 = hgnc_symbol,
        Item_3_chr    = chromosome_name
      )
  }
  if (k >= 4L) {
    output <- merge(
      bm, output,
      by.x = "ensembl_gene_id", by.y = "Item_4", all.y = TRUE
    ) %>%
      rename(
        Item_4        = ensembl_gene_id,
        gene_symbol_4 = hgnc_symbol,
        Item_4_chr    = chromosome_name
      )
  }
  if (k >= 5L) {
    output <- merge(
      bm, output,
      by.x = "ensembl_gene_id", by.y = "Item_5", all.y = TRUE
    ) %>%
      rename(
        Item_5        = ensembl_gene_id,
        gene_symbol_5 = hgnc_symbol,
        Item_5_chr    = chromosome_name
      )
  }
  
  if (k == 2L) {
    output$itemset_id <- paste(output$Item_1, output$Item_2, sep = "/")
  } else if (k == 3L) {
    output$itemset_id <- paste(output$Item_1, output$Item_2, output$Item_3, sep = "/")
  } else if (k == 4L) {
    output$itemset_id <- paste(output$Item_1, output$Item_2, output$Item_3, output$Item_4, sep = "/")
  } else if (k == 5L) {
    output$itemset_id <- paste(output$Item_1, output$Item_2, output$Item_3, output$Item_4, output$Item_5, sep = "/")
  }
  
  output
}

## ------------------------------------------------------------------
## 4. Inheritance pattern (cases)
## ------------------------------------------------------------------

## x: RareComb result (Item_1, Item_2)
## d: long table with (s, gene_id, Group, transmission_pattern)
inheritance_pattern <- function(x, d, parallel = TRUE, verbose = FALSE) {
  if (verbose) message("Computing inheritance patterns in cases...")
  
  if (parallel) {
    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    
    a <- foreach(
      i = 1:nrow(x),
      .combine  = rbind,
      .packages = c("dplyr", "tibble")
    ) %dopar% {
      item_1 <- x[i, "Item_1"] %>% unlist()
      item_2 <- x[i, "Item_2"] %>% unlist()
      
      s <- intersect(
        d %>% filter(Group == "Case", gene_id == item_2) %>% pull(s),
        d %>% filter(Group == "Case", gene_id == item_1) %>% pull(s)
      )
      tibble(Item_2 = item_2, Item_1 = item_1, s = s)
    }
    
    stopCluster(cl)
  } else {
    a <- tibble()
    for (i in seq_len(nrow(x))) {
      item_1 <- x[i, "Item_1"] %>% unlist()
      item_2 <- x[i, "Item_2"] %>% unlist()
      
      s <- intersect(
        d %>% filter(Group == "Case", gene_id == item_2) %>% pull(s),
        d %>% filter(Group == "Case", gene_id == item_1) %>% pull(s)
      )
      a <- bind_rows(a, tibble(Item_2 = item_2, Item_1 = item_1, s = s))
    }
  }
  
  a <- merge(
    a,
    d %>% select(gene_id, s, Item_2_from = transmission_pattern),
    by.x = c("Item_2", "s"), by.y = c("gene_id", "s"), all.x = TRUE
  )
  a <- merge(
    a,
    d %>% select(gene_id, s, Item_1_from = transmission_pattern),
    by.x = c("Item_1", "s"), by.y = c("gene_id", "s"), all.x = TRUE
  ) %>%
    unique()
  
  ## collapse duplicate transmission info per (Item_1, Item_2, s)
  if (a[duplicated(a[, 1:4], fromLast = TRUE), ]$Item_1_from %>% length() != 0) {
    a[duplicated(a[, 1:4], fromLast = TRUE), ]$Item_1_from <-
      paste0(
        a[duplicated(a[, 1:4], fromLast = TRUE), ]$Item_1_from, ", ",
        a[duplicated(a[, 1:4], fromLast = FALSE), ]$Item_1_from
      )
  }
  if (a[duplicated(a[, c(1:3, 5)], fromLast = TRUE), ]$Item_2_from %>% length() != 0) {
    a[duplicated(a[, c(1:3, 5)], fromLast = TRUE), ]$Item_2_from <-
      paste0(
        a[duplicated(a[, c(1:3, 5)], fromLast = TRUE), ]$Item_2_from, ", ",
        a[duplicated(a[, c(1:3, 5)], fromLast = FALSE), ]$Item_2_from
      )
  }
  a <- a[!duplicated(a[, 1:3], fromLast = FALSE), ]
  
  a <- a %>%
    mutate(
      Inherit_pattern = case_when(
        grepl("father", paste(Item_2_from, Item_1_from)) &
          grepl("mother", paste(Item_2_from, Item_1_from)) &
          !grepl("DNV", paste(Item_2_from, Item_1_from)) &
          !grepl("unkown_phase", paste(Item_2_from, Item_1_from)) ~ "Maternal + Paternal",
        grepl("father", paste(Item_2_from, Item_1_from)) &
          !grepl("mother", paste(Item_2_from, Item_1_from)) &
          !grepl("DNV", paste(Item_2_from, Item_1_from)) &
          !grepl("unkown_phase", paste(Item_2_from, Item_1_from)) ~ "Both Paternal",
        !grepl("father", paste(Item_2_from, Item_1_from)) &
          grepl("mother", paste(Item_2_from, Item_1_from)) &
          !grepl("DNV", paste(Item_2_from, Item_1_from)) &
          !grepl("unkown_phase", paste(Item_2_from, Item_1_from)) ~ "Both Maternal",
        !grepl("father", paste(Item_2_from, Item_1_from)) &
          !grepl("mother", paste(Item_2_from, Item_1_from)) &
          grepl("DNV", paste(Item_2_from, Item_1_from)) &
          !grepl("unkown_phase", paste(Item_2_from, Item_1_from)) ~ "Both DNV",
        grepl("father", paste(Item_2_from, Item_1_from)) &
          !grepl("mother", paste(Item_2_from, Item_1_from)) &
          grepl("DNV", paste(Item_2_from, Item_1_from)) &
          !grepl("unkown_phase", paste(Item_2_from, Item_1_from)) ~ "DNV + Paternal",
        !grepl("father", paste(Item_2_from, Item_1_from)) &
          grepl("mother", paste(Item_2_from, Item_1_from)) &
          grepl("DNV", paste(Item_2_from, Item_1_from)) &
          !grepl("unkown_phase", paste(Item_2_from, Item_1_from)) ~ "DNV + Maternal",
        grepl("father", paste(Item_2_from, Item_1_from)) &
          grepl("mother", paste(Item_2_from, Item_1_from)) &
          grepl("DNV", paste(Item_2_from, Item_1_from)) &
          !grepl("unkown_phase", paste(Item_2_from, Item_1_from)) ~ "Maternal + Paternal + DNV",
        TRUE ~ "unkown phase"
      )
    )
  
  merge(x, a, by = c("Item_2", "Item_1"), all.x = TRUE) %>%
    unique()
}

## ------------------------------------------------------------------
## 5. Inheritance pattern (controls)
## ------------------------------------------------------------------

inheritance_pattern_cntl <- function(x, d, parallel = TRUE, verbose = FALSE) {
  if (verbose) message("Computing inheritance patterns in controls...")
  
  if (parallel) {
    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    
    a <- foreach(
      i = 1:nrow(x),
      .combine  = rbind,
      .packages = c("dplyr", "tibble")
    ) %dopar% {
      item_1 <- x[i, "Item_1"] %>% unlist()
      item_2 <- x[i, "Item_2"] %>% unlist()
      
      s <- intersect(
        d %>% filter(Group == "Control", gene_id == item_2) %>% pull(s),
        d %>% filter(Group == "Control", gene_id == item_1) %>% pull(s)
      )
      tibble(Item_2 = item_2, Item_1 = item_1, s = s)
    }
    
    stopCluster(cl)
  } else {
    a <- tibble()
    for (i in seq_len(nrow(x))) {
      item_1 <- x[i, "Item_1"] %>% unlist()
      item_2 <- x[i, "Item_2"] %>% unlist()
      
      s <- intersect(
        d %>% filter(Group == "Control", gene_id == item_2) %>% pull(s),
        d %>% filter(Group == "Control", gene_id == item_1) %>% pull(s)
      )
      a <- bind_rows(a, tibble(Item_2 = item_2, Item_1 = item_1, s = s))
    }
  }
  
  a <- merge(
    a,
    d %>% select(gene_id, s, Item_2_from = transmission_pattern),
    by.x = c("Item_2", "s"), by.y = c("gene_id", "s"), all.x = TRUE
  )
  a <- merge(
    a,
    d %>% select(gene_id, s, Item_1_from = transmission_pattern),
    by.x = c("Item_1", "s"), by.y = c("gene_id", "s"), all.x = TRUE
  ) %>%
    unique()
  
  if (a[duplicated(a[, 1:4], fromLast = TRUE), ]$Item_1_from %>% length() != 0) {
    a[duplicated(a[, 1:4], fromLast = TRUE), ]$Item_1_from <-
      paste0(
        a[duplicated(a[, 1:4], fromLast = TRUE), ]$Item_1_from, ", ",
        a[duplicated(a[, 1:4], fromLast = FALSE), ]$Item_1_from
      )
  }
  if (a[duplicated(a[, c(1:3, 5)], fromLast = TRUE), ]$Item_2_from %>% length() != 0) {
    a[duplicated(a[, c(1:3, 5)], fromLast = TRUE), ]$Item_2_from <-
      paste0(
        a[duplicated(a[, c(1:3, 5)], fromLast = TRUE), ]$Item_2_from, ", ",
        a[duplicated(a[, c(1:3, 5)], fromLast = FALSE), ]$Item_2_from
      )
  }
  a <- a[!duplicated(a[, 1:3], fromLast = FALSE), ]
  
  a <- a %>%
    mutate(
      Inherit_pattern = case_when(
        grepl("father", paste(Item_2_from, Item_1_from)) &
          grepl("mother", paste(Item_2_from, Item_1_from)) &
          !grepl("DNV", paste(Item_2_from, Item_1_from)) &
          !grepl("unkown_phase", paste(Item_2_from, Item_1_from)) ~ "Maternal + Paternal",
        grepl("father", paste(Item_2_from, Item_1_from)) &
          !grepl("mother", paste(Item_2_from, Item_1_from)) &
          !grepl("DNV", paste(Item_2_from, Item_1_from)) &
          !grepl("unkown_phase", paste(Item_2_from, Item_1_from)) ~ "Both Paternal",
        !grepl("father", paste(Item_2_from, Item_1_from)) &
          grepl("mother", paste(Item_2_from, Item_1_from)) &
          !grepl("DNV", paste(Item_2_from, Item_1_from)) &
          !grepl("unkown_phase", paste(Item_2_from, Item_1_from)) ~ "Both Maternal",
        !grepl("father", paste(Item_2_from, Item_1_from)) &
          !grepl("mother", paste(Item_2_from, Item_1_from)) &
          grepl("DNV", paste(Item_2_from, Item_1_from)) &
          !grepl("unkown_phase", paste(Item_2_from, Item_1_from)) ~ "Both DNV",
        grepl("father", paste(Item_2_from, Item_1_from)) &
          !grepl("mother", paste(Item_2_from, Item_1_from)) &
          grepl("DNV", paste(Item_2_from, Item_1_from)) &
          !grepl("unkown_phase", paste(Item_2_from, Item_1_from)) ~ "DNV + Paternal",
        !grepl("father", paste(Item_2_from, Item_1_from)) &
          grepl("mother", paste(Item_2_from, Item_1_from)) &
          grepl("DNV", paste(Item_2_from, Item_1_from)) &
          !grepl("unkown_phase", paste(Item_2_from, Item_1_from)) ~ "DNV + Maternal",
        grepl("father", paste(Item_2_from, Item_1_from)) &
          grepl("mother", paste(Item_2_from, Item_1_from)) &
          grepl("DNV", paste(Item_2_from, Item_1_from)) &
          !grepl("unkown_phase", paste(Item_2_from, Item_1_from)) ~ "Maternal + Paternal + DNV",
        is.na(Item_2_from) | is.na(Item_1_from) ~ "unknown phase",
        TRUE ~ "Unknown phase"
      )
    )
  
  merge(x, a, by = c("Item_2", "Item_1"), all.x = TRUE) %>%
    unique()
}

## ------------------------------------------------------------------
## 6. Binomial test wrapper
## ------------------------------------------------------------------

binom_test <- function(x, n, p, alternative = "greater") {
  res <- stats::binom.test(x = x, n = n, p = p, alternative = alternative)
  list(p = res$p.value)
}


## ------------------------------------------------------------------
## 7. Plotting
## ------------------------------------------------------------------
theme_v1 <- theme(
  text=element_text(
    size = 10),
  legend.background = element_blank(),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.key=element_blank(),
  axis.text.x = element_text(size = 10, color = 'black'),
  axis.text.y = element_text(size = 10, color = 'black'),
  axis.title.y = element_text(size=10),
  axis.title.x = element_text(size=10),
  plot.title = element_text(size = 12, vjust = 1.5),
  panel.grid.major = element_line(color = "black", size = .1, linetype = 3),
  panel.grid.minor = element_line(color = "black", size = .1, linetype = 3),
  panel.background = element_blank(),
  panel.border = element_blank(),
  axis.line.y = element_line(colour = "black", size = .35),
  axis.line.x.bottom = element_line(colour = "black", size = .35),
  axis.ticks = element_blank(),
  strip.background = element_rect(colour = 'black'
                                  , fill = 'gray88'
                                  , linetype = 0,
                                  size = .3),
  strip.text = element_text(size=10),
  legend.key.size = unit(0.38, 'cm'),
  legend.margin=margin(t = .1, r = .1, b = .1, l = .1, unit = 'cm'))



plot_inherit_pattern <- function(x){
  x %>% count(Inherit_pattern) %>%
    # mutate(Inherit_pattern = factor(Inherit_pattern, levels = c('Maternal + Paternal', 'Both Paternal', 'Both Maternal', 'DNV + Paternal', 'DNV + Maternal', 'Unknown phase'))) %>%
    mutate(percent = paste0(prop.table(n) %>% round(., digits = 2)*100) %>% as.integer()) %>% 
    arrange(percent) %>%
    ggplot(aes(reorder(Inherit_pattern, -percent), percent, fill=Inherit_pattern)) +
    geom_bar(stat='identity') + 
    geom_text(aes(label=n), vjust = -0.4, size = 2.5) +
    xlab('Inheritance pattern') + 
    ylab('Proportion ') +
    scale_x_discrete(breaks = NULL) +
    scale_fill_manual(values = c("#4063a3", "#34b6c6", "#79ad41", "#ddc000", "#d7aca1", 'grey')) +
    theme(legend.title = element_blank()) +
    theme_v1 +
    guides(fill=guide_legend(title="")) 
}



## ------------------------------------------------------------------
## 8. Required data
## ------------------------------------------------------------------
# Genecode genes
gm = readRDS(file = 'gencode.v44.annotation.rds')
cdsGenes = gm %>% filter(gene_type == 'protein_coding') %>% pull(gene_id)

# ASD-related genes
fu2022_tada_file <- "path/to/Fu2022_TADA_geneset.xlsx"
fu = read_excel(fu2022_tada_file)
fu72 = fu %>% filter(ASD72 == 1) %>% pull(gene_gencodeV33)
fu72_id = fu %>% filter(ASD72 == 1) %>% pull(gene_id)

sfari_file       <- "path/to/sfari_geneset.csv"
sfari_dat <- read.csv(sfari_file) %>%
  dplyr::filter(gene.score %in% c(1, 2))
sfari_dat$ensembl.id[sfari_dat$gene.symbol == "ADA"]  <- "ENSG00000105976"
sfari_dat$ensembl.id[sfari_dat$gene.symbol == "MET"]  <- "ENSG00000196839"
sfari_dat$ensembl.id[sfari_dat$gene.symbol == "PHB1"] <- "ENSG00000167085"
sfari_dat$ensembl.id[sfari_dat$gene.symbol == "AR"]   <- "ENSG00000169083"
sf = sf %>% filter(ensembl.id %in% cdsGenes)
sf_id = sf$`ensembl.id`