#########################################################
## RareComb pipeline: multi-cohort pathogenic variant pairs
##  - Per-cohort pathogenic filtering (no rare/MAF filtering)
##  - Cohort merging
##  - RareComb boolean input matrix
##  - Running pyRareComb (example commands)
##  - Downstream analysis of RareComb pairwise results
#########################################################

library(data.table)
library(dplyr)
library(tidyr)
library(purrr)

## user-defined helper functions
source("path/to/function_definition.R")

#########################################################
## 0. Paths / configuration
#########################################################

## Cohort-level variant tables:
SSC_file <- "/path/to/SSC_pathogenic_input.rds"
SPARK_WGS_file <- "/path/to/SPARK_WGS_pathogenic_input.rds"
SPARK_WES_file <- "/path/to/SPARK_WES_pathogenic_input.rds"
Korean_WGS_file <- "/path/to/Korean_WGS_pathogenic_input.rds"
Korean_WES_file <- "/path/to/Korean_WES_pathogenic_input.rds"

## Sample metadata containing merged information across all cohorts
sample_info_file <- "/path/to/EAS_sample_list.txt"

## Working directory
work_dir     <- "/path/to/project_root"
rarecomb_dir <- file.path(work_dir, "rarecomb")
dir.create(rarecomb_dir, showWarnings = FALSE, recursive = TRUE)

## RareComb input/output
rarecomb_input_csv <- file.path(
  rarecomb_dir,
  "EAS_AllCohorts.PathogenicOnly.rarecomb_input_matrix.csv.gz"
)
rarecomb_prefix     <- "EAS_AllCohorts_PathogenicOnly"
rarecomb_output_dir <- rarecomb_dir

## Expected RareComb output files
rarecomb_pair_file <- file.path(
  rarecomb_dir,
  sprintf("length2_min5_max0.25_fpgrowth.%s.csv", rarecomb_prefix)
)
rarecomb_triple_file <- file.path(
  rarecomb_dir,
  sprintf("length3_min5_max0.25_fpgrowth.%s.csv", rarecomb_prefix)
)

#########################################################
## 1. Pathogenic missense / PTV filter (no rare/MAF filter)
#########################################################

filter_pathogenic_variants <- function(dat,
                                       cadd_cutoff = 20,
                                       mis_mean_cutoff = 0.7) {
  dat %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      mis_mean_score = mean(
        revel_rankscore,
        metarnn_rankscore,
        primateai_rankscore,
        fathmm_converted_rankscore,
        mpc_rankscore,
        provean_converted_rankscore,
        vest4_rankscore,
        bayesdel_addaf_rankscore,
        gmvp_rankscore,
        na.rm = TRUE
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(
      VariantType == "protein_truncating" |
        (VariantType == "missense" &
           cadd_phred >= cadd_cutoff &
           mis_mean_score >= mis_mean_cutoff)
    )
}

#########################################################
## 2. Load and filter each cohort
#########################################################
SSC_dat <- readRDS(SSC_file) %>%
  filter_pathogenic_variants() %>%
  dplyr::mutate(cohort = ifelse("cohort" %in% names(.), cohort, "SSC"))

SPARK_WGS_dat <- readRDS(SPARK_WGS_file) %>%
  filter_pathogenic_variants() %>%
  dplyr::mutate(cohort = ifelse("cohort" %in% names(.), cohort, "SPARK_WGS"))

SPARK_WES_dat <- readRDS(SPARK_WES_file) %>%
  filter_pathogenic_variants() %>%
  dplyr::mutate(cohort = ifelse("cohort" %in% names(.), cohort, "SPARK_WES"))
  
Korean_WGS_dat <- readRDS(Korean_WGS_file) %>%
  filter_pathogenic_variants() %>%
  dplyr::mutate(cohort = ifelse("cohort" %in% names(.), cohort, "Korean_WGS"))

Korean_WES_dat <- readRDS(Korean_WES_file) %>%
  filter_pathogenic_variants() %>%
  dplyr::mutate(cohort = ifelse("cohort" %in% names(.), cohort, "Korean_WES"))

#########################################################
## 3. Merge cohorts
#########################################################

setwd(work_dir)

variant_all <- dplyr::bind_rows(SSC_dat, SPARK_WGS_dat, SPARK_WES_dat, Korean_WGS_dat, Korean_WES_dat)

saveRDS(
  variant_all,
  file = file.path(
    rarecomb_dir,
    "EAS_AllCohorts.PathogenicOnly.long_format.rds"
  )
)

#########################################################
## 4. Generate RareComb input matrix (boolean gene × sample)
#########################################################
rarecomb_input <- variant_all
mt_bool <- make_boolean_mt(rarecomb_input)

fwrite(
  mt_bool,
  file      = rarecomb_input_csv,
  row.names = FALSE,
  quote     = FALSE,
  col.names = TRUE,
  sep       = ","
)

#########################################################
## 5. Running pyRareComb (example command)
#########################################################

system(
  "python run_PyRarecomb.py \
   EAS_AllCohorts.PathogenicOnly.rarecomb_input_matrix.csv.gz \
   2 5 0.25 fpgrowth 0.2 EAS_AllCohorts_PathogenicOnly rarecomb_output/"
)

#########################################################
## 6. Downstream analysis of RareComb pairwise results
#########################################################

## 6-1. Load results
dat_input   <- variant_all
dat_info    <- read.delim(sample_info_file)
comb_pairs  <- read.csv(rarecomb_pair_file)
comb_triple <- read.csv(rarecomb_triple_file)

## number of cases / controls
n_case <- dat_input %>%
  dplyr::select(s, Group) %>%
  dplyr::distinct() %>%
  dplyr::filter(Group == "Case") %>%
  nrow()

n_ctrl <- dat_input %>%
  dplyr::select(s, Group) %>%
  dplyr::distinct() %>%
  dplyr::filter(Group == "Control") %>%
  nrow()

## 6-2. Case/control filters + binomial test + relative risk
comb_pairs_filtered <- comb_pairs %>%
  dplyr::filter(Case_Adj_Pval_BH < 0.2, Cont_pvalue_more > 0.2) %>%
  annot_symbol(., k = 2) %>%
  dplyr::mutate(itemset_id = paste0(Item_1, "/", Item_2)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    binom_p = binom_test(
      x  = Case_Obs_Count_Combo,
      n  = Case_Obs_Count_Combo + Cont_Obs_Count_Combo,
      p  = n_case / (n_case + n_ctrl),
      alternative = "greater"
    )$p,
    RR = Case_Obs_Prob_Combo / Cont_Obs_Prob_Combo
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    Case_Adj_Pval_BH < 0.2,
    Cont_pvalue_more > 0.2,
    binom_p < 0.05
  )

## 6-3. Inheritance patterns in cases
pair_inheritance_case <- inheritance_pattern(
  x = comb_pairs_filtered,
  d = dat_input,
  parallel = TRUE
) %>%
  dplyr::distinct() %>%
  dplyr::left_join(dat_info, by = c("s" = "sample_id")) %>%
  dplyr::distinct()

## 6-4. Summaries of inheritance distribution by gene pair
mat_family_ids <- dat_info %>%
  dplyr::filter(Group == "Case", ROLE == "parent", sex == "Female") %>%
  dplyr::pull(family_id)

pat_family_ids <- dat_info %>%
  dplyr::filter(Group == "Case", ROLE == "parent", sex == "Male") %>%
  dplyr::pull(family_id)

pair_inherit_summary <- pair_inheritance_case %>%
  dplyr::select(family_id, itemset_id, Inherit_pattern, s, ROLE) %>%
  dplyr::group_by(itemset_id) %>%
  dplyr::summarise(
    nOccur = dplyr::n_distinct(s),
    nFam   = dplyr::n_distinct(family_id),
    nBipar = dplyr::n_distinct(s[Inherit_pattern %in% c(
      "Maternal + Paternal",
      "DNV + Maternal",
      "DNV + Paternal"
    )]),
    nMat   = dplyr::n_distinct(s[Inherit_pattern == "Both Maternal"]),
    nPat   = dplyr::n_distinct(s[Inherit_pattern == "Both Paternal"]),
    nCasePatParent = dplyr::n_distinct(s[
      Inherit_pattern == "Both Paternal" &
        family_id %in% pat_family_ids
    ]),
    nCaseMatParent = dplyr::n_distinct(s[
      Inherit_pattern == "Both Maternal" &
        family_id %in% mat_family_ids
    ]),
    nCaseParentUnphased = dplyr::n_distinct(s[
      Inherit_pattern == "unknown phase" &
        ROLE == "parent"
    ])
  ) %>%
  dplyr::mutate(
    p_Bipar = (nBipar + nCaseParentUnphased + nCasePatParent + nCaseMatParent) / nOccur,
    p_Mat   = nMat / nOccur,
    p_Pat   = nPat / nOccur,
    p_MorP  = (nMat + nPat) / nOccur
  )

pair_annotated <- comb_pairs_filtered %>%
  dplyr::left_join(pair_inherit_summary, by = "itemset_id")

## 6-5. Candidate gene pairs (example: biparental fraction ≥ 0.6)
pair_candidates <- pair_annotated %>%
  dplyr::filter(p_Bipar >= 0.6)

pair_candidates_case <- pair_inheritance_case %>%
  dplyr::filter(itemset_id %in% pair_candidates$itemset_id) %>%
  dplyr::distinct()

## 6-6. Inheritance patterns in controls
pair_inheritance_ctrl <- inheritance_pattern_cntl(
  x = comb_pairs_filtered,
  d = dat_input,
  parallel = TRUE
) %>%
  dplyr::left_join(dat_info, by = c("s" = "sample_id")) %>%
  dplyr::filter(!is.na(s)) %>%
  dplyr::distinct()

pair_candidates_ctrl <- pair_inheritance_ctrl %>%
  dplyr::filter(itemset_id %in% pair_candidates$itemset_id) %>%
  dplyr::distinct()

#########################################################
## 7. Save output objects
#########################################################

save(
  pair_candidates,
  pair_candidates_case,
  pair_inheritance_ctrl,
  pair_candidates_ctrl,
  comb_pairs_filtered,
  pair_inheritance_case,
  pair_annotated,
  dat_input,
  file = file.path(
    rarecomb_dir,
    "eas_candidate_pairs.rds"
  )
)
