#########################################################
# Whole-genome variant processing pipeline
#########################################################

import hail as hl
import time
import os
from functools import reduce


def quality_control_wgs(vcf_path, output_dir, chr_name, lcr_bed_path):
    start = time.time()
    
    # 1) Import raw VCF as MatrixTable
    mt = hl.import_vcf(vcf_path, force_bgz=True, reference_genome='GRCh38', skip_invalid_loci=True)

    # 2) Filter out low-complexity regions (LCR) and very long alleles (length ≥ 50bp)
    lcr_bed = hl.import_bed(lcr_bed_path, reference_genome='GRCh38', min_partitions=100)
    mt = mt.filter_rows(~hl.is_defined(lcr_bed[mt.locus]))
    mt = mt.filter_rows((mt.alleles[0].length() < 50) & (mt.alleles[1].length() < 50))
    
    # 3) Recompute total depth (DP) from AD  - For only Korean ASD WGS
    mt = mt.annotate_entries(DP = hl.sum(mt.AD))

    # 4) Normalize local likelihood fields (AD, LPL, LA) for multi-allelic handling - For only Korean ASD WGS
    mt = mt.annotate_entries(
        AD  = hl.case()
                .when(mt.GT.is_hom_ref(), hl.vds.local_to_global(hl.array([mt.DP]), hl.array([0]), n_alleles=hl.len(mt.alleles), fill_value=0, number='R'))
                .default(mt.AD),
        LPL = hl.case()
                .when(mt.GT.is_hom_ref(), hl.missing('array<int32>'))
                .default(mt.LPL),
        LA  = hl.case()
                .when(mt.GT.is_hom_ref(),      hl.array([0]))
                .when(mt.GT.is_het_ref(),      hl.array([hl.int32(0), hl.int32(mt.LAA[0])]))
                .when(mt.GT.is_hom_var(),      hl.array([hl.int32(0), hl.int32(mt.LAA[0])]))
                .when(mt.GT.is_het_non_ref(),  hl.array([hl.int32(0), hl.int32(mt.LAA[0]), hl.int32(mt.LAA[1])]))
                .or_missing()
    )

    # 5) Convert local PL (LPL/LA) to global diploid PLs for all alleles  - For only Korean ASD WGS
    mt = mt.annotate_entries(PL = hl.vds.local_to_global(mt.LPL, mt.LA, n_alleles=hl.len(mt.alleles), fill_value=hl.max(mt.LPL), number='G'))

    # 6) Split multi-allelic sites and left-normalize to minimal representation
    mt = hl.split_multi_hts(mt, left_aligned=True)
    mt = mt.annotate_rows(min_rep = hl.min_rep(mt.locus, mt.alleles)).key_rows_by().drop('locus', 'alleles')
    mt = mt.annotate_rows(locus = mt.min_rep.locus, alleles = mt.min_rep.alleles).key_rows_by('locus', 'alleles')

    # 7) Re-define PL for hom-ref genotypes  - For only Korean ASD WGS
    mt = mt.annotate_entries(
        PL = hl.case()
                .when(mt.GT.is_hom_ref(), [0, mt.GQ, hl.max(mt.GQ, 3*(mt.AD[0] - mt.AD[1]))])
                .default(mt.PL)
    )

    # 8) Annotate variant class (SNP / indel) and compute variant-level QC metrics
    mt = mt.annotate_rows(variant_class = hl.case().when(hl.is_snp(mt.alleles[0], mt.alleles[1]), 'snp').when(hl.is_indel(mt.alleles[0], mt.alleles[1]), 'indel').or_missing())
    mt = hl.variant_qc(mt)

    # 9) Compute allelic balance (AB) from AD
    mt = mt.annotate_entries(AB = mt.AD[1] / hl.sum(mt.AD))

    # 10) Genotype-level QC filters (DP/GQ thresholds and AB-based filters for het/hom-var)
    mt = mt.filter_entries(
        (mt.GT.is_het()     & (0.2 <= mt.AB) & (mt.AB <= 0.8)  & (mt.DP >= 10) & (mt.GQ >= 20)) |
        (mt.GT.is_hom_ref() &                                 (mt.DP >= 10) & (mt.GQ >= 20))     |
        (mt.GT.is_hom_var() & (mt.AB >= 0.95)                & (mt.DP >= 10) & (mt.GQ >= 20))
    )

    # 11) Run variant_qc after genotype filtering and apply variant-level QC filters
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows((mt.variant_qc.call_rate >= 0.1) & (mt.variant_qc.p_value_hwe > 1e-12))

    # 12) Write QCed MatrixTable to disk
    mt.write(os.path.join(output_dir, f'WGS_chr{chr_name}_QC.mt'), overwrite=True)
    print(f'Done QC for chr{chr_name} : Execution time = {time.time() - start} sec')
    return mt




def call_denovo_wgs(mt, ped_path, gnomad_ht_path, output_dir, chr_name):
    start = time.time()
    
    # 1) Read external allele frequency resource (gnomAD) and annotate AF to rows
    gnomad_ht = hl.read_table(gnomad_ht_path)
    mt = mt.annotate_rows(gnomAD_AF = gnomad_ht[mt.locus, mt.alleles].AF)

    # 2) Restrict to very rare (or absent) variants based on gnomAD AF (< 0.01%)
    mt_dnv = mt.filter_rows((mt.gnomAD_AF < 0.0001) | hl.is_missing(mt.gnomAD_AF))

    # 3) Load pedigree information and construct trio relationships
    trio_fam = hl.Pedigree.read(ped_path)

    # 4) Run Hail's trio-based de novo caller with depth/AB/GQ-based thresholds
    # dnv_ht = modified_de_novo( #  - For only Korean ASD WGS
    dnv_ht = hl.de_novo(
        mt                 = mt_dnv,
        pedigree           = trio_fam,
        pop_frequency_prior= mt.gnomAD_AF,
        max_parent_ab      = 0.1,
        min_child_ab       = 0.3,
        min_dp_ratio       = 0.3,
        min_gq             = 20
        #min_p              = 0.05 #  - For only Korean ASD WGS
    )

    # 5) Keep only HIGH and MEDIUM confidence de novo candidates
    dnv_ht = dnv_ht.filter((dnv_ht.confidence == "HIGH") | (dnv_ht.confidence == "MEDIUM"))

    # 6) Write de novo results as Hail Table
    dnv_ht.write(os.path.join(output_dir, f'WGS_chr{chr_name}_DNV.ht'), overwrite=True)
    print(f'Done DNV calling for chr{chr_name} : Execution time = {time.time() - start} sec')
    return dnv_ht




def process_rare_inherited_wgs(mt, gnomad_ht_path, ped_path, sample_info_path, internal_ac_ht_path, output_dir, chr_name):
    start = time.time()

    # 1) Read external AF resource (gnomAD) and annotate allele frequency to rows
    gnomad_ht = hl.read_table(gnomad_ht_path)
    mt = mt.annotate_rows(gnomAD_AF = gnomad_ht[mt.locus, mt.alleles].AF)

    # 2) Genotype-level QC for inherited variants (AB, PL, DP, GQ filters by genotype class)
    mt = mt.filter_entries(
        hl.case()
          .when(mt.GT.is_het(),     (0.22 <= mt.AB) & (mt.AB <= 0.78) & (mt.PL[0] >= 25) & (mt.DP >= 10) & (mt.GQ >= 20))
          .when(mt.GT.is_hom_var(), (mt.AB >= 0.95)                    & (mt.PL[0] >= 25) & (mt.DP >= 10) & (mt.GQ >= 20))
          .default(                                             (mt.DP >= 10) & (mt.GQ >= 20))
    )

    # 3) Variant-level QC (call rate and HWE)
    mt = mt.filter_rows((mt.variant_qc.call_rate >= 0.1) & (mt.variant_qc.p_value_hwe > 1e-12))

    # 4) Define rare and ultra-rare variant subsets based on external AF
    mt_rare       = mt.filter_rows((mt.gnomAD_AF < 0.001)   | hl.is_missing(mt.gnomAD_AF))
    mt_ultra_rare = mt.filter_rows((mt.gnomAD_AF < 0.00001) | hl.is_missing(mt.gnomAD_AF))

    # 5) Annotate sample-level metadata (role, family ID)
    info = hl.import_table(sample_info_path, delimiter='\t', key='s')
    mt_rare = mt_rare.annotate_cols(role = info[mt_rare.s].role,
                                    fam_id = info[mt_rare.s].fam_id)
    
    # 6) Annotate internal allele count (AC) from external cohort resource - For only EAS samples
    internal_ac_ht = hl.read_table(internal_ac_ht_path)
    mt_rare = mt_rare.annotate_rows(
        internal_AC = hl.or_else(internal_ac_ht[mt_rare.locus, mt_rare.alleles].internal_AC, 0))
    mt_ultra_rare = mt_ultra_rare.annotate_rows(
        internal_AC = hl.or_else(internal_ac_ht[mt_ultra_rare.locus, mt_ultra_rare.alleles].internal_AC, 0))

    # 7) Apply cohort-specific frequency filters using internal AC from external resource - For only EAS samples
    mt_rare       = mt_rare.filter_rows(mt_rare.internal_AC < 18)
    mt_ultra_rare = mt_ultra_rare.filter_rows(mt_ultra_rare.internal_AC < 2)
    
    # 8) Checkpoint rare / ultra-rare sets BEFORE Mendelian error filtering
    mt_rare.write(
        os.path.join(output_dir, f'WGS_chr{chr_name}_rare_pre_mendelian.mt'),
        overwrite=True
    )
    mt_ultra_rare.write(
        os.path.join(output_dir, f'WGS_chr{chr_name}_ultra_rare_pre_mendelian.mt'),
        overwrite=True
    )

    # 9) Restrict to trio samples for inheritance pattern analysis - For only trio-based dataset (excluding non-trio cohorts)
    pedigree = hl.Pedigree.read(ped_path)
    trios = hl.array(sum([[x.s, x.mat_id, x.pat_id] for x in pedigree.complete_trios()], []))

    # --------------------------
    # 9-1) Mendelian filtering for rare variants
    # --------------------------
    mt_trio_rare = mt_rare.filter_cols(trios.contains(mt_rare.s))

    tm_rare  = hl.trio_matrix(mt_trio_rare, pedigree, complete_trios=True)
    p_tm_rare = hl.experimental.phase_trio_matrix_by_transmission(tm_rare)

    p_tm_rare = p_tm_rare.annotate_entries(
        transmission_pattern = (
            hl.case()
              # transmission-based patterns from phased proband genotype
              .when(hl.is_defined(p_tm_rare.proband_entry.PBT_GT) & (p_tm_rare.proband_entry.PBT_GT == hl.call(1, 0, phased=True)), 'father')
              .when(hl.is_defined(p_tm_rare.proband_entry.PBT_GT) & (p_tm_rare.proband_entry.PBT_GT == hl.call(0, 1, phased=True)), 'mother')
              .when(hl.is_defined(p_tm_rare.proband_entry.PBT_GT) & (p_tm_rare.proband_entry.PBT_GT == hl.call(1, 1, phased=True)), 'both')
              # parental-only carrier pattern (no phased transmission available for proband)
              .when(
                  hl.is_missing(p_tm_rare.proband_entry.PBT_GT) &
                  (
                      (hl.is_defined(p_tm_rare.father_entry.GT) & p_tm_rare.father_entry.GT.is_non_ref()) |
                      (hl.is_defined(p_tm_rare.mother_entry.GT) & p_tm_rare.mother_entry.GT.is_non_ref())
                  ),
                  'Parental'
              )
              .or_missing()
        )
    )

    valid = ['father', 'mother', 'both', 'Parental']
    mt_rare_filtered = p_tm_rare.filter_entries(hl.array(valid).contains(p_tm_rare.transmission_pattern))

    # 9-2) Mendelian filtering for ultra-rare variants (same logic)
    mt_trio_ultra = mt_ultra_rare.filter_cols(trios.contains(mt_ultra_rare.s))

    tm_ultra  = hl.trio_matrix(mt_trio_ultra, pedigree, complete_trios=True)
    p_tm_ultra = hl.experimental.phase_trio_matrix_by_transmission(tm_ultra)

    p_tm_ultra = p_tm_ultra.annotate_entries(
        transmission_pattern = (
            hl.case()
              .when(hl.is_defined(p_tm_ultra.proband_entry.PBT_GT) & (p_tm_ultra.proband_entry.PBT_GT == hl.call(1, 0, phased=True)), 'father')
              .when(hl.is_defined(p_tm_ultra.proband_entry.PBT_GT) & (p_tm_ultra.proband_entry.PBT_GT == hl.call(0, 1, phased=True)), 'mother')
              .when(hl.is_defined(p_tm_ultra.proband_entry.PBT_GT) & (p_tm_ultra.proband_entry.PBT_GT == hl.call(1, 1, phased=True)), 'both')
              .when(
                  hl.is_missing(p_tm_ultra.proband_entry.PBT_GT) &
                  (
                      (hl.is_defined(p_tm_ultra.father_entry.GT) & p_tm_ultra.father_entry.GT.is_non_ref()) |
                      (hl.is_defined(p_tm_ultra.mother_entry.GT) & p_tm_ultra.mother_entry.GT.is_non_ref())
                  ),
                  'Parental'
              )
              .or_missing()
        )
    )

    mt_ultra_rare_filtered = p_tm_ultra.filter_entries(hl.array(valid).contains(p_tm_ultra.transmission_pattern))

    # 10) Write rare / ultra-rare inherited variant MatrixTables to disk
    mt_rare_filtered.write(
        os.path.join(output_dir, f'WGS_chr{chr_name}_rare_inherited.mt'),
        overwrite=True
    )
    mt_ultra_rare_filtered.write(
        os.path.join(output_dir, f'WGS_chr{chr_name}_ultra_rare_inherited.mt'),
        overwrite=True
    )

    print(f'Done rare & ultra-rare inherited processing for chr{chr_name} : Execution time = {time.time() - start} sec')
    return mt_rare_filtered




def prepare_for_prs_wgs(mt, gwas_sum_path, liftover_chain_path, output_dir, chr_name):
    start = time.time()

    # 1) Variant-level QC for common variants
    mt = mt.filter_rows(
        (mt.variant_qc.call_rate >= 0.95) &
        (mt.variant_qc.p_value_hwe >= 1e-6) &
        (mt.variant_qc.AF[1] > 0.05))

    # 2) Register liftover chain and perform GRCh38 → GRCh37 coordinate liftOver
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg38.add_liftover(liftover_chain_path, rg37)

    mt = mt.annotate_rows(new_locus = hl.liftover(mt.locus, 'GRCh37', include_strand=True))
    mt = mt.filter_rows(hl.is_defined(mt.new_locus) & ~mt.new_locus.is_negative_strand)
    mt = mt.key_rows_by(locus = mt.new_locus.result, alleles = mt.alleles)

    # 3) Reformat long alleles (>13bp) for compatible representation
    mt = mt.annotate_rows(
        new_alleles = hl.array([
            (mt.alleles[0][:13] + '+' + hl.str(mt.alleles[0].length() - 13)) if (mt.alleles[0].length() > 13) else mt.alleles[0],
            (mt.alleles[1][:13] + '+' + hl.str(mt.alleles[1].length() - 13)) if (mt.alleles[1].length() > 13) else mt.alleles[1],
        ]))

    # 4) Remove strand-ambiguous variants (A/T, C/G) to avoid strand issues
    mt = mt.filter_rows(~hl.is_strand_ambiguous(mt.new_alleles[0], mt.new_alleles[1]))

    # 5) Read GWAS summary statistics (on GRCh37) with locus + allele keys
    gwas = hl.import_table(
        gwas_sum_path,
        types={'locus': hl.tlocus(reference_genome='GRCh37')},
        impute=True,
        delimiter='\t',
        force_bgz=True,
        key=['locus', 'A1', 'A2'])

    # 6) Compute cohort AC/AN and align with GWAS AF for frequency harmonization
    rows = mt.rows().key_by('locus', 'alleles')
    rows = rows.annotate(
        AC1     = rows.variant_qc.AC[0],
        AC2     = rows.variant_qc.AC[1],
        N       = rows.variant_qc.AN,
        gwas_AF = gwas[rows.locus, rows.alleles[0], rows.alleles[1]].F)

    # 7) Calculate chi-square statistic to flag AF-discordant variants
    rows = rows.annotate(
        chi_stat = (
            (rows.AC1 - rows.N * rows.gwas_AF) ** 2 / (rows.N * rows.gwas_AF)
        ) + (
            (rows.AC2 - rows.N * (1 - rows.gwas_AF)) ** 2 / (rows.N * (1 - rows.gwas_AF))
        ))

    # 8) Remove variants with outlier AF discrepancy (chi-square > mean + 1 SD)
    stats = rows.aggregate(hl.agg.stats(rows.chi_stat))
    threshold = stats['mean'] + stats['stdev']
    mt = mt.filter_rows(rows[mt.row_key].chi_stat < threshold)

    # 9) Export PRS-ready genotype data in PLINK format (for PRScs / PRS pipeline)
    hl.export_plink(
        mt,
        os.path.join(output_dir, f"WGS_chr{chr_name}_PRS_ready"),
        varid = mt.SNP,
        ind_id = mt.s)

    print(f'Done PRS preparation for chr{chr_name} : Execution time = {time.time() - start} sec')
    return mt




def annotate_and_prioritize_variants(variant_table, vep_config_path, mpc_ht_path, constraint_ht_path, output_dir):
    # 1) Run VEP to obtain functional annotations for all variants
    vt = hl.vep(variant_table, vep_config_path)

    # 2) Extract key transcript-level metadata (gene symbol, transcript ID, gene ID)
    vt = vt.annotate(
        gene_symbol   = vt.vep.transcript_consequences.gene_symbol[0],
        transcript_id = vt.vep.transcript_consequences.transcript_id[0],
        gene_id       = vt.vep.transcript_consequences.gene_id[0])

    # 3) Annotate variant-level MPC scores by locus / alleles / transcript
    mpc = hl.read_table(mpc_ht_path).key_by('locus', 'alleles', 'ENST')
    vt = vt.key_by(vt.locus, vt.alleles, vt.transcript_id)
    vt = vt.annotate(MPC = mpc[vt.locus, vt.alleles, vt.transcript_id].MPC)

    # 4) Annotate gene-level constraint metrics (e.g. pLI, LOEUF) from external constraint resource
    cons = hl.read_table(constraint_ht_path).key_by('gene', 'transcript')
    vt = vt.key_by(vt.gene_symbol, vt.transcript_id)
    vt = vt.annotate(constraint = cons[vt.gene_symbol, vt.transcript_id])
    vt = vt.key_by(vt.locus, vt.alleles)

    # 5) Define protein-truncating variants (PTV) using consequence + LOF annotation from VEP
    ptv = (
        hl.array(['frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained'])
          .contains(vt.vep.most_severe_consequence) &
        (vt.vep.transcript_consequences.lof[0] == 'HC') &
        (
            (vt.vep.transcript_consequences.lof_flags[0] == 'SINGLE_EXON') |
            hl.is_missing(vt.vep.transcript_consequences.lof_flags[0])
        ))

    # 6) Group remaining functional classes: missense, inframe indel, synonymous
    missense      = hl.array(['missense_variant', 'stop_lost', 'start_lost', 'protein_altering_variant']).contains(vt.vep.most_severe_consequence)
    inframe_indel = hl.array(['inframe_deletion', 'inframe_insertion']).contains(vt.vep.most_severe_consequence)
    synonymous    = hl.array(['synonymous_variant', 'stop_retained_variant', 'incomplete_terminal_codon_variant']).contains(vt.vep.most_severe_consequence)

    # 7) Assign a coarse variant category label (VariantType) for downstream prioritization
    vt = vt.annotate(
        VariantType = hl.case()
            .when(ptv,           'protein_truncating')
            .when(missense,      'missense')
            .when(inframe_indel, 'inframe_indel')
            .when(synonymous,    'synonymous')
            .or_missing())

    # 8) Restrict to protein-coding transcripts with a defined variant category
    vt = vt.filter(
        (vt.vep.transcript_consequences.biotype[0] == 'protein_coding') &
        hl.is_defined(vt.VariantType))

    # 9) Write annotated and categorized variant table to disk
    vt.write(os.path.join(output_dir, 'WGS_annotated.ht'), overwrite=True)
    return vt


#########################################################
# Driver script
#########################################################

if __name__ == "__main__":
    hl.init(log="/tmp/hail_pipeline.log")

    vcf_path            = "/path/to/joint_called.vcf.bgz"
    gnomad_ht_path      = "/path/to/gnomad.af.ht"
    lcr_bed_path        = "/path/to/LCR-hs38.bed"
    ped_path            = "/path/to/family.ped"
    sample_info_path    = "/path/to/sample_info.txt"
    internal_ac_ht_path = "/path/to/EAS_internal_AC.ht"
    liftover_chain_path = "/path/to/grch38_to_grch37.over.chain.gz"
    gwas_sum_path       = "/path/to/gwas_sumstats.txt.bgz"
    vep_config_path     = "/path/to/vep_config.json"
    mpc_ht_path         = "/path/to/MPC38.ht"
    constraint_ht_path  = "/path/to/constraint_metrics.ht"

    output_dir = "/path/to/output_dir/"
    chr_name   = "1"

    # QC
    mt_qc = quality_control_wgs(
        vcf_path     = vcf_path,
        output_dir   = output_dir,
        chr_name     = chr_name,
        lcr_bed_path = lcr_bed_path
    )

    # De novo variants
    dnv_ht = call_denovo_wgs(
        mt            = mt_qc,
        ped_path      = ped_path,
        gnomad_ht_path= gnomad_ht_path,
        output_dir    = output_dir,
        chr_name      = chr_name
    )

    # Rare / ultra-rare inherited variants
    mt_inherited_rare = process_rare_inherited_wgs(
        mt                 = mt_qc,
        gnomad_ht_path     = gnomad_ht_path,
        ped_path           = ped_path,
        sample_info_path   = sample_info_path,
        internal_ac_ht_path= internal_ac_ht_path,
        output_dir         = output_dir,
        chr_name           = chr_name
    )

    # PRS preparation
    mt_prs = prepare_for_prs_wgs(
        mt                 = mt_qc,
        gwas_sum_path      = gwas_sum_path,
        liftover_chain_path= liftover_chain_path,
        output_dir         = output_dir,
        chr_name           = chr_name
    )

    # Functional annotation and prioritization (example: de novo variants)
    annotated_ht = annotate_and_prioritize_variants(
        variant_table    = dnv_ht,
        vep_config_path  = vep_config_path,
        mpc_ht_path      = mpc_ht_path,
        constraint_ht_path=constraint_ht_path,
        output_dir       = output_dir
    )

    print("Pipeline completed.")