import os
import hail as hl
import sys
import json



project_root = os.path.dirname(os.path.dirname(__file__))
print(project_root)


BUCKET = "gs://elgh-ddd"
#Define chromosome here
tmp_dir="/Users/pa10/Programming/google-code/google/tmp"

if __name__ == "__main__":
    #need to create spark cluster first before intiialising hail
    #Define the hail persistent storage directory
    hl.init(default_reference="GRCh38", tmp_dir=tmp_dir)


    chrX_haploid_vcf_path= "f{BUCKET}/chrX-haploid.vcf.bgz"
    chrX_diploid_vcf_path= "f{BUCKET}/chrX-diploid.vcf.bgz"
    chrY_haploid_vcf_path= "f{BUCKET}/chrY-haploid.vcf.bgz"
    chrY_diploid_vcf_path= "f{BUCKET}/chrY-diploid.vcf.bgz"

    print("Read in all VCFs:")
    mtX_hap = hl.import_vcf(chrX_haploid_vcf_path, force_bgz=True, reference_genome='GRCh38')
    mtX_dip = hl.import_vcf(chrX_diploid_vcf_path, force_bgz=True, reference_genome='GRCh38')
    mtY_hap = hl.import_vcf(chrY_haploid_vcf_path, force_bgz=True, reference_genome='GRCh38')
    mtY_dip = hl.import_vcf(chrY_diploid_vcf_path, force_bgz=True, reference_genome='GRCh38')
    ############################
    CHROMOSOME = "chrX"
    print(f"Reading {CHROMOSOME} diploid VCF calls")
    mtX_dip_split = hl.split_multi_hts(mtX_dip, keep_star = False)
    mtX_hap_split = hl.split_multi_hts(mtX_hap, keep_star = False)


    ####### sex imputation
    print("Sex imputation:")
    mtx_unphased = mtX_dip_split.select_entries(GT=hl.unphased_diploid_gt_index_call(mtX_dip_split.GT.n_alt_alleles()))
    imputed_sex = hl.impute_sex(mtx_unphased.GT)


    #Annotate samples male or female:
    #mt = mt_split.annotate_cols(
    #    sex=hl.cond(imputed_sex[mtX_hap_split.s].is_female, "female","male")
    #
    #

    #Haploid VCF: filter to male samples only
    print("Haploid chrX VCF: filter to male samples only")
    mtX_hap_males = mtX_hap_split.filter_cols(imputed_sex[mtX_hap_split.s].is_female, keep=False)
    #Diploid vcf: one matrixtable only Female, one matrixtable only male
    print("Diploid chrX vcf: create new matrixtable one for males one for female samples")
    mtX_dip_males = mtX_dip_split.filter_cols(imputed_sex[mtX_dip_split.s].is_female, keep = False)
    mtX_dip_females = mtX_dip_split.filter_cols(imputed_sex[mtX_dip_split.s].is_female, keep = True)

    #Get par region
    rg = hl.get_reference('GRCh38')
    par = rg.par

    #For haploid male chrX variants remove all falling within par
    mtX_hap_males_NONPAR = hl.filter_intervals(mtX_hap_males, par, keep=False)

    #For diploid male samples chrX variants KEEP PAR regions only
    mtX_dip_males_PAR = hl.filter_intervals(mtX_dip_males, par, keep = True)

    #The hap and dip VCF files differ in entry fielse. Dropping these fields to allow union_rows and union_cols to proceed
    fields_to_drop = ['PGT', 'PID', 'PS']

    mtX_dip_males_PAR_dropf = mtX_dip_males_PAR.drop(*fields_to_drop)
    mtX_dip_females_dropf = mtX_dip_females.drop(*fields_to_drop)
    mtX_union_males = mtX_hap_males_NONPAR.union_rows(mtX_dip_males_PAR_dropf)


    mt_final = mtX_union_males.union_cols(mtX_dip_females_dropf)
    print("Output matrixtable:")
    mt_final = mt_final.checkpoint(
        f"{BUCKET}/elgh-ddd-chrX_surgery.mt", overwrite=True)
    print("Export VCF")
    hl.export_vcf(mt_final, f"{BUCKET}/elgh-ddd-chrX-surgery.vcf.bgz")
