
'''
author: Pavlos Antoniou
date: 22/07/19
'''

import os
import hail as hl
import pyspark
import json
import sys




project_root = os.path.dirname(os.path.dirname(__file__))
print(project_root)

s3credentials = os.path.join(project_root, "config_files/s3_credentials.json")
print(s3credentials)

storage = os.path.join(project_root , "config_files/storage.json")

thresholds = os.path.join(project_root, "config_files/thresholds.json")

with open(f"{s3credentials}", 'r') as f:
    credentials = json.load(f)

with open(f"{storage}", 'r') as f:
    storage = json.load(f)

with open(f"{thresholds}", 'r') as f:
    thresholds = json.load(f)

def import_vcf(vcf_path,partitions):
    mt = hl.import_vcf(vcf_path,
                            force_bgz=True,
                            reference_genome='GRCh38',
                            array_elements_required=False,
                            min_partitions=partitions)

    if mt.n_partitions() > partitions:
        mt = mt.naive_coalesce(partitions)
    return mt


def fix_genotpe(mt):

    pl_expr = (
        hl.zip_with_index(mt.PL).flatmap(lambda x: hl.range(x[0] + 1).map(lambda i: hl.cond(i == x[0], x[1], 999))))

    gt_allele = mt.GT.n_alt_alleles()
    gt_and_pl = (hl.switch(mt.GT.ploidy).when(1, (hl.call(gt_allele, gt_allele), pl_expr)).default((mt.GT, mt.PL)))

    mt = mt.annotate_entries(GT=gt_and_pl[0], PL=gt_and_pl[1])

    return mt

if __name__ == "__main__":
    #need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    #Define the hail persistent storage directory
    tmp_dir = os.path.join(os.environ["HAIL_HOME"], "tmp")
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    #s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])
    
    
    sample_QC_nonHail = hl.import_table(storage["intervalwes"]["s3"]["sex_samples"], impute=True,delimiter=",")


    partitions=450
    haploid_vcf_path= storage["intervalwes"]["s3"]["haploid_vcf"]
    diploid_vcf_path= storage["intervalwes"]["s3"]["diploid_vcf"]

    print("Read in all VCFs:")
    mt_hap=import_vcf(haploid_vcf_path,partitions)
    mt_dip=import_vcf(diploid_vcf_path,partitions)

    mt_hap= mt_hap.checkpoint(f"{tmp_dir}/intervalwes/haploidXY.mt", overwrite=True)
    mt_dip= mt_dip.checkpoint(f"{tmp_dir}/intervalwes/diploidXY.mt", overwrite=True)

    mtX_hap=mt_hap.filter_rows(mt_hap.locus.contig == 'chrX')
    mtY_hap=mt_hap.filter_rows(mt_hap.locus.contig == 'chrY')
    mtX_dip=mt_dip.filter_rows(mt_dip.locus.contig == 'chrX')
    mtY_dip=mt_dip.filter_rows(mt_dip.locus.contig == 'chrY')
    mtX_hap= mtX_hap.checkpoint(f"{tmp_dir}/intervalwes/haploidX.mt", overwrite=True)
    mtX_dip= mtX_dip.checkpoint(f"{tmp_dir}/intervalwes/diploidX.mt", overwrite=True)
    mtY_hap= mtY_hap.checkpoint(f"{tmp_dir}/intervalwes/haploidY.mt", overwrite=True)
    mtY_dip= mtY_dip.checkpoint(f"{tmp_dir}/intervalwes/diploidY.mt", overwrite=True)
    print("Fixing chrX and chrY GT haploid")
    mtX_hap=fix_genotpe(mtX_hap)
    mtY_hap=fix_genotpe(mtY_hap)
    mtX_hap = mtX_hap.checkpoint(f"{tmp_dir}/intervalwes/chrX_haploid_GT_fixed.mt", overwrite=True)
    mtY_hap = mtY_hap.checkpoint(f"{tmp_dir}/intervalwes/chrY_haploid_GT_fixed.mt", overwrite=True)


    #####################################################################
    ###################### chromosome Y  ##############################
    #####################################################################    
    print("chrY  surgery")
    mtY_hap_split = hl.split_multi_hts(mtY_hap, keep_star=False)
    call1 = hl.eval(hl.call(0, 0))
    mtY_hap_sex_annot = mtY_hap_split.annotate_cols(sex=sample_QC_nonHail.key_by("ega_sample_name")[mtY_hap_split.s].gender)
    mtY_hap_sex_annot = mtY_hap_sex_annot.annotate_entries(
        GT=hl.cond(mtY_hap_sex_annot.sex == 2, call1, mtY_hap_sex_annot.GT))
    


    info_fields_to_drop = ['AS_BaseQRankSum',
                           'AS_FS',
                           'AS_InbreedingCoeff',
                           'AS_MQ',
                           'AS_MQRankSum',
                           'AS_QD',
                           'AS_ReadPosRankSum',
                           'AS_SOR',
                           'RAW_MQandDP',
                           'DB']
    #info1 = mtY_hap_sex_annot.info.drop(*info_fields_to_drop)
    #mtY_hap_sex_annot = mtY_hap_sex_annot.annotate_rows(info=info1)
    print("Checkpoint chrY")
    mtY_hap_sex_annot= mtY_hap_sex_annot.checkpoint(f"{tmp_dir}/intervalwes/chrY_haploid_final.mt", overwrite=True)
    print("Export chrY VCF")
    hl.export_vcf(mtY_hap_sex_annot, f"{tmp_dir}/intervalwes/chrY-surgery_final.vcf.bgz")

    #####################################################################
    ###################### chromosome X  ##############################
    #####################################################################
    CHROMOSOME = "chrX"
    print(f"Split multi_hts {CHROMOSOME} diploid VCF calls")
    mtX_dip_split = hl.split_multi_hts(mtX_dip, keep_star = False)
    mtX_hap_split = hl.split_multi_hts(mtX_hap, keep_star = False)
    mtX_dip_split = mtX_dip_split.checkpoint(f"{tmp_dir}/intervalwes/chrX_diploid_split.mt", overwrite=True)
    mtX_hap_split = mtX_hap_split.checkpoint(f"{tmp_dir}/intervalwes/chrX_haploid_split.mt", overwrite=True)

    #New sex imputation code from here"

    ####### sex imputation
    print("Sex from Sample QC:")
    mtX_dip_split_annotated = mtX_dip_split.annotate_cols(sex=sample_QC_nonHail.key_by("ID")[mtX_dip_split.s].SEX)
    mtX_hap_split_annotated= mtX_hap_split.annotate_cols(sex=sample_QC_nonHail.key_by("ID")[mtX_hap_split.s].SEX)
    #mtx_unphased = mtX_dip_split.select_entries(GT=hl.unphased_diploid_gt_index_call(mtX_dip_split.GT.n_alt_alleles()))
    #imputed_sex = hl.impute_sex(mtx_unphased.GT)



    #Haploid VCF: filter to male samples only
    print("Haploid chrX VCF: filter to male samples only")
    #mtX_hap_males = mtX_hap_split.filter_cols(imputed_sex[mtX_hap_split.s].is_female, keep=False)
    mtX_hap_males = mtX_hap_split_annotated.filter_cols(mtX_hap_split_annotated.sex==1, keep=True)
    info1 = mtX_hap_males.info.drop(*info_fields_to_drop)
    mtX_hap_males = mtX_hap_males.annotate_rows(info=info1)
    mtX_hap_males = mtX_hap_males.checkpoint(f"{tmp_dir}/intervalwes/chrX_haploid_males.mt", overwrite=True)

    #Diploid vcf: one matrixtable only Female, one matrixtable only male
    print("Diploid chrX vcf: create new matrixtable one for males one for female samples")
    mtX_dip_males = mtX_dip_split_annotated.filter_cols(mtX_dip_split_annotated.sex==1, keep = True)
    mtX_dip_females = mtX_dip_split_annotated.filter_cols(mtX_dip_split_annotated.sex==1, keep=False)
    mtX_dip_males = mtX_dip_males.checkpoint(f"{tmp_dir}/intervalwes/chrX_diploid_males.mt", overwrite=True)
    mtX_dip_females = mtX_dip_females.checkpoint(f"{tmp_dir}/intervalwes/chrX_diploid_females.mt", overwrite=True)


    #Get par region
    rg = hl.get_reference('GRCh38')
    par = rg.par

    #For haploid male chrX variants remove all falling within par
    mtX_hap_males_NONPAR = hl.filter_intervals(mtX_hap_males, par, keep=False)

    #For diploid male samples chrX variants KEEP PAR regions only
    mtX_dip_males_PAR = hl.filter_intervals(mtX_dip_males, par, keep = True)

    #The hap and dip VCF files differ in entry fielse. Dropping these fields to allow union_rows and union_cols to proceed

    fields_to_drop = ['PGT', 'PID','PS']
    mtX_dip_males_PAR_dropf = mtX_dip_males_PAR.drop(*fields_to_drop)
    mtX_dip_females_dropf = mtX_dip_females.drop(*fields_to_drop)

    mtX_union_males = mtX_hap_males_NONPAR.union_rows(mtX_dip_males_PAR_dropf)

    # fields_to_drop_secondmt = ['ClippingRankSum', 'RAW_MQ']
    # info1 = mtX_hap_males_NONPAR.info.drop(*info_fields_to_drop)
    # info2 = mtX_dip_males_PAR_dropf.info.drop(*fields_to_drop_secondmt)

    # mtX_hap_males_NONPAR = mtX_hap_males_NONPAR.annotate_rows(info=info1)
    # mtX_dip_males_PAR_dropf = mtX_dip_males_PAR_dropf.annotate_rows(info=info2)
    # mtX_dip_males_PAR_dropf = mtX_dip_males_PAR_dropf.annotate_entries(SB=mtX_dip_males_PAR_dropf.SB.map(lambda x: x))
    # mtX_union_males = mtX_hap_males_NONPAR.union_rows(mtX_dip_males_PAR_dropf)

    mtX_dip_females_dropf = mtX_dip_females_dropf.annotate_entries(SB=mtX_dip_females_dropf.SB.map(lambda x: x))
    mt_final = mtX_union_males.union_cols(mtX_dip_females_dropf)

    #Sort new merged samples chrX sample order like it was before (mtX_hap was one of the original mts)
    samples=mtX_hap.s.collect()
    name_to_index = {col['s']: col['col_idx'] for col in 
                 mt_final.add_col_index().cols().key_by().select('col_idx', 's').collect()}
    mt_final = mt_final.choose_cols([name_to_index[s] for s in samples])

    print("Output matrixtable:")
    mt_final = mt_final.checkpoint(
        f"{tmp_dir}/intervalwes/chrX_surgery_FINAL.mt", overwrite=True)
    #print("Export VCF")
    #hl.export_vcf(mt_final, f"{tmp_dir}/elgh-ddd/elgh-ddd-chrX-surgery_final.vcf.bgz")
   
    #join chrX and chrY together
    XYfinal=mt_final.union_rows(mtY_hap_sex_annot)

    #Join with WES autosomes
    wes=hl.read_matrix_table(f"{tmp_dir}/intervalwes/WES_vqslod_split-multi_checkpoint.mt")
    wes_autosomes=wes.filter_rows((wes.locus.contig == 'chrX') & (wes.locus.contig == 'chrY'), keep=False )
    fields_to_drop = ['PGT', 'PID','PS']
    wes_autosomes=wes_autosomes.drop(*fields_to_drop)
    #multi split wes:

    #join with wes
    wes_autosomes=wes.filter_rows((wes.locus.contig == 'chrX') & (wes.locus.contig == 'chrY'), keep=False )
    wes_new=wes_autosomes.union_rows(XYfinal)

    wes_final=wes_new.checkpoint(f"{tmp_dir}/intervalwes/WES_AFTER_surgery.mt", overwrite=True)