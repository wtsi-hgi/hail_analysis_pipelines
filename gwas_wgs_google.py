import os
from datetime import datetime
import yaml
import hail as hl
import re
from bokeh.plotting import output_file, save
import json

project_root = os.path.dirname(os.path.dirname(__file__))
print(project_root)


storage = os.path.join(project_root , "config-files/storage.json")

thresholds = os.path.join(project_root, "config-files/thresholds.json")

with open(f"{storage}", 'r') as f:
    storage = json.load(f)

with open(f"{thresholds}", 'r') as f:
    thresholds = json.load(f)

BUCKET = storage['intervalwgs']['google']['bucket']

tmp_dir=storage['intervalwgs']['google']['tmp_dir']


CHROMOSOMES = ["chr2",
               "chr3",
               "chr4",
               "chr5",
               "chr6",
               "chr7",
               "chr8",
               "chr9",
               "chr10",
               "chr11",
               "chr12",
               "chr13",
               "chr14",
               "chr15",
               "chr16",
               "chr17",
               "chr18",
               "chr19",
               "chr20",
               "chr21",
               "chr22",
               "chrX",
               "chrY"
               ]
tmp_dir=storage['intervalwgs']['google']['tmp_dir']


covariates =[
    "neut_p_gwas_normalised",
    "eo_p_gwas_normalised",
    "mono_p_gwas_normalised",
   "lymph_p_gwas_normalised",
    "baso_p_gwas_normalised",
   "ret_p_gwas_normalised",
   "hlr_p_gwas_normalised",
   "hct_gwas_normalised",
   "pct_gwas_normalised",
   "hgb_gwas_normalised",
   "rbc_gwas_normalised",
   "wbc_gwas_normalised",
   "mpv_gwas_normalised",
   "plt_gwas_normalised",
   "rdw_gwas_normalised",
   "pdw_gwas_normalised",
   "mcv_gwas_normalised",
   "mch_gwas_normalised",
   "mchc_gwas_normalised",
   "ret_gwas_normalised",
   "hlr_gwas_normalised",
   "neut_gwas_normalised",
   "mono_gwas_normalised",
   "baso_gwas_normalised",
   "eo_gwas_normalised",
   "lymph_gwas_normalised",
   "irf_gwas_normalised",
   "myeloid_wbc_gwas_normalised",
   "gran_gwas_normalised",
   "eo_baso_sum_gwas_normalised",
   "neut_eo_sum_gwas_normalised",
   "baso_neut_sum_gwas_normalised",
   "gran_p_myeloid_wbc_gwas_normalised",
   "eo_p_gran_gwas_normalised",
   "neut_p_gran_gwas_normalised",
   "baso_p_gran_gwas_normalised"

]




if __name__ == "__main__":
    #need to create spark cluster first before intiialising hail
    #Define the hail persistent storage directory
    hl.init(default_reference="GRCh38", tmp_dir=tmp_dir)

    VQSLOD_snps = hl.import_table(f"{BUCKET}/qc-files/VQSLOD_snps.bgz",
                                  types={"Locus": "locus<GRCh38>", "VQSLOD": hl.tfloat64})
    VQSLOD_indels = hl.import_table(f"{BUCKET}/qc-files/VQSLOD_indels.bgz",
                                    types={"Locus": "locus<GRCh38>", "VQSLOD": hl.tfloat64})
    sample_QC_nonHail = hl.import_table(f"{BUCKET}/qc-files/INTERVAL_WGS_Sample_QC_28-05-2019_fixed_header.txt",
                                        impute=True)
    gws_gwa_map = hl.import_table(f"{BUCKET}/qc-files/WGS-2-GWA_omicsMap.txt", impute=True)

    full_blood_count = hl.import_table(f"{BUCKET}/qc-files/Int50kSampleId.tsv", impute=True, key='Wgs_RAW_bl')

    print('Joining annotations')
    ja = sample_QC_nonHail.key_by('ID')
    gws_gwa_map = gws_gwa_map.key_by('unique_WGS_ID')

    ja = ja.annotate(gws_gwa_map=gws_gwa_map[ja['ID']])
    fbc = full_blood_count.key_by('Wgs_RAW_bl')
    ja = ja.annotate(fbc=full_blood_count[ja.gws_gwa_map.Wgs_RAW_bl])

    mt_chr1 = hl.read_matrix_table(f"{BUCKET}/matrixtables/chr1/chr1-sampleqc-variantqc-filtered-FINAL.mt")

    for CHROMOSOME in CHROMOSOMES:
        mt = hl.read_matrix_table(f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}-sampleqc-variantqc-filtered-FINAL.mt")
        mt_chr1 = mt_chr1.union_rows(mt)

    CHROMOSOME = "WGS"

    print("Annotating matrixtable with fbc:")
    mt = mt_chr1.annotate_cols(sample_qc_and_phenotype=ja[mt_chr1.s])

    print("Linear regression")
    # TIM NOTE: I changed the below to show how linear_regression_rows can process groups of phenotypes in a vectorized way
    gwas = hl.linear_regression_rows(
        y=[
            mt.sample_qc_and_phenotype.fbc.neut_p_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.eo_p_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.mono_p_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.lymph_p_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.baso_p_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.ret_p_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.hlr_p_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.hct_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.pct_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.hgb_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.rbc_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.wbc_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.mpv_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.plt_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.rdw_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.pdw_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.mcv_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.mch_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.mchc_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.ret_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.hlr_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.neut_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.mono_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.baso_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.eo_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.lymph_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.irf_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.myeloid_wbc_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.gran_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.eo_baso_sum_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.neut_eo_sum_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.baso_neut_sum_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.gran_p_myeloid_wbc_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.eo_p_gran_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.neut_p_gran_gwas_normalised,
            mt.sample_qc_and_phenotype.fbc.baso_p_gran_gwas_normalised
            #                   mt.sample_qc_and_phenotype.fbc.wbc_gwas_normalised,
            #                   mt.sample_qc_and_phenotype.fbc.mpv_gwas_normalised],
        ],
        x=mt.GT.n_alt_alleles(), covariates=[1.0], pass_through=[mt.rsid])

    # gwas = hl.linear_regression_rows(y=[[mt_filtered.sample_qc_and_phenotype.wbc_gwas_normalised, ...], [family2...]],
    print("Linear regression CHECKPOINT")
    # TIM NOTE: checkpoint here to prevent multiple execution (write to a file, read that file)
    gwas = gwas.checkpoint(f"{BUCKET}/gwas/{CHROMOSOME}-gwasfbc-checkpoint", overwrite=True)
    # gwas = gwas.checkpoint(s3"{tmp_dir}/gwas_wbc_chr19_checkpoint.mt")
    print("Linear regression output table")
    gwas.export(f"{BUCKET}/gwas/gwas-{CHROMOSOME}-export.tsv.bgz", header=True)

    print("Plotting")

    for i in range(0, 36):
        print(f"Plotting {i}:{covariates[i]}")
        p = hl.plot.manhattan(gwas.p_value[i], title=f"Interval WGS GWAS Manhattan Plot: {covariates[i]}")
        output_file(f"{i}.WGS-manhattan-{covariates[i]}.html")
        save(p)
        hl.hadoop_copy(f"{i}.WGS-manhattan-{covariates[i]}.html", f"{BUCKET}/gwas/plots/")

       # p = hl.plot.qq(gwas.p_value[i], title=f"Interval WGS GWAS QQ Plot: {covariates[i]}")
        # export_png(p,f"{BUCKET}/output-tables/wgs-qq-{covariates[i]}.html")
       # output_file(f"{i}.WGS-qq-{covariates[i]}.html")
       # save(p)
       # hl.hadoop_copy(f"{i}.WGS-qq-{covariates[i]}.html", f"{BUCKET}/gwas/plots/")

    mt = mt.checkpoint(f"{BUCKET}/matrixtables/{CHROMOSOME}/{CHROMOSOME}-sampleqc-variantqc-filtered-FINAL.mt", overwrite = True)