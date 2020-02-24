import os
from datetime import datetime
import yaml
import hail as hl
import re
from bokeh.plotting import output_file, save
import json

project_root = os.path.dirname(os.path.dirname(__file__))
print(project_root)


snp_vqsr_threshold= -0.6647
indel_vqsr_threshold=-1.2537



BUCKET = "gs://interval-wes"
#Define chromosome here
tmp_dir="/Users/pa10/Programming/google-code/google/tmp"


if __name__ == "__main__":
    #need to create spark cluster first before intiialising hail
    #Define the hail persistent storage directory
    hl.init(default_reference="GRCh38", tmp_dir=tmp_dir)
    mt=hl.read_matrix_table('gs://interval-wes/WES_full-final-QCed.mt')
    db = hl.experimental.DB() 
    mt = db.annotate_rows_db(mt, "DANN", "CADD", "Ensembl_homo_sapiens_reference_genome", "Ensembl_homo_sapiens_low_complexity_regions", "gencode", "gnomad_lof_metrics", "gnomad_exome_sites", "gnomad_genome_sites", "clinvar_gene_summary", "clinvar_variant_summary", "dbNSFP_variants", "dbNSFP_genes", "gerp_scores", "gerp_elements")
    mt=mt.checkpoint(f"{BUCKET}/WES_full-final-QCed-annotated.mt", overwrite=True)
