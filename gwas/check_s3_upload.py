import pandas as pd
import os
import json
import pathlib
from pathlib import Path
import re
import pprint as pp
# open path
# s3cmd ls s3://interval-15x-matrixtables/WGS_final_march_2020_dbsnp_v53.mt/entries/rows/parts/ > entries_rows.txt
# /opt/sanger.ac.uk/hgi/hail/tmp/entries_rows.txt


def find_missing(lst):
    return [x for x in range(lst[0], lst[-1]+1)
            if x not in lst]


#filepaths = "/Users/pa10/Programming/hail_pipeline_github/entries_rows.txt"
filepaths = "/opt/sanger.ac.uk/hgi/hail/tmp/scripts/hail-pipelines-internal/entries_rows.txt"
#allpaths = [p for p in pathlib.Path(plots_path).iterdir() if p.is_file()]
f = open(filepaths, "r")
lines = f.readlines()
s3indices = []
for x in lines:
    m = re.search(
        r'(.*)s3\:\/\/interval-15x-matrixtables\/WGS\_final\_march\_2020\_dbsnp\_v53\.mt\/entries\/rows\/parts/part-(\d+)-(\d+)-(\d+)-(.*)', x)
    if m:
        # print(m.group(2))
        s3indices.append(int(m.group(2)))

f.close()

missing_files = find_missing(s3indices)
pp.pprint(missing_files)
print(len(missing_files))
# print(allpaths)
# print(filename)


#data = list(filter(None, data))
# print(data)
# with open('/Users/pa10/Programming/gwas_display_app/public/data/data.json', 'w') as fout:
#    json.dump(data, fout)
