#! /bin/bash

for pat in `aws s3 ls --recursive s3://human-pangenomics/working/HPRC | grep ".pat.fa.gz" | grep hifiasm_v0.14_raw | cut -f 5 -d ' '`; do echo "s3://human-pangenomics/"${pat}; done > hifiasm_v0.14_raw_pat_list.txt

while read line; do echo ${line%%.pat.fa.gz}.mat.fa.gz; done < hifiasm_v0.14_raw_pat_list.txt > hifiasm_v0.14_raw_mat_list.txt

while read line; do line2=${line%%/assemblies*}; line3=${line2##s3://human-pangenomics/working/HPRC/}; line4=${line3##s3://human-pangenomics/working/HPRC_PLUS/}; echo $line4; done < hifiasm_v0.14_raw_pat_list.txt > hifiasm_v0.14_raw_sample_list.txt

paste hifiasm_v0.14_raw_sample_list.txt hifiasm_v0.14_raw_pat_list.txt hifiasm_v0.14_raw_mat_list.txt |  awk 'BEGIN {OFS = "\t"} { if ($1 == "HG002") {$1 = "HG002-full-0.14"}; {print}}' > hifiasm_v0.14_raw_list.txt

rm -rf hifiasm_v0.14_raw_pat_list.txt
rm -rf hifiasm_v0.14_raw_mat_list.txt

