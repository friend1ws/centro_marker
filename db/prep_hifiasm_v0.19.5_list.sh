#! /bin/bash

for pat in `aws s3 ls --recursive s3://human-pangenomics/submissions | grep ".pat.fa.gz" | grep hifiasm_v0.19.5 | cut -f 5 -d ' '`; 
do 
    echo "s3://human-pangenomics/"${pat}
done > hifiasm_v0.19.5_pat_list.txt

echo -n > hifiasm_v0.19.5_mat_list.txt 
while read line; 
do 
    echo ${line%%.pat.fa.gz}.mat.fa.gz >> hifiasm_v0.19.5_mat_list.txt
done < hifiasm_v0.19.5_pat_list.txt 


echo -n > hifiasm_v0.19.5_sample_list.txt
while read line; 
do 
    line2=${line%%/assemblies*} 
    line3=${line2##s3://human-pangenomics/submissions/53FEE631-4264-4627-8FB6-09D7364F4D3B--ASM-COMP/}
    echo $line3 >> hifiasm_v0.19.5_sample_list.txt
done < hifiasm_v0.19.5_pat_list.txt 

paste hifiasm_v0.19.5_sample_list.txt hifiasm_v0.19.5_pat_list.txt hifiasm_v0.19.5_mat_list.txt > hifiasm_v0.19.5_list.txt

# rm -rf hifiasm_v0.14_raw_pat_list.txt
# rm -rf hifiasm_v0.14_raw_mat_list.txt

