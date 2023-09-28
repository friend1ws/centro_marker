#! /bin/bash

<<_
for pat in `aws s3 ls --recursive s3://human-pangenomics/submissions | grep "_trio.paternal.fa.gz" | grep verkko_1.3.1 | cut -f 5 -d ' '`; 
do 
    echo "s3://human-pangenomics/"${pat}
done > verkko_1.3.1_pat_list.txt

echo -n > verkko_1.3.1_mat_list.txt 
while read line; 
do 
    echo ${line%%.paternal.fa.gz}.maternal.fa.gz >> verkko_1.3.1_mat_list.txt
done < verkko_1.3.1_pat_list.txt 
_

# s3://human-pangenomics/submissions/53FEE631-4264-4627-8FB6-09D7364F4D3B--ASM-COMP/HG002/assemblies/verkko_1.3.1/trio/assembly_verkko_v1.3_trio.paternal.fa.gz

echo -n > verkko_1.3.1_sample_list.txt
while read line; 
do 
    line2=${line%%/assemblies*} 
    line3=${line2##s3://human-pangenomics/submissions/53FEE631-4264-4627-8FB6-09D7364F4D3B--ASM-COMP/}
    echo $line3 >> verkko_1.3.1_sample_list.txt
done < verkko_1.3.1_pat_list.txt 

paste verkko_1.3.1_sample_list.txt verkko_1.3.1_pat_list.txt verkko_1.3.1_mat_list.txt > verkko_1.3.1_list.txt

# rm -rf hifiasm_v0.14_raw_pat_list.txt
# rm -rf hifiasm_v0.14_raw_mat_list.txt

