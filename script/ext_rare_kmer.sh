#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log -o log
#$ -l s_vmem=32G 

module load /usr/local/package/modulefiles/singularity/3.7.0

CEN_ID=$1

if [ "$CEN_ID" -eq 23 ]; then
    CEN_ID="X"
fi

if [ ! -d ../output/rare_kmer/chr${CEN_ID}/hifiasm_v0.14_raw ]
then
    mkdir -p ../output/rare_kmer/chr${CEN_ID}/hifiasm_v0.14_raw
fi

if [ ! -d ../output/rare_kmer/chr${CEN_ID}/hifiasm_v0.19.5 ]
then
    mkdir -p ../output/rare_kmer/chr${CEN_ID}/hifiasm_v0.19.5
fi


<<_

for cfile in `ls ../output/cen_seq/hifiasm_v0.14_raw/*/chr${CEN_ID}/*.chm13.cen.filt.fa`; 
do 
    if [[ $cfile == *HG002-full-0.14* ]]
    then
        continue
    fi

    if [ -s $cfile ]; 
    then 
        cp $cfile ../output/rare_kmer/chr${CEN_ID}/hifiasm_v0.14_raw
        cp ${cfile%chm13.cen.filt.fa}sd_info.txt ../output/rare_kmer/chr${CEN_ID}/hifiasm_v0.14_raw
    fi 
done

for cfile in `ls ../output/cen_seq/hifiasm_v0.19.5/*/chr${CEN_ID}/*.chm13.cen.filt.fa`;
do
    if [ -s $cfile ];
    then
        cp $cfile ../output/rare_kmer/chr${CEN_ID}/hifiasm_v0.19.5
        cp ${cfile%chm13.cen.filt.fa}sd_info.txt ../output/rare_kmer/chr${CEN_ID}/hifiasm_v0.19.5
    fi 

done

python3 rare_kmer_parse.py ../output/rare_kmer/chr${CEN_ID} ../output/rare_kmer/chr${CEN_ID}/merged.rare_kmer.txt


singularity exec ~/image/kmer_utils_0.1.2.sif python3 filt_kmer.py ../output/rare_kmer/chr${CEN_ID}/merged.rare_kmer.txt ../output/rare_kmer/chr${CEN_ID} ../output/rare_kmer/chr${CEN_ID}/merged.rare_kmer.pruned.txt
_

python3 add_rare_kmer_info.py ../output/rare_kmer/chr${CEN_ID}/merged.rare_kmer.pruned.txt ../output/rare_kmer/chr${CEN_ID} ../output/rare_kmer/chr${CEN_ID}/merged.rare_kmer.pruned.annot.txt




