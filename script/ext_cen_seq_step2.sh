#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log -o log
#$ -l s_vmem=16G -pe def_slot 2

set -xe

module load /usr/local/package/modulefiles/singularity/3.7.0

TASK_ID=$SGE_TASK_ID
# TASK_ID=1
PLATFORM=$1
CEN_ID=$2

OUT_DIR=../output/cen_seq/${PLATFORM}
LIST_FILE=../db/${PLATFORM}_list.txt


SAMPLE=`head -n ${TASK_ID} ${LIST_FILE} | tail -n 1 | cut -f 1` 

if [ ! -d ${OUT_DIR}/${SAMPLE} ]
then
    mkdir -p ${OUT_DIR}/${SAMPLE}
fi




tcen_chr=`head -n ${CEN_ID} ../db/cen_list.proc.bed | tail -n 1 | cut -f 1`
tcen_start=`head -n ${CEN_ID} ../db/cen_list.proc.bed | tail -n 1 | cut -f 2`
tcen_end=`head -n ${CEN_ID} ../db/cen_list.proc.bed | tail -n 1 | cut -f 3`

tcen_id=${tcen_chr#chr}

if [ ! -d ${OUT_DIR}/${SAMPLE}/${tcen_chr} ]
then
    mkdir -p ${OUT_DIR}/${SAMPLE}/${tcen_chr}
fi 


<<_
python3 ext_cen.py ${OUT_DIR}/${SAMPLE}/${SAMPLE}.pat.chm13.paf ${tcen_chr} ${tcen_start} ${tcen_end} > ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.pat.chm13.cen.namelist.txt 
python3 ext_cen.py ${OUT_DIR}/${SAMPLE}/${SAMPLE}.mat.chm13.paf ${tcen_chr} ${tcen_start} ${tcen_end} > ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.mat.chm13.cen.namelist.txt 

~/bin/seqtk/seqtk subseq ${OUT_DIR}/${SAMPLE}/${SAMPLE}.pat.fa.gz ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.pat.chm13.cen.namelist.txt > ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.pat.chm13.cen.fa
~/bin/seqtk/seqtk subseq ${OUT_DIR}/${SAMPLE}/${SAMPLE}.mat.fa.gz ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.mat.chm13.cen.namelist.txt > ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.mat.chm13.cen.fa 


singularity exec ~/image/stringdecomposer_1.1.2--py310h30d9df9_1.sif stringdecomposer -t 2 ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.pat.chm13.cen.fa ../db/HORmonData20211006/MonomersFinal/cen${tcen_id}_monomers.fa -o ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.pat.stringdecomposer
singularity exec ~/image/stringdecomposer_1.1.2--py310h30d9df9_1.sif stringdecomposer -t 2 ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.mat.chm13.cen.fa ../db/HORmonData20211006/MonomersFinal/cen${tcen_id}_monomers.fa -o ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.mat.stringdecomposer

_
python3 get_sd_info.py ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.pat.stringdecomposer/final_decomposition.tsv > ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.pat.sd_info.txt
python3 get_sd_info.py ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.mat.stringdecomposer/final_decomposition.tsv > ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.mat.sd_info.txt


echo -n > ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.pat.chm13.cen.filt.fa
pcount=`awk -F'\t' 'BEGIN {count=0} { if ($6 == "True") count++} END {print count}' ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.pat.sd_info.txt`
if [ $pcount = 1 ]
then

    while read line
    do
        tcontig=`echo $line | cut -f 1 -d ' '`
        tstart=`echo $line | cut -f 2 -d ' '`
        tend=`echo $line | cut -f 3 -d ' '`
        tstrand=`echo $line | cut -f 5 -d ' '`
        tcomplete=`echo $line | cut -f 6 -d ' '` 

        if [[ $tcomplete != "True" ]]
        then
            continue
        fi

        tstart=$((tstart - 100000))
        tend=$((tend + 100000))

        treg=${tcontig}:${tstart}-${tend}
        if [ $tcomplete = "True" ]
        then
            if [ $tstrand = "+" ]
            then
                ~/bin/samtools-1.14/samtools faidx ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.pat.chm13.cen.fa ${treg} >> ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.pat.chm13.cen.filt.fa
            else
                ~/bin/samtools-1.14/samtools faidx ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.pat.chm13.cen.fa -i ${treg} >> ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.pat.chm13.cen.filt.fa
            fi
        fi
    done < ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.pat.sd_info.txt
fi


echo -n > ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.mat.chm13.cen.filt.fa
mcount=`awk -F'\t' 'BEGIN {count=0} { if ($6 == "True") count++} END {print count}' ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.mat.sd_info.txt`
if [ $mcount = 1 ]
then

    while read line
    do
        tcontig=`echo $line | cut -f 1 -d ' '`
        tstart=`echo $line | cut -f 2 -d ' '`
        tend=`echo $line | cut -f 3 -d ' '`
        tstrand=`echo $line | cut -f 5 -d ' '`
        tcomplete=`echo $line | cut -f 6 -d ' '`

        if [[ $tcomplete != "True" ]]
        then
            continue
        fi

        tstart=$((tstart - 100000))
        tend=$((tend + 100000))

        treg=${tcontig}:${tstart}-${tend}
        if [ $tcomplete = "True" ]
        then
            if [ $tstrand = "+" ]
            then
                ~/bin/samtools-1.14/samtools faidx ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.mat.chm13.cen.fa ${treg} >> ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.mat.chm13.cen.filt.fa
            else
                ~/bin/samtools-1.14/samtools faidx ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.mat.chm13.cen.fa -i ${treg} >> ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.mat.chm13.cen.filt.fa
            fi
        fi
    done < ${OUT_DIR}/${SAMPLE}/${tcen_chr}/${SAMPLE}.mat.sd_info.txt
fi



