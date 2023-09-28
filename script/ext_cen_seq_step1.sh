#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log -o log
#$ -l s_vmem=4G -pe def_slot 8

set -xe

module load /usr/local/package/modulefiles/singularity/3.7.0

TASK_ID=$SGE_TASK_ID
PLATFORM=$1
# TASK_ID=$2

OUT_DIR=../output/cen_seq/${PLATFORM}
LIST_FILE=../db/${PLATFORM}_list.txt


SAMPLE=`head -n ${TASK_ID} ${LIST_FILE} | tail -n 1 | cut -f 1` 
PAT_PATH=`head -n ${TASK_ID} ${LIST_FILE} | tail -n 1 | cut -f 2`
MAT_PATH=`head -n ${TASK_ID} ${LIST_FILE} | tail -n 1 | cut -f 3`

if [ ! -d ${OUT_DIR}/${SAMPLE} ]
then
    mkdir -p ${OUT_DIR}/${SAMPLE}
fi


echo -e "${SAMPLE}\t${PAT_PATH}\t${MAT_PATH}"


~/.local/bin/aws s3 cp ${PAT_PATH} ${OUT_DIR}/${SAMPLE}/${SAMPLE}.pat.fa.gz
~/.local/bin/aws s3 cp ${MAT_PATH} ${OUT_DIR}/${SAMPLE}/${SAMPLE}.mat.fa.gz

~/bin/minimap2/minimap2 -t 8 -x asm5 ~/db/reference_genome/chm13.draft_v1.1.fasta ${OUT_DIR}/${SAMPLE}/${SAMPLE}.pat.fa.gz > ${OUT_DIR}/${SAMPLE}/${SAMPLE}.pat.chm13.paf
~/bin/minimap2/minimap2 -t 8 -x asm5 ~/db/reference_genome/chm13.draft_v1.1.fasta ${OUT_DIR}/${SAMPLE}/${SAMPLE}.mat.fa.gz > ${OUT_DIR}/${SAMPLE}/${SAMPLE}.mat.chm13.paf


