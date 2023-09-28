#! /bin/bash

if [ ! -f chm13.draft_v1.1.cenAnnotation.bed ]
then
    aws s3 cp s3://human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v1.1.cenAnnotation.bed ./
fi

python3 proc_cen_list.py chm13.draft_v1.1.cenAnnotation.bed > cen_list.proc.bed


