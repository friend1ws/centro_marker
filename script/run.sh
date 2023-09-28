#! /bin/bash

qsub -t 1-12:1 ext_cen_seq.sh hifiasm_v0.19.5

qsub -t 1-10:1 ext_cen_seq.sh verkko_1.3.1 

qsub -t 1-47:1 ext_cen_seq.sh hifiasm_v0.14_raw


