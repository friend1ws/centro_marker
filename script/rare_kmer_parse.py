#! /usr/bin/env python3

import sys, os


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))


   
def gather_kmer(fasta_file, kmer_size = 27): 

    kmer2count = {}

    def rare_kmer_check(qid, qseq):

        for i in range(len(qseq) - kmer_size + 1):
            kmer = qseq[i:(i + kmer_size)]
            rkmer = reverse_complement(kmer)

            # if kmer == "TATTGTGTGTATTCAACTCACAGAGTTGAAC" or rkmer == "TATTGTGTGTATTCAACTCACAGAGTTGAAC":
            #     print(f'{qid} {i} {kmer}')

            if kmer not in kmer2count: kmer2count[kmer] = 0
            kmer2count[kmer] = kmer2count[kmer] + 1

            if rkmer not in kmer2count: kmer2count[rkmer] = 0
            kmer2count[rkmer] = kmer2count[rkmer] + 1


    temp_rid = None
    temp_rseq = None        
    with open(fasta_file, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if temp_rid is not None:
                    rare_kmer_check(temp_rid, temp_rseq)

                temp_rid = line.lstrip('>')
                temp_rseq = ''
            else:
                temp_rseq = temp_rseq + line

        if temp_rid is not None:
            rare_kmer_check(temp_rid, temp_rseq)

    return(kmer2count)


import glob, os, sys

input_dir = sys.argv[1]
output_file = sys.argv[2]
kmer_size = sys.argv[3]

all_fasta_files = glob.glob(input_dir + "/*/*.chm13.cen.filt.fa")
black_list_kmer = {}
unique_kmer_list_raw = {}
for fasta_file in sorted(all_fasta_files):
    hap = os.path.basename(fasta_file).replace(".chm13.cen.filt.fa", '')

    temp_kmer2count = gather_kmer(fasta_file, int(kmer_size))
    print(f'Read haplotype: {hap}, kmer count: {len(temp_kmer2count)}') 
    
    
    for kmer in temp_kmer2count:
        if temp_kmer2count[kmer] > 1: black_list_kmer[kmer] = 1

    # import pdb; pdb.set_trace()

    for kmer in temp_kmer2count:    
        if kmer not in black_list_kmer: unique_kmer_list_raw[kmer] = 1

    print(f'Current unique kmer count: {len(unique_kmer_list_raw)}')

unique_kmer_list = {}
for kmer in unique_kmer_list_raw:
    if kmer not in black_list_kmer: unique_kmer_list[kmer] = 1

print(f'Final unique kmer count: {len(unique_kmer_list)}')
with open(output_file, 'w') as hout:
    for kmer in unique_kmer_list:
        print(kmer, file = hout)


