#! /usr/bin/env python3

import sys, os
import networkx as nx


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))


   
def kmer_edge_check(fasta_file, G, kmer_list, kmer_size = 27): 


    def seq_kmer_edge_check(qid, qseq):

        prev_kmer = None
        prev_rkmer = None
        for i in range(len(qseq) - kmer_size + 1):
            kmer = qseq[i:(i + kmer_size)]
            rkmer = reverse_complement(kmer)

            if prev_kmer is not None and prev_rkmer is not None:
                if prev_kmer in kmer_list and kmer in kmer_list:
                    G.add_edge(prev_kmer, kmer)
                if prev_rkmer in kmer_list and rkmer in kmer_list:
                    G.add_edge(prev_rkmer, rkmer)

            prev_kmer = kmer
            prev_rkmer = rkmer


    temp_rid = None
    temp_rseq = None        
    with open(fasta_file, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if temp_rid is not None:
                    seq_kmer_edge_check(temp_rid, temp_rseq)

                temp_rid = line.lstrip('>')
                temp_rseq = ''
            else:
                temp_rseq = temp_rseq + line

        if temp_rid is not None:
            seq_kmer_edge_check(temp_rid, temp_rseq)



import glob, os, sys

kmer_file = sys.argv[1]
input_dir = sys.argv[2]
output_file = sys.argv[3]


all_fasta_files = glob.glob(input_dir + "/*/*.chm13.cen.filt.fa")

kmer_list = {}
kmer_size = None
with open(kmer_file, 'r') as hin:
    for line in hin:
        kmer = line.rstrip('\n')
        kmer_list[kmer] = 1
        kmer_size = len(kmer)

G = nx.Graph()
for kmer in kmer_list:
    rkmer = reverse_complement(kmer)
    if rkmer == kmer: continue
    if rkmer in kmer_list:
        G.add_edge( kmer, rkmer )


for fasta_file in sorted(all_fasta_files):
    hap = os.path.basename(fasta_file).replace(".chm13.cen.filt.fa", '')
    # if hap in blacklist_hap: continue

    kmer_edge_check(fasta_file, G, kmer_list, kmer_size)
    print(f'Read haplotype: {hap}')
    

hout = open(output_file, 'w') 
con = 0
for cG in nx.connected_components(G):
    temp_kmer = None
    for kmer in cG:
        temp_kmer = kmer

    print(f'{temp_kmer}', file = hout)

hout.close()


