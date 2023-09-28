#! /usr/bin/env python3


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))

def kmer_check(qid, qseq):

    for i in range(len(qseq) - kmer_size + 1):

        kmer_is_inv = False
        if is_inv_contig and i >= inv_start and i <= inv_end: kmer_is_inv = True
        
        kmer = qseq[i:(i + kmer_size)]
        rkmer = reverse_complement(kmer)

        if kmer in kmer2info:
            kmer2info[kmer].append(f'{sample1},{platform},{contig_len},{i},+,{is_inv_contig},{kmer_is_inv}')

        if rkmer in kmer2info:
            kmer2info[rkmer].append(f'{sample1},{platform},{contig_len},{i},-,{is_inv_contig},{kmer_is_inv}')


import glob, os, sys
import statistics

kmer_file = sys.argv[1]
contig_dir = sys.argv[2]
output_file = sys.argv[3]

kmer_size = None

kmer2info = {}
with open(kmer_file, 'r') as hin:
    for line in hin:
        kmer2info[line.rstrip('\n')] = []
        kmer_size = len(line.rstrip('\n'))

all_fasta_files = glob.glob(contig_dir + "/*/*.chm13.cen.filt.fa")
for ffile in sorted(all_fasta_files):
    sample1 = os.path.basename(ffile).replace(".chm13.cen.filt.fa", '')
    sample2 = sample1.replace(".pat", '').replace(".mat", '')
    platform = os.path.basename(os.path.dirname(ffile))
    sd_file = ffile.replace(".chm13.cen.filt.fa", '') + ".sd_info.txt"

    is_inv_contig, inv_start, inv_end = False, None, None
    contig_len, strand = None, None
    with open(sd_file, 'r') as hin:
        F = hin.readline().rstrip('\n').split('\t')
        contig_len = int(F[2]) - int(F[1])
        strand = F[4]
        if F[6] == "True":
            is_inv_contig = True
            if strand == '+': 
                inv_start, inv_end = int(eval(F[9])[0][0]) - int(F[1]), int(eval(F[9])[0][1]) - int(F[1])
            else:
                inv_start, inv_end = int(F[2]) - int(eval(F[8])[0][1]), int(F[2]) - int(eval(F[8])[0][0])

    print(sample1, platform, contig_len, strand, is_inv_contig, inv_start, inv_end)

    temp_rid = None
    temp_rseq = None
    with open(ffile, 'r') as hin:

        for line in hin:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if temp_rid is not None:
                    kmer_check(temp_rid, temp_rseq)

                temp_rid = line.lstrip('>')
                temp_rseq = ''
            else:
                temp_rseq = temp_rseq + line

        if temp_rid is not None:
            kmer_check(temp_rid, temp_rseq)

        
with open(output_file, 'w') as hout:
    for kmer in kmer2info:

        rel_pos, num, num_kmer_is_inv, num_is_inv_contig = [], 0, 0, 0
        for FFF in kmer2info[kmer]:
            contig_name, platform, contig_len, pos, strand, is_inv_contig, kmer_is_inv = FFF.split(',')
            rel_pos.append(float(pos) / float(contig_len))
            num = num + 1
            if kmer_is_inv == "True":
                num_kmer_is_inv = num_kmer_is_inv + 1
            if is_inv_contig == "True":
                num_is_inv_contig = num_is_inv_contig + 1

        all_kmer_is_inv = False
        if num == num_kmer_is_inv:
            all_kmer_is_inv = True

        all_is_inv_contig = False
        if num == num_is_inv_contig:
            all_is_inv_contig = True

        rel_pos_median = statistics.median(rel_pos)
        rel_pos_var = statistics.variance(rel_pos) if len(rel_pos) > 1 else None

        tinfo = ';'.join(kmer2info[kmer])
        print(f'{kmer}\t{rel_pos_median}\t{rel_pos_var}\t{num}\t{num_kmer_is_inv}\t{num_is_inv_contig}\t{all_kmer_is_inv}\t{all_is_inv_contig}\t{tinfo}', file = hout)



