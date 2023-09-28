#! /usr/bin/env python3

import sys

input_file = sys.argv[1]

cen_margin = 5000
cen_thres = 90.0


class SD_info:

    def __init__(self, contig_name):
        self.contig_name = contig_name
        self.cen_start_tmp = None
        self.cen_end_tmp = None
        self.cen_start_plus_tmp = None
        self.cen_end_plus_tmp = None
        self.cen_start_minus_tmp = None
        self.cen_end_minus_tmp = None

        self.cen_region_list = []
        self.cen_region_plus_list = []
        self.cen_region_minus_list = []

        self.contig_len = None
        self.is_complete = None
        self.is_inversion = None

        self.cen_start = None
        self.cen_end = None
        self.strand = None


    def flush(self):

        if self.cen_end_tmp is not None and abs(self.contig_len - self.cen_end_tmp) < cen_margin:
            self.cen_region_list.append((self.cen_start_tmp, self.cen_end_tmp))

        if self.cen_end_plus_tmp is not None and abs(self.contig_len - self.cen_end_plus_tmp) < cen_margin:
            self.cen_region_plus_list.append((self.cen_start_plus_tmp, self.cen_end_plus_tmp))

        if self.cen_end_minus_tmp is not None and abs(self.contig_len - self.cen_end_minus_tmp) < cen_margin:
            self.cen_region_minus_list.append((self.cen_start_minus_tmp, self.cen_end_minus_tmp))


        max_reg_len = 0
        for reg in self.cen_region_list:
            if reg[1] - reg[0] > max_reg_len:
                self.cen_start = reg[0]
                self.cen_end = reg[1]
                max_reg_len = reg[1] - reg[0]

        self.is_complete = False
        if len(self.cen_region_list) == 1:
            if (max_reg_len > 500000 and self.cen_start > 5 * cen_margin and self.contig_len - self.cen_end > 5 * cen_margin): self.is_complete = True
                 
        
        if self.is_complete:

            is_start_match_plus, is_end_match_plus, is_start_match_minus, is_end_match_minus = False, False, False, False
            cen_size_plus, cen_size_minus = 0, 0

            for reg in self.cen_region_plus_list:
                if abs(self.cen_start - reg[0]) < cen_margin: is_start_match_plus = True
                if abs(self.cen_end - reg[1]) < cen_margin: is_end_match_plus = True
                cen_size_plus = cen_size_plus + reg[1] - reg[0]
            
            for reg in self.cen_region_minus_list:
                if abs(self.cen_start - reg[0]) < cen_margin: is_start_match_minus = True
                if abs(self.cen_end - reg[1]) < cen_margin: is_end_match_minus = True
                cen_size_minus = cen_size_minus + reg[1] - reg[0]

            if is_start_match_plus and is_end_match_plus: self.strand = '+'
            if is_start_match_minus and is_end_match_minus: self.strand = '-'

            self.is_inversion = False
            if self.strand == '+' and cen_size_minus > 500000: self.is_inversion = True
            if self.strand == '-' and cen_size_plus > 500000: self.is_inversion = True           

        print(f'{self.contig_name}\t{self.cen_start}\t{self.cen_end}\t{self.contig_len}\t{self.strand}\t{self.is_complete}\t{self.is_inversion}\t{str(self.cen_region_list)}\t{str(self.cen_region_plus_list)}\t{str(self.cen_region_minus_list)}')





sd_info = SD_info(None)

with open(input_file, 'r') as hin:

    for line in hin:
        F = line.rstrip('\n').split('\t')

        if F[0] != sd_info.contig_name:
            if sd_info.contig_name is not None:
                sd_info.flush()

            sd_info = SD_info(F[0])

        strand = '+' if not F[1].endswith("'") else '-'


        if float(F[4]) > cen_thres:

            if sd_info.cen_start_tmp is None: 
                sd_info.cen_start_tmp = int(F[2])
            if sd_info.cen_start_tmp is not None:
                sd_info.cen_end_tmp = int(F[3])

            if strand == '+':

                if sd_info.cen_start_plus_tmp is None:
                    sd_info.cen_start_plus_tmp = int(F[2])
                if sd_info.cen_start_plus_tmp is not None:
                    sd_info.cen_end_plus_tmp = int(F[3])
    
                if sd_info.cen_start_minus_tmp is not None:
                    if int(F[2]) - sd_info.cen_end_minus_tmp > cen_margin:
                        sd_info.cen_region_minus_list.append((sd_info.cen_start_minus_tmp, sd_info.cen_end_minus_tmp))
                        sd_info.cen_start_minus_tmp = None
                        sd_info.cen_end_minus_tmp = None

            else:

                if sd_info.cen_start_minus_tmp is None:
                    sd_info.cen_start_minus_tmp = int(F[2])
                if sd_info.cen_start_minus_tmp is not None:
                    sd_info.cen_end_minus_tmp = int(F[3])

                if sd_info.cen_start_plus_tmp is not None:
                    if int(F[2]) - sd_info.cen_end_plus_tmp > cen_margin:
                        sd_info.cen_region_plus_list.append((sd_info.cen_start_plus_tmp, sd_info.cen_end_plus_tmp))
                        sd_info.cen_start_plus_tmp = None
                        sd_info.cen_end_plus_tmp = None

 
        else:

            if sd_info.cen_start_tmp is not None:
                if  int(F[2]) - sd_info.cen_end_tmp > cen_margin: 
                    sd_info.cen_region_list.append((sd_info.cen_start_tmp, sd_info.cen_end_tmp))
                    sd_info.cen_start_tmp = None
                    sd_info.cen_end_tmp = None

            if sd_info.cen_start_minus_tmp is not None:
                if int(F[2]) - sd_info.cen_end_minus_tmp > cen_margin:
                    sd_info.cen_region_minus_list.append((sd_info.cen_start_minus_tmp, sd_info.cen_end_minus_tmp))
                    sd_info.cen_start_minus_tmp = None
                    sd_info.cen_end_minus_tmp = None

            if sd_info.cen_start_plus_tmp is not None:
                if int(F[2]) - sd_info.cen_end_plus_tmp > cen_margin:
                    sd_info.cen_region_plus_list.append((sd_info.cen_start_plus_tmp, sd_info.cen_end_plus_tmp))
                    sd_info.cen_start_plus_tmp = None
                    sd_info.cen_end_plus_tmp = None


   
        sd_info.contig_len = int(F[3])


    if sd_info is not None:

        sd_info.flush()


