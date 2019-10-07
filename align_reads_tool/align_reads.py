#!/usr/bin/env python
# -*- coding: UTF-8

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from collections import defaultdict

d_ref = {}
ref_len = {}
pos_d = defaultdict(dict)

k_ref = {}
cut_len    = 22
search_len = int(cut_len/2) + 1

for record in SeqIO.parse("ref.fa", "fasta"):
    s = record.seq
    ref_len[record.id] = len(s)
    for i in range(0, len(s), cut_len):
        d_ref[str(s[i:i + cut_len])] = record.id + ':' + str(i)
        k_ref[record.id + ':' + str(i)] = str(s[i:i + cut_len])

p_set = set(d_ref.values())
for i in d_ref:
    str_K = d_ref[i]
    start_p = int(str_K.split(':')[1])
    ref_name = str_K.split(':')[0]

    region_s = ""
    if ref_name + ':' + str(start_p-1 * cut_len) in p_set and ref_name + ':' + str(start_p+1* cut_len) in p_set:
        region_s = k_ref[ref_name + ':' + str(start_p-1 * cut_len)] + i + k_ref[ref_name + ':' + str(start_p+1* cut_len)]
        for j in range(0, len(region_s), search_len): 
            pos_d[str_K][start_p - 1 * cut_len + j] = region_s[j:j+search_len]
    elif ref_name + ':' + str(start_p-1 * cut_len) not in p_set:
        region_s = i + k_ref[ref_name + ':' + str(start_p+1 * cut_len)]
        for j in range(0, len(region_s), search_len): 
            pos_d[str_K][start_p + j] = region_s[j:j+search_len]
    elif ref_name + ':' + str(start_p+1* cut_len) not in p_set:
        region_s = k_ref[ref_name + ':' + str(start_p-1* cut_len)] + i
        for j in range(0, len(region_s), search_len): 
            pos_d[str_K][start_p - 1 * cut_len + j] = region_s[j:j+search_len]

def rev_comp(s):
    r_s = ""

    d = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}
    for i in s:
        r_s = d[i] + r_s

    return r_s

for i, seq, qual in FastqGeneralIterator( open( "test1.fastq" ) ):
    rc_seq = rev_comp(str(seq))
    for j in range(0, len(seq) - cut_len):
        x = str(seq[j:j+cut_len])
        r_x = rev_comp(x)

        ref_p = ""
        if x in d_ref :
            ref_p = d_ref[x]
        elif r_x in d_ref:
            ref_p = d_ref[r_x]
        
        if ref_p != "":
            long_fragments = pos_d[ref_p]
            tmp_list = [] 
            for f,s in sorted(long_fragments.items()):
                sub_seq = pos_d[ref_p][f]
                if sub_seq in seq or sub_seq in rc_seq:
                    tmp_list.append(f)
            print("{0}\t{1}\t{2}\t{3}".format(i, ref_p, tmp_list[0], sub_seq))
            break
