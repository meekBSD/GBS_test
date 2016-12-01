#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from collections import defaultdict
import os

# define windows_size  and slide region
# need bamfile for SQ header for chromosomes length

window_size = 100000
step = 10000


SNP_file = open("test_FreeBayes_ninty_SNP.txt", 'r')
chromo = defaultdict(dict)

for i in SNP_file:
    if i.startswith("Chr"):
        continue
    else:
        splits = i.rstrip().split("\t")
        chr_name = splits[0]
        Position = splits[1]
        steps_in_win = window_size/step
        chr_win_start = (int(int(Position)//step) - steps_in_win +1) * step
        for win in range(chr_win_start, chr_win_start + window_size, step): 
            if win >= 0 and win in chromo[chr_name]:
                chromo[chr_name][win] += 1
            elif win >= 0:
                chromo[chr_name][win] = 1
SNP_file.close()

'''
print ("The number and names of chromosomes in this specie is:")
print (len(chromo))


for i in chromo:
    print (i)
    print (len(chromo[i]))
    for j in chromo[i]:
        pass
        #print (j)
'''

chr_info = {}

chr_reading = os.popen("samtools view -H Aln_PE_17_sorted.bam")
for i in chr_reading.read().split("\n"):
    if i.startswith("@SQ"):
        info = i.split("\t")
        ch_n = info[1].split(":")[1]
        ch_l = info[2].split(":")[1]
        chr_info[ch_n] = ch_l

for i in chr_info:
    chromosome_length = int(chr_info[i])
    if chromosome_length >= window_size:
        for j in range(0, int(chr_info[i]), step):
            if j + window_size <= chromosome_length + step:
                locus = j + window_size/2
                print ("{0}\t{1}-{2}\t{3}\t{4}".format(i, j, j+window_size ,locus,  chromo[i].get(j, 0)))
