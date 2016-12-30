#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from collections import defaultdict
from random import randint
import os, sys
import math
import argparse

def prepare_DataFor_Circ(SNP_filter, window, List_F, Spe_Abbr, bam_F, bed_F ):
    chromo = defaultdict(dict)

    CHR_N = {}
    with open(List_F, 'r') as CT:
        for i in CT:
            C = i.rstrip().split("\t")
            CHR_N[C[1]] = Spe_Abbr + C[0]

    maf_c = open("maf_norm.txt", 'w')
    pmv_c = open("pmv_norm.txt", 'w')

    maf_d = defaultdict(list)
    pmv_d = defaultdict(list)

    with open(SNP_filter , 'r') as SNP_summary:
        for i in SNP_summary:
            if i.startswith("Chr"):
                continue
            else:
                s = i.rstrip().split("\t")
                pos = int(s[1])
                chr_r = int( pos / window) * window
                maf_tmp = 1/(float(s[7]) + 0.001)
                mp_key = int(pos/(window/10)) * window/10
                if s[0] in CHR_N:
                    #maf_c.write("{0} {1} {2} {3}\n".format(CHR_N[s[0]], pos, pos , math.log(maf_tmp, 2)))
                    #pmv_c.write("{0} {1} {2} {3}\n".format(CHR_N[s[0]], pos, pos , 100/(float(s[6]) + 1)))
                    cn = CHR_N[s[0]]
                    maf_d[cn+"_"+str(mp_key)].append(math.log(maf_tmp,2))
                    pmv_d[cn+"_"+str(mp_key)].append(100/(float(s[6]) + 1))
                if chr_r in chromo[s[0]]:
                    chromo[s[0]][chr_r] += 1
                else:
                    chromo[s[0]][chr_r] = 1

    for i in sorted(maf_d.items(), key = lambda x: x[0].split("_")[0]):
        cs = i[0].split("_")
        chr = cs[0]
        start = cs[1]
        end = int(start) + window/10 -1
        aver = sum(i[1])/len(i[1])
        maf_c.write("{0} {1} {2} {3}\n".format(chr, start, end , aver))

    for i in sorted(pmv_d.items(), key = lambda x: x[0].split("_")[0]):
        cs = i[0].split("_")
        chr = cs[0]
        start = cs[1]
        end = int(start) + window/10 -1
        aver = sum(i[1])/len(i[1])
        pmv_c.write("{0} {1} {2} {3}\n".format(chr, start, end , aver))

    maf_c.close()
    pmv_c.close()
    chr_info = {}

    #  get chromosome info and create genome file
    chr_reading = os.popen("samtools view -H "+ bam_F)
    for i in chr_reading.read().split("\n"):
        if i.startswith("@SQ"):
            info = i.split("\t")
            ch_n = info[1].split(":")[1]
            ch_l = info[2].split(":")[1]
            chr_info[ch_n] = ch_l

    colors = [["vvlred","vlred","lred","red","dred","vdred","vvdred"],
          ["vvlpred", "vlpred", "lpred", "pred", "dpred", "vdpred", "vvdpred"],
          ["vvlgreen", "vlgreen" ,"lgreen", "green", "dgreen", "vdgreen", "vvdgreen"],
          ["vvlpgreen", "vlpgreen", "lpgreen", "pgreen", "dpgreen", "vdpgreen", "vvdpgreen"],
          ["vvlblue", "vlblue", "lblue", "blue", "dblue","vdblue","vvdblue"],
          ["vvlpblue", "vlpblue", "lpblue", "pblue", "dpblue", "vdpblue", "vvdpblue"],
          ["vvlpurple", "vlpurple", "lpurple", "purple", "dpurple", "vdpurple", "vvdpurple"],
          ["vvlppurple", "vlppurple", "lppurple", "ppurple", "dppurple", "vdppurple", "vvdppurple"],
          ["vvlorange", "vlorange", "lorange", "orange", "dorange", "vdorange", "vvdorange"],
          ["vvlporange", "vlporange", "lporange", "porange", "dporange", "vdporange", "vvdporange"],
          ["vvlyellow", "vlyellow", "lyellow", "yellow", "dyellow", "vdyellow", "vvdyellow"],
          ["vvlpyellow", "vlpyellow", "lpyellow", "pyellow", "dpyellow", "vdpyellow", "vvdpyellow"]]

    genome_C = open("SPE_GENOME.txt", 'w')
    sort_chromosome = sorted(CHR_N.items(), key = lambda x: int(x[1].lstrip(Spe_Abbr)))
    for n, i in enumerate(sort_chromosome):
        genome_C.write("chr - {0} {1} 0 {2} {3}\n".format(i[1], n+1, chr_info[i[0]], colors[randint(0,11)][randint(0,6)]))
    genome_C.close()

    tag_c = defaultdict(dict)
    for i in chr_info:
        CL = int(chr_info[i])   # CL is chromosome Length
        if CL >= window:
            for j in range(0, CL + 1, window):
                tag_c[i][j] = 0

    with open(bed_F , 'r') as Bed_tags:
        for i in Bed_tags:
            sp = i.rstrip().split("\t", 3)
            if len(sp) > 3:
                key_pos = (int(sp[1])//window) * window
                site = key_pos + window
                start = int(sp[1])
                end = int(sp[2])
                if 2 * site >=  start + end and key_pos in tag_c[sp[0]]:
                    tag_c[sp[0]][key_pos] += 1
                elif site in tag_c[sp[0]] and key_pos+window in tag_c[sp[0]]:
                    tag_c[sp[0]][key_pos + window] += 1

    print ("Tag number analysis starts:")

    tag_OUT = open("tag_density.txt", 'w')

    for i in CHR_N:
        for j in sorted(tag_c[i]):
        #tag_norm = math.log(tag_c[i][j]/1000 + 1, 10)
            tag_norm = (tag_c[i][j]/1000 + 1)/ 10**1.5
            tag_OUT.write("{0} {1} {2} {3}\n".format(CHR_N[i], j, j+window-1, math.log(tag_norm +1, 2)))

    tag_OUT.close()

    SNP_OUT = open("snp_density.txt", 'w')
    for i in CHR_N:
        CL = int(chr_info[i])   # CL is chromosome Length
        if CL >= window:
            for j in range(0, CL, window):
                if j in chromo[i]:
                    end = j + window
                    snp_norm = math.log(chromo[i][j] + 1, 10)
                    SNP_OUT.write("{0}\t{1}\t{2}\t{3}\n".format(CHR_N[i], j, end-1 , snp_norm))

    SNP_OUT.close()

if __name__ == "__main__":

    USAGE = "python %s.py -i vcf"%("this script")

    parser = argparse.ArgumentParser(description = USAGE)
    parser.add_argument("-i", "--input", action = "store", required = True, help = "the tab delimited file containing SNP position")
    parser.add_argument("-w", "--win", type = int, default = 100000, help = "the window of chromosome region")
    parser.add_argument("-L", "--list_c", action ="store",  help = "the list of chromosome names")
    parser.add_argument("-A", "--Abbr", action = "store", default = "Chr", help = "Abbreviation of species.")
    parser.add_argument("-b", "--bamfile", action = "store", required=True, help = "the bamfile")
    parser.add_argument("-d", "--bedfile", action = "store", required=True, help = "the bedfile of merged bam")

    args = parser.parse_args()
    if (not os.path.exists(args.bamfile)) or (not os.path.exists(args.bedfile)):
        print ("""bam and bed file required, please the path of these two file, 
    if you have bam file but not bed file, use the following command to generated it.
         /usr/bin/bamToBed -i some.bam > some.bed""")
        sys.exit("Quit for some file missing.")
    prepare_DataFor_Circ(args.input, args.win, args.list_c, args.Abbr, args.bamfile, args.bedfile )
    #CHR_Tran = "/home/yangjy/16T/GENOME/GEN160781/GEN160781/Ref_Genome/List.txt"
