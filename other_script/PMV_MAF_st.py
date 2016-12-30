#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import math
import argparse

### the output files is ready for Plot_e.R ###

def data_ready(input_f):
    a = list(range(0,51))
    PMV_D = {
        80: 0 ,
        70: 0 ,
        60: 0 ,
        40: 0 ,
        20: 0 ,
        10: 0 ,
         5: 0 ,
         1: 0}

    for i in PMV_D:
        d = {}
        for j in a:
            d[j] = 0
        PMV_D[i] = d

    maf_d = {}
    for num in a:
        maf_d[num] = 0

    het_d = {}
    for i in range(0,101):
        het_d[i] = 0

    ## Reade SNP filteration result and get maf Het observation and Percent missing value , store them in corresponding dict.

    with open(input_f, 'r') as stat_handle:
        for i in stat_handle:
            if i.startswith("Chr"):
                continue
            else:
                sp_line = i.rstrip().split("\t")
                maf = int(float(sp_line[7]) * 100)   # could use math.floor in replace of int() function
                maf_d[maf] += 1
                pmv = float(sp_line[6])
                het = int( float(sp_line[9]) * 100)
                het_d[het] += 1
                for k in PMV_D:
                    if pmv <= k:
                        PMV_D[k][maf] += 1

    return (PMV_D,  maf_d,  het_d)
if __name__ == "__main__":

    USAGE = "python %s.py -i vcf"%("this script ")

    parser = argparse.ArgumentParser(description = USAGE)
    parser.add_argument("-i", "--input", action = "store", required = True, help = "the vcffile containing variants, including SNP and INDEL")
    parser.add_argument("-o", "--output", action = "store", default = "test_PMV_SNP", help = "the output file containing variants, including SNP and INDEL")

    args = parser.parse_args()

    stat_1, stat_2, stat_3 = data_ready(args.input)
    ## write SNP statistic results according PMV and MAF
    PMV_out = open(args.output + ".txt", 'w')
    PMV_out.write("MAF_value")

    sort_v = sorted(list(stat_1.keys()))
    for v in sort_v:
        PMV_out.write("\tP{0}".format(v))
    PMV_out.write("\n")

    for i in sorted(stat_1[1]):
        PMV_out.write("{0}".format(i * 0.01))
        for j in sort_v:
            PMV_out.write("\t{0}".format( stat_1[j][i] ))
        PMV_out.write("\n")

    PMV_out.close()

    ## write SNP number of different maf thresholds
    handle_maf = open("maf_stats.txt", 'w')
    handle_maf.write("MAF\tSNP_NUM\n")
    for i in sorted(stat_2.items(), key = lambda x: x[0]):
        handle_maf.write("{0}\t{1}\n".format(i[0]*0.01, i[1]))
    handle_maf.close()

    ## write SNP number of different Herterozygosity
    handle_het = open("het_stats_v.txt", 'w')
    handle_het.write("HET\tSNP_NUM\n")

    for i in sorted(list(stat_3.keys())):
        if i % 2 == 0 and i < 98:      
            handle_het.write("{0}\t{1}\n".format( i * 0.01, stat_3[i] + stat_3[i+1]))
        elif i == 98:
            handle_het.write("{0}\t{1}\n".format( i * 0.01, stat_3[i] + stat_3[i+1] + stat_3[i+2]))
        elif i == 100:
            handle_het.write("{0}\t{1}\n".format( i * 0.01, 0))
    handle_het.close()

## End
