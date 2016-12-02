#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import math

a = []

for i in range(0,51):
    a.append(i)

PMV_D = {
        80: 0 ,
        70: 0 ,
        60: 0 ,
        40: 0 ,
        20: 0 ,
        10: 0 ,
         5: 0 ,
         1: 0
    }

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
stat_handle = open("Fin_output_SNP.txt", "r")

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

stat_handle.close()

## write SNP statistic results according PMV and MAF
PMV_out = open("test_PMV_SNP.txt", 'w')
PMV_out.write("MAF_value")

sort_v = sorted(list(PMV_D.keys()))
for v in sort_v:
    PMV_out.write("\t{0}".format(v))
PMV_out.write("\n")

for i in sorted(PMV_D[1]):
    PMV_out.write("{0}".format(i * 0.01))
    for j in sort_v:
        PMV_out.write("\t{0}".format( PMV_D[j][i] ))
    PMV_out.write("\n")

PMV_out.close()

## write SNP number of different maf thresholds
handle_maf = open("maf_stats.txt", 'w')
handle_maf.write("MAF\tSNP_NUM\n")
for i in sorted(maf_d.items(), key = lambda x: x[0]):
    handle_maf.write("{0}\t{1}\n".format(i[0]*0.01, i[1]))
handle_maf.close()

## write SNP number of different Herterozygosity
handle_het = open("het_stats_v.txt", 'w')
handle_het.write("HET\tSNP_NUM\n")

#for i in sorted(het_d.items(), key = lambda x: x[0]):   
#    handle_het.write("{0}\t{1}\n".format(i[0]*0.01, i[1]))

for i in sorted(list(het_d.keys())):
    if i % 2 == 0 and i < 100:      
       handle_het.write("{0}\t{1}\n".format( i * 0.01, het_d[i] + het_d[i+1]))

handle_het.close()

## End
