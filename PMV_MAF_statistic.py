#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import math

# use the following awk sentence to validate the results
# awk -F '\t' '{if ($7<=70 && ($8 <0.02 && $8 >=0.01)) print}' Fin_output_SNP.txt | wc -l  
# awk -F '\t' '{if ($7<=70 && ($8 <0.01 && $8 >=0.00)) print}' Fin_output_SNP.txt | wc -l        
# awk -F '\t' '{if ($7<=80 && ($8 <0.01 && $8 >=0.00)) print}' Fin_output_SNP.txt | wc -l  

a = []

for i in range(0,51):
    a.append(i)
    #a.append(str(round(0.01 * i, 2)))

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

# MAF_D = dict(zip(a, [0] * 50))
# for i in PMV_D:
#     PMV_D[i] = copy.deepcopy(MAF_D)    # should import the copy module

stat_handle = open("Fin_output_SNP.txt", "r")

for i in stat_handle:
    if i.startswith("Chr"):
        continue
    else:
        sp_line = i.rstrip().split("\t")
        maf = int(float(sp_line[7]) * 100)   # could use math.floor in replace of int() function
        
        pmv = float(sp_line[6])

        for k in PMV_D:
            if pmv <= k:
                PMV_D[k][maf] += 1

stat_handle.close()

PMV_out = open("test_PMV_SNP.txt", 'w')
PMV_out.write("PMV_value")
for freq in a:
    PMV_out.write("\t{0}".format(freq))
PMV_out.write("\n")

for i in sorted(list(PMV_D.keys())):
    PMV_out.write("{0}".format(i))
    for j in sorted(PMV_D[i]):
        PMV_out.write("\t{0}".format( PMV_D[i][j] ))
    PMV_out.write("\n")

PMV_out.close()
