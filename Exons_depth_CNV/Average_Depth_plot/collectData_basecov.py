#!/usr/bin/env python
# -coding: UTF-8 -*-

from glob import glob
import os
sfile = "samples.txt"

bedFile = "T_exons.bed"

d_base = {}
hb = open(  bedFile  ,  'r'  )
for line in hb:
    x = line.rstrip().split('\t')
    start = int(x[1])
    end   = int(x[2])

    for i in range(start , end + 1 ):
        d_base[x[0] + ":" + str(i)] = x[3]
hb.close()

handle_S = open(  sfile ,  'r'  )
for line in handle_S:
    sample = line.rstrip()
    outC = open("Exons_" + sample + "_cov.xls", 'w')

    all_baseCovs= glob('/Data/NGS_analysis/'+ sample + '/2.mapping/*base.coverage' )
    baseCovFile = ''
    for i in all_baseCovs:
        sssName = os.path.basename(i)
        if sssName.split('.')[0] != sample:
            baseCovFile = i
        else:
            continue 
    if baseCovFile == '':
        continue

    bcov_h = open( baseCovFile  ,  'r'  ) 
    bcov_h.readline()
    for line in bcov_h:
        #if not line.startswith("chrom")

        x = line.rstrip().split('\t')
        k =  x[0] + ':' + x[1]
        if k in d_base:
            outC.write("{0}\t{1}\t{2}\t{3}\n".format(x[0], x[1], x[3], d_base[k]))
    outC.close()
    bcov_h.close()



