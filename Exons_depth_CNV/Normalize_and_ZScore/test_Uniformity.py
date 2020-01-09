#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from glob import glob
import numpy as np
from collections import defaultdict
import copy
import os

coverageDir = "coverage_files"
genderF = open( 'samples.xls'  ,  'r'  )


female_samples =  []

male_samples   =  []

for line in genderF:
    x = line.rstrip().split('\t')

    if x[1] == 'F':
        female_samples.append(x[0])

    else:
        male_samples.append(x[0])

genderF.close()

## exon length of gene
exonLen_dict = {}
Exon_file = open(  'exons_length.txt'  ,  'r'  )
for line in Exon_file:
    k = line.rstrip().split('\t')
    exonLen_dict[k[0]]  =  int( k[1] )
Exon_file.close()

covfiles = glob(coverageDir + '/' + "*base.coverage")


d_Uniform = defaultdict(list)
for f in covfiles:
    fname = os.path.basename(f)
    sampleName = fname.split('.')[0]

    if sampleName in male_samples:
        dd_depth = copy.deepcopy(defaultdict(list))
        fs = open(f, 'r')
        for line in fs:
            if not line.startswith('chrom'):
                linecol = line.rstrip().split('\t')
                dd_depth[linecol[2]].append(int( linecol[3] ) )
                
        fs.close()        
        for exon in dd_depth:
            exonLen = exonLen_dict[exon]
            depS = np.array(dd_depth[exon])
            d_avg = np.mean(depS)
        
            above_AvgNum = 0
            for d in depS:
                if d >= d_avg * 0.5:
                    above_AvgNum += 1
                
            d_Uniform[exon].append(above_AvgNum * 100.0/exonLen)    
            #print('{0}\t{1}\t{2}\t{3}\t{4:.2f}'.format(f, exon, exonLen, above_AvgNum, above_AvgNum * 100.0/exonLen))

cutoff = 90
high_Exons = set()
for i in d_Uniform:
    
    U = d_Uniform[i]
    #print('%s\t%s' % ( i, d_Uniform[i] ) )    

    if all([ x > cutoff for x in U]):
        high_Exons.add(i)
        print("{0}\t{1}\t{2:.2f}".format(i, len(U), np.mean(U)))
        
a_out = open("bad_uniformity_Exons.xls", 'w')
for i in exonLen_dict:
    if i not in high_Exons:
        a_out.write("{0}\n".format(i))
        
a_out.close()
