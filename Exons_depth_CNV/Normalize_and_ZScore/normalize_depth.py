#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from glob import glob
import os


sampleList = []

covfiles = glob( "coverage_files/*.base.coverage" )
outDir   = 'standardized_Depth/'

a = open(  'samples_mean_depth.xls'  ,  'r'  )
for line in a:
    x = line.rstrip().split('\t')
    sampleList.append(x[0])

a.close()

avg_Dep = []

sampleD_dict = {}
for i in covfiles:
    sample = os.path.basename(i).split('.')[0]

    if sample in sampleList:
        X = []

        fh = open(i, 'r')
        for line in fh:
            if not line.startswith('chrom'):
                px = line.rstrip().split('\t')
                X.append(int(px[3]))

        fh.close()
        avg_Dep.append(sum(X) * 1.0/ len(X))
        sampleD_dict[sample]  =  sum(X) * 1.0/ len(X)
        print("{0}\t{1}".format(sample, sum(X) * 1.0/ len(X)))

MaxD = max(avg_Dep)

for i in covfiles:
    sample = os.path.basename(i).split('.')[0]

    if sample in sampleList:
        outA = open(outDir + sample + '.n_cov.txt' , 'w')
        fh = open(i, 'r')
        for line in fh:
            if not line.startswith('chrom'):
                px        = line.rstrip().split('\t')
                new_depth = int(px[3]) * MaxD / sampleD_dict[sample]
                outA.write("{0}\t{1}\t{2}\t{3}\n".format(px[0], px[1], px[2], new_depth) )
            else:
                outA.write(line)
                
        outA.close()
        fh.close()
                    


