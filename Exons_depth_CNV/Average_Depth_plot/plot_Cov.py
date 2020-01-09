#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from collections import defaultdict


geneNames =  []
h1 = open( 'gene_regions.txt'  ,  'r'  )

for line in h1:
    x = line.rstrip()
    geneNames.append(x)
h1.close()
    
exonCovFiles = glob('Exons_*_cov.xls') 

for g in geneNames:
    d1 = defaultdict(list)
    for f in exonCovFiles:
        sampleID = f.split('_')[1]
        fh = open(f, 'r')
        for line in fh:
            x = line.rstrip().split('\t')
            if x[3] == g:
                d1[sampleID].append('\t'.join(x[1:3]))
    
    fh.close()
    
    plt.rcParams['figure.figsize'] = (15, 4.5)
    colors = ["#2F4F4F", "#CD7F32", "#527F76", "#4E2F2F", "#32CD32", "#32CD99", "#7F00FF", "#DB7093", "#FF7F00", "007FFF"]
    
    ci = 0
    for i in d1:
        pv = d1[i]
        gd = {}
        for j in pv:
            xj = j.split("\t")
            gd[int(xj[0])] = int(xj[1])
        genome_pos = [int(x.split('\t')[0]) for x in pv]
        
        ## store X values
        e_start = min(genome_pos)
        e_end   = max(genome_pos)
        X = np.arange(e_start, e_end + 1)
        
        if e_end - e_start > 650000:
            break
        
        ## store Y values
        Y = []
        for pos in X:
            if pos not in genome_pos:
                Y.append(0)
            else:
                Y.append(gd[pos])    

        plt.plot(X,Y,label=i,linewidth=0.24,color= colors[ci],marker='o',  markerfacecolor='blue',markersize=0.4) 
        ci += 1

    plt.xlabel('Genome Coordinate') 
    plt.ylabel('Sequencing Depth') 
    plt.title('%s\nCoverage'% g) 
    plt.legend(bbox_to_anchor=(1.05, 0), loc=3, borderaxespad=0) 
    #plt.legend() 
    #plt.show()
    plt.savefig("out_%s.png" % g, bbox_inches="tight")
    plt.close() 


