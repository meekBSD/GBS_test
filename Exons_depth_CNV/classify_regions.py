#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys

exon_bases = set()

gtfF = open(sys.argv[1], 'r' )

for line in gtfF:
    x = line.rstrip().split('\t')
    if x[2] == 'exon':
        for j in range(int(x[3]), int(x[4]) +1):
            exon_bases.add(j)

gtfF.close()

depFile = open(sys.argv[2], 'r')

for i in depFile:
    x =i.rstrip().split('\t')[1]
    if int(x) in exon_bases:
        print(i.rstrip() + '\texon' )
    else:
        print(i.rstrip() + '\tintron')    

depFile.close()

