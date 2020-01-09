#!/usr/bin/env python
# -*- coding: UTF-8 -*-


Elist = []

k_h = open( "../exons_Uniformity.txt", 'r' )
for line in k_h:
    if not line.startswith("Exon_ID"):
        mx = line.rstrip().split('\t')
        Elist.append(mx[0])

k_h.close()

outZ = open( 'select_Z_score.xls', 'w' )
a = open( 'test_Z-Score.xls', 'r')
for line in a :
    x = line.rstrip().split('\t')
    if not line.startswith("target") and x[0] in Elist:
        if all([i != 'NA' for i in x[1:30] ]):
            z_List = [float(i) for i in x[1:30]]

            if any([abs(z) >= 3 for z in z_List]) and float(x[31]) > 15:
                outZ.write(line)

outZ.close()
a.close()
        
## filter Exons by Z_score percentage for each Exon region
        
Elen = dict()

preStat = open(  '../exons_length.txt' ,  'r'  )
for line in preStat:
    x = line.rstrip().split('\t')
    Elen[x[0]] = int(x[1])

preStat.close()

dz = {}

testO = open( 'select_Z_score.xls',  'r' )
for line in testO:
    kx = line.rstrip().split('\t')

    dz[kx[0]] = dz.get(kx[0] , 0) + 1
testO.close()


for i in dz:
    print("{0}\t{1}\t{2}\t{3}".format(i, dz[i], Elen[i], dz[i] * 100.0/Elen[i]))
