#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from collections import defaultdict
import numpy as np
import os, re

covFiles = [ ]

###  ll -d /mnt/disk01/project/*

for line in open('targets_samples.txt', 'r'):
    sX = line.rstrip().split('\t')
     
    ctdir = sX[0]+"/outS.ctdna_fraction/"
    baseF = sX[0] + "/VD.mapping/"  +sX[1] +".base.coverage"
    if not os.path.exists(ctdir) or len(os.listdir(ctdir) ) == 0:
        #continue
        if os.path.exists(baseF):
            covFiles.append(  baseF  )   
    else:
        continue
        #covFiles.append(sX[0] + "/2.mapping/"  +sX[1] +".base.coverage")   

d = {}
n = 1


R1 = {}

pat_ens = re.compile("ensembl_gn=(.*?);.*")
fE = open( '/mnt/disk01/project/Ngene.bed', 'r')
for line in fE:
    if not line.startswith('@') and not line.startswith("#") :
        x = line.rstrip().split('\t')
        geneInfo = x[3]  # x[-1]

        mat= re.match( pat_ens, geneInfo)  # x[-1]

        if mat != None:
            
            geneName = mat.group(1)
        else:
            geneName = geneInfo.split(';')[0]

    #geneName = x[-1].split(':')[-1]
    #n =  x[-1].split(':')[-2]
        if geneName not in d:
            d[geneName] = 1
            n = 1
        sn = "{0:03d}".format(n)
        R1[geneName + '_Region_' + sn] = '\t'.join(x[:3])
        n += 1

fE.close()

genomePosDict = {}
depthD = defaultdict(list)

for cfile in covFiles:
    c_h = open(cfile, 'r')
    for line in c_h:
        if not line.startswith('chrom'):
            LX = line.rstrip().split('\t')
 
            depthD[LX[0] + ':' +  LX[1]].append(int(LX[3]))

    c_h.close()


for i in R1:
    RX = R1[i]
    k = RX.split('\t')
    
    start = int(k[1])
    end   = int(k[2])  +  1
    for p in range(start, end):
        genomePosDict[k[0]+':' + str(p)] = i

Kp = {}
KN = 1
sampleDep = defaultdict(dict)
exon_pNum = defaultdict(list)
outMean = open("results_mean_PD.xls", 'w')
for i in sorted(depthD.keys()):

    D = depthD[i]
    m = np.mean(D)
    sd = np.std(D, ddof=1)
    if i not in genomePosDict:
        continue
    R_name = genomePosDict[i]
   
    if R_name not in Kp:
        Kp[R_name] = 1
        KN = 1

    for sn, j in enumerate(D):
        sampleDep[sn][R_name] = sampleDep[sn].get(R_name,0)+ j

    
    exon_pNum[R_name].append(i)
    KN += 1

outMean.write("RegionSN\tgenome_interv")
for i in covFiles:
    outMean.write("\t{0}".format(i.split('/')[-3]))
outMean.write('\n')

for i in sorted(exon_pNum):
    outMean.write('Region:{0}\t{1}'.format(i, R1[i].replace("\t",":")))
    for j in sampleDep:
        mean_d = (sampleDep[j][i] * 1.0)/ len(exon_pNum[i])
        outMean.write('\t{0}'.format(mean_d))
    outMean.write('\n')
outMean.close()
