#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import subprocess
import os

## check FFPE samples, get previous targets files
f1 = open("samples.txt", "r")

previousD = "/mnt/disk01/project"
targets_List = []
f_t = open(  "targets_samples.txt", "w"   )


for line in f1:
    k = line.rstrip().split("/")
    baseFileN = k[-1].split(".")[0]
    if baseFileN == k[-3] and 'FL' in k[-3]:
        previousTarFile = previousD + baseFileN + "S01/VD.mapping/" + baseFileN + "S01.target.coverage"
        #
        if os.path.exists(previousTarFile):
            #print("{0}\t{1}\t{2}".format( k[-1], k[-3] ,  previousTarFile  )   )
            targets_List.append(line.rstrip())
            f_t.write("{0}\t{1}\n".format(  previousD + baseFileN + "S01", baseFileN + "S01"  ))

f_t.close()
f1.close()

## merge all target mean coverage
if not os.path.exists("results_mean_PD.xls"):
    subprocess.call("python calc_sampleMean.py", shell=True)

## find low coverage and check coverage values in new probe sets
targetsPos_D = {}
fh = open(  "results_mean_PD.xls",  "r"  )

fw = open(  "low_cov_previous.xls", "w")
for line in fh:
    if not line.startswith("RegionSN"):
        x = line.rstrip().split("\t")
        dv = [float(i) for i in x[2:]]
        gx = x[1].split(":")

        if all(n < 300 for n in dv):
            for j in range(int(gx[1]), int(gx[2])+1):
                targetsPos_D[gx[0]+":" + str(j)]= x[1]
            fw.write(line)    

fw.close()
fh.close()

outD = {}
c_tar = {}
for i in targetsPos_D:
    c_tar[targetsPos_D[i]] = c_tar.get(targetsPos_D[i],0) + 1

for f in targets_List:
    outD = {}
    ftf = open(f, "r")
    for line in ftf:
        if not line.startswith("chrom"):
            x = line.rstrip().split("\t")
            s = x[0] + ":" + x[1]
            if s in targetsPos_D:
                outD[targetsPos_D[s]] = outD.get(targetsPos_D[s], 0) + int(x[3]) 

    ftf.close()


    outtf = open(os.path.basename(f) + ".newcov.txt", "w")
    for i in outD:
        outtf.write("{0}\t{1}\t{2}\t{3}\n".format(i, c_tar[i], outD[i], outD[i]*1.0/c_tar[i]) ) 
    outtf.close()

## merge results




