#!/usr/bin/env python
# -*- coding: UTF-8 -*-


from glob import glob
from collections import defaultdict
import os
import numpy as np
import pandas as pd


IN_DIR  = "coverage_files/"
OUT_DIR = "Zscore_results/"

depthFiles = glob( IN_DIR + "*base.coverage")
samples = [ "1821314", "1821316" , "1821322"  , "1821326"  , "1821328" , "1821330" , "1821334" , "1821336" , "1821338" , "1821340" , "1821342" , "1821346" , "1821348" , "1821350" ,  "1821352" ,  "1821354" ,  "1821358"  ,  "1821360" , "1821362"  ,  "1821374"  , "1821376" ,  "1821378" , "1821382"  ,  "1821386" ,  "1821388"  ,  "1821390"  ,  "1821396"  ,  "1821400"  ,  "1860866" ]

frames = []

for i in depthFiles:
    sample = os.path.basename(i).split('.')[0]
    
    if sample in samples:
    
        dataF = pd.read_csv(i , sep = '\t' ,  header=0)
        data_n = dataF.iloc[:,2:4]
        data_n.columns = [ "target"  ,  "cov_" + sample  ]
    
        frames.append(data_n.reset_index(drop=True))
        #print(data_n.head(4))
    
all_D = pd.concat(frames, sort=False, axis=1)
df = all_D.loc[:,~all_D.columns.duplicated()]

df.loc[:,"mean_COV"] = df.loc[:,"cov_1821314": "cov_1860866"].mean(axis=1, numeric_only=True)
df.loc[:,"std_COV"]  = df.loc[:,"cov_1821314": "cov_1860866"].std(axis=1, numeric_only=True)

df.to_csv( OUT_DIR  +  "raw_cov.xls", sep='\t', index = False)

os.chdir(OUT_DIR)

outZfile = open("test_Z-Score.xls", 'w')

handle_A = open("raw_cov.xls", "r")
for line in handle_A:
    if line.startswith("target"):
        outZfile.write(line.rstrip().replace("cov", "zscore") + "\tavg_zscore\n")
    else:
        
        x = line.rstrip().split('\t')
        outZfile.write(x[0] + '\t')
        
        std = float(x[31])
        z_list = []
        if std > 0:
            for j in x[1:30]:
                z_score = (float(j) - float(x[30]))/ std
                outZfile.write( "{0:.4f}\t".format(z_score))
                z_list.append(z_score)
        else:
            for j in x[1:30]:
                z_score = "NA"
                outZfile.write( "{0}\t".format(z_score))
                z_list.append(0)
        
        Arr_z = np.array(z_list)
        avg_zs = np.mean(Arr_z)
        outZfile.write( "{0}\t{1}\t{2:.4f}\n".format(x[30], x[31], avg_zs))
        
outZfile.close()
handle_A.close()

