#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO
from glob import glob
import numpy as np
import os
import argparse

def get_length_stat(f):
    the_Seq_Rec = SeqIO.parse( f , "fasta")
    sample_name = ".".join(f.split(".")[:-1])
    len_Of_Seq = []
    Tot_Len = 0

    for record in the_Seq_Rec:
        length = len(record.seq)
        len_Of_Seq.append(length)
        Tot_Len += length
    arr = np.array(len_Of_Seq)

    mx_len = np.max(arr)
    mean_len = np.mean(arr)
    General_stat = open("Gen_"+sample_name+"_stat.txt", "w")
    General_stat.write("The Total number of Tags in sample {0} is: {1}\n".format(sample_name, arr.size))
    General_stat.write("The Min length of Tags in sample {0} is: {1}\n".format(sample_name, np.min(arr)))
    General_stat.write("The Max length of Tags in sample {0} is: {1}\n".format(sample_name, mx_len))
    General_stat.write("The Mean length of Tags in sample {0} is: {1}\n".format(sample_name, mean_len))
    General_stat.close()

    hist, bin_edges = np.histogram(arr, bins = np.arange(0, mx_len + 20, 20), density = False)
#print (hist)
#print (bin_edges)

    Matri_Len = {}

    for n,i in enumerate(bin_edges):
        try:
            Matri_Len[str(i) +"-" + str(i+20)] = hist[n]
            #print ("{0}-{1}\t{2}".format(i,i+20, hist[n]))
        except IndexError:
            pass
            #print ("{0}-{1}\texceed Number range".format(i, i+20))
    stat_result = open(sample_name + "_stat.txt", 'w')
    stat_result.write("The distribution of length is:\n")
    Sort_Mat = sorted(Matri_Len.iteritems(), key = lambda x: int(x[0].split("-")[0]))
    for n,v in Sort_Mat:
        stat_result.write("{0}\t{1}\n".format(n, v))
    stat_result.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()     #simplifys the wording of using argparse as stated in the python tutorial
    parser.add_argument("-f", "--outputDir", type=str, action='store',  dest='FastaDir', help="input the bed file containing re digest length") # allows input of the file for reading
    parser.add_argument("-s", type=str, action='store',  dest='Summary_f',required = True, help="Provide a filename for the summary of tag numbers results")
    args = parser.parse_args()
    
    os.chdir(args.FastaDir)
    fas = glob("Unique_Sample_*.fa") 
    
    for i in fas:
        get_length_stat(i)
    
    StdPro = os.popen("awk /The Total number of Tags/ Gen_*_stat.txt")
    res = open(args.Summary_f, "w")
    
    res.write(StdPro.read())
    res.close()
