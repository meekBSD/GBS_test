#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from datetime import datetime
import numpy as np
import argparse
import os
import sys


USAGE = r"""
             %prog -e CTGCAG -i genomic.fna -o out_file
             example 1: C:\Users\user\AppData\Local\Programs\Python\Python35\python.exe %prog.py -i C:\Users\user\Downloads\GCF_sample_genomic.fna -e GT^AC -l 600 -o loci_digested
                 GTAC 
                 
             example 2: C:\Users\user\AppData\Local\Programs\Python\Python35\python.exe %prog.py -i sample_genome.fasta -l 800 -w 20 -s 3000
                 CTGCAG 
          """

def read_fasta(fp):
    name, seq = (None, [])
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield (name, ''.join(seq))
            name, seq = (line, [])
        else:
            seq.append(line.upper())
    if name:
        yield (name, ''.join(seq))

## Define a function get the length of Restriction Endonuclease digest fragments and store them in two list
def calc_Frag_reEndo(file_in, MaxLen):
    Shorter_Max = []
    Longer_thanMax = []
    #file_in = open(bed_file, 'r')
    for each_Line in file_in:
        L_field = each_Line.rstrip().split("\t")[-1]
        fragment_L = int(L_field)
        if fragment_L <= MaxLen:
            Shorter_Max.append(fragment_L)
        else:
            Longer_thanMax.append(fragment_L)
            
    return (Shorter_Max, Longer_thanMax)

## Create length distribution raw data and store them in dict for histogram graph                                                                             
def Length_Dist(a, ax, win):
    arr = np.array(a)
    hist, bin_eg = np.histogram(arr, bins = np.arange(0,ax+1,win), density =False)
    MatriMax_D = {}
    for n,i in enumerate(bin_eg):
        try:
            MatriMax_D[str(i) +"-" + str(i + win)] = hist[n]
        except IndexError:
            pass
    return MatriMax_D

def get_lines(file_obj):
    items = [""]*2
    for number, line in enumerate(file_obj):
        slot = number % 2
        locus = line.rstrip()
        items[slot] = locus
        if number >=1 and items[0].split("\t")[0] == items[1].split("\t")[0]:
            yield (items)

if __name__ == "__main__":
    
    print (datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print ("Data analysis start.")
    parser = argparse.ArgumentParser(description = USAGE)   #simplifys declare usage of this script

    parser.add_argument('-i', "--input", action = "store", type= str, help = "the fasta file of Genome")
    parser.add_argument("-e","--restriction", action = "store", type= str, dest = "RE", default = "CTGCA^G", help = "the restriction endonuclease you choose")
    parser.add_argument("-l","--Frag_L", type= int, default = 600, help = "the maximal fragment length retained in the bed file, default is 600, if set as 0, output all regions")
    parser.add_argument("-m","--Frag_m", type= int, default = 50, help = "the minimal fragment length retained in the bed file, default is 50")
    
    parser.add_argument("-s", type=int, action='store',  dest='ArgM', default = 2000, help="Maximal size of digested fragments")
    parser.add_argument("-w", type=int, action='store',  dest='ArgW', default = 100, help="gradient number, the default value is 100")
    parser.add_argument('-o', "--output", action = "store", default = "DigestedBed", type = argparse.FileType("w"), help="prefix of bed file name containing re digest length")
    parser.add_argument("-d", "--distribution", default = "frag_length.txt", help="Directs the length of fragments to a results file")

    args = parser.parse_args()

    if args.Frag_L < args.Frag_m:
        sys.exit("the maximal length is small than minimal length, please check the number of parameter '-l'. ")

    FA = args.input
    fa_dir = os.path.dirname(FA)
    print(fa_dir)

    cut_site = args.RE.find("^")
    enz_seq =args.RE.replace("^","")
    print (enz_seq)

    Max_Length = args.ArgM
    window = args.ArgW

    with open(FA, "r") as f_handle:
        for contig_n, seq in read_fasta(f_handle):
            chr_n = contig_n.strip(">").split()[0]
            #print (args.RE)
            args.output.write("{0}\t{1}\n".format(chr_n, 0))
            enz_locus = -1
            while True:
             # find enzyme sequence location from the start
                enz_locus = seq.find(enz_seq, enz_locus + 1)
                # break if not found
                if enz_locus == -1:
                    args.output.write("{0}\t{1}\n".format(chr_n, len(seq)))
                    break
                args.output.write("{0}\t{1}\n".format(chr_n, enz_locus + cut_site))

    args.output.close()

    txt_f = open(os.path.join(fa_dir, args.output.name + ".txt"), 'w')
    bed_file = open(os.path.join(fa_dir, args.output.name + ".bed"), 'w')

    with open(args.output.name, 'r') as fobj:               
        for loci in get_lines(fobj):
            chr_name = loci[0].split("\t")[0]
            loci_s = int(loci[0].split("\t")[1])
            loci_e = int(loci[1].split("\t")[1])
            length = abs(loci_e - loci_s)
            if loci_s > loci_e:
                txt_f.write("{0}\t{1}\t{2}\t{3}\n".format(chr_name, loci_e, loci_s, length))
                if args.Frag_L == 0:
                    bed_file.write("{0}\t{1}\t{2}\t.\t.\t+\n".format(chr_name, loci_e, loci_s, length))
                elif length <= args.Frag_L and length >= args.Frag_m:
                    bed_file.write("{0}\t{1}\t{2}\t.\t.\t+\n".format(chr_name, loci_e, loci_s, length))
            else:
                txt_f.write("{0}\t{1}\t{2}\t{3}\n".format(chr_name, loci_s, loci_e, length))
                if args.Frag_L == 0:
                    bed_file.write("{0}\t{1}\t{2}\t.\t.\t+\n".format(chr_name, loci_s, loci_e, length))
                elif length <= args.Frag_L and length >= args.Frag_m:
                    bed_file.write("{0}\t{1}\t{2}\t.\t.\t+\n".format(chr_name, loci_s, loci_e, length))
    txt_f.close()
    bed_file.close()

  
    # open the file contain endonuclease digest fragments start, end and length numbers
    f_all = open(os.path.join(fa_dir, args.output.name + ".txt") , "r")
    
    # invoke the function calc_Frag_reEndo and return results to S_L and L_L
    S_L, L_L = calc_Frag_reEndo(f_all, Max_Length)
    
    # invoke the function Length_Dist and return the results as M_distribution
    M_distribution = Length_Dist(S_L, Max_Length, window)
    
    # Sort the dict M_distribution according the fragments length
    Sort_Mat = sorted(M_distribution.items(), key = lambda x: int(x[0].split("-")[0]))
    
    f_all.close()

    print ("The distribution of length is:")
    result = open(os.path.join(fa_dir,args.distribution), 'w')
    for n,v in Sort_Mat:
        print ("{0}\t{1}".format(n, v))
        result.write("{0}\t{1}\n".format(n, v))
    result.close()
    print ("The number of sequences longer than %d is: %d"%(Max_Length, len(L_L)))
    print (datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print ("Analysis Done.")

        
