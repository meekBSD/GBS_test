# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from datetime import datetime
import argparse
import numpy as np

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
            #print ("{0}-{1}\t{2}".format(i,i+20, hist[n]))
        except IndexError:
            pass
            #print ("{0}-{1}\texceed Number range".format(i, i+20))
    return MatriMax_D

if __name__ == "__main__":
    
    print (datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print ("Data analysis start.")

    parser = argparse.ArgumentParser()     #simplifys the wording of using argparse as stated in the python tutorial
    parser.add_argument("-i", "--input", type=str, action='store',  dest='bedfile', help="input the bed file containing re digest length") # allows input of the file for reading
    parser.add_argument("-m", type=int, action='store',  dest='ArgM', default = 2000, help="Range border number in the right")
    parser.add_argument("-w", type=int, action='store',  dest='ArgW', default = 100, help="section number, the default value is 100")
    parser.add_argument("-o", "--output", help="Directs the length of fragments to a results file")
    args = parser.parse_args()

    Max_Length = args.ArgM
    window = args.ArgW
    
    # open the file contain endonuclease digest fragments start, end and length numbers
    f_all = open(args.bedfile , "r")
    
    # invoke the function calc_Frag_reEndo and return results to S_L and L_L
    S_L, L_L = calc_Frag_reEndo(f_all, Max_Length)
    
    # invoke the function Length_Dist and return the results as M_distribution
    M_distribution = Length_Dist(S_L, Max_Length, window)
    
    # Sort the dict M_distribution according the fragments length
    Sort_Mat = sorted(M_distribution.items(), key = lambda x: int(x[0].split("-")[0]))
    
    f_all.close()

    print ("The distribution of length is:")
    result = open(args.output, 'w')
    for n,v in Sort_Mat:
        print ("{0}\t{1}".format(n, v))
        result.write("{0}\t{1}\n".format(n, v))
    result.close()
    print ("The number of sequences longer than %d is: %d"%(Max_Length, len(L_L)))
    print (datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print ("Analysis Done.")
