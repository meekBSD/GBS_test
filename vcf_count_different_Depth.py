#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from collections import defaultdict
import os
import argparse
import vcf
import copy

def coverage_stats(f, max, interval):     # also can set parameter max = 500 and interval = 50, but it can't be used in argparse conveniently.
    depth = {}
    for i in range(50, max, interval):
        depth[i] = 0

    depth[10] = 0
    
    all_sample_dep = {}
    
    a = vcf.Reader(filename = f)
    for samp in a.samples:
        all_sample_dep[samp] = copy.deepcopy(depth)
    
    AD = list(depth.keys())
    AD.sort()
    for record in a:
        for s in record.samples:
            s_n = s.sample
            Dep_curr = s['DP']
            for n,num_r in enumerate(AD[:-1]):                  # num_r is the bottom number range margin                
                if Dep_curr >= num_r and Dep_curr < AD[n+1]:    # AD[n+1] is top number margin
                    all_sample_dep[s_n][num_r] += 1
                else:
                    continue

    return (depth, all_sample_dep)

if __name__ == "__main__":

    USAGE = "python %s.py -i vcf"%("filter_vcf")

    parser = argparse.ArgumentParser(description = USAGE)

    parser.add_argument("-i", "--input", action = "store", required = True, help = "the vcffile containing variants, including SNP and INDEL")
    parser.add_argument("-m", "--MaxDepth", type = int, default = 500, help = "the Depth in INFO field of VCFfile")
    parser.add_argument("-t", "--interv", type = int , default = 50, help = "the depth of variant in each_sample")
    parser.add_argument("-o", "--Pre", action = "store", default = "SNP_Cov", help = "the vcffile containing variants, including SNP and INDEL")

    args = parser.parse_args()

    # check if the vcffile exists
    if os.path.exists(args.input):
        Dep_Gradient, Samples_Depth = coverage_stats(args.input, args.MaxDepth, args.interv)
        print ("The samples in this vcf are: ")
        for each_sample in Samples_Depth:
            print ("The depth of {0}:".format(each_sample))
            Coverage_Dict = Samples_Depth[each_sample]
            for j in sorted(Coverage_Dict.items(), key = lambda x: x[0]):
                print ("{0}\t{1}".format(j[0], j[1]))
        output_file = open(args.Pre, 'w')

        sample_list = []
        output_file.write("Minimal_Cov")
        for i in Samples_Depth:
            output_file.write("\t{0}".format(i))
            sample_list.append(i)
        output_file.write("\n")

        Depths = [x[0] for x in sorted(Dep_Gradient.items(), key = lambda x: x[0])]
        for Gra in Depths:
            output_file.write("{0}".format(Gra))
            for s in sample_list:
                output_file.write("\t{0}".format(Samples_Depth[s][Gra]))
            output_file.write("\n")            
        output_file.close()
        print ("The sample number in vcf is {0}.\n Analysis finished.".format(len(Samples_Depth)))
    else:
        print ("The vcf file name is required, but it is not found in providing path.")
