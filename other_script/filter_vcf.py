#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from collections import defaultdict
import argparse
import os
import sys
import vcf
import copy

def coverage_stats(f, max, interval):     # also can set parameter max = 500 and interval = 50, but it can't be used in argparse conveniently.
    depth = {}
    for i in range(50, max, interval):
        depth[i] = 0
    depth[10] = 0

    all_sample_dep = {}
    
    AD = list(depth.keys())
    AD.sort()
    a = vcf.Reader(filename = f)
    for samp in a.samples:
        all_sample_dep[samp] = copy.deepcopy(depth)
    
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

# Abbreviation of total depth is TotDP, Abbreviation of filter standard of Quality is filter_st

def filter_SNP(vcf_file, TotDP, filter_st, prefix):

    vcf_reader = vcf.Reader(filename = vcf_file)
    vcf_W = vcf.Writer(open(prefix + "_out_filtered_var.vcf", 'w'), vcf_reader)
    filtered_SNP = open(prefix + "_output_SNP.txt", 'w')
    filtered_INDEL = open(prefix + "_output_IND.txt" ,'w')
    filtered_SNP.write("Chr\tPOS\tRec_ID\tREF\tALT\tQUAL\tPMV\tMAF\tHET\tHet_Obs\n")
    filtered_INDEL.write("Chr\tPOS\tRec_ID\tREF\tALT\tQUAL\tPMV\tMAF\tHET\tHet_Obs\n")

    Var_Num = 0
    SNP_Num = 0
    IND_Num = 0

    Qual_Of_Samples = defaultdict(list)
    C = float(len(vcf_reader.samples))

    PMV_count = {
        80: 0 ,
        70: 0 ,
        60: 0 ,
        40: 0 ,
        20: 0 ,
        10: 0 ,
         5: 0 ,
         1: 0
    }

    for record in vcf_reader:

        if 'AF1' in record.INFO:
            af = record.INFO['AF1']       # af is allele frequency
        elif 'AF' in record.INFO:
            af = record.INFO['AF'][0]
        #PM_Num = round(len(record.samples) * filter_st)
        #if len([s for s in record.samples if s['DP'] > 4]) >= PM_Num:

        #if all(s['DP'] >=4 for s in record.samples):
        if record.QUAL > filter_st:
            PMV_r = len([s for s in record.samples if s['DP'] >2]) * 100 / C
            PMV = 100 - PMV_r            
            for k in PMV_count:
                if round(PMV_r) > 100 - k:
                    PMV_count[k] += 1

            GT_a = [i['GT'] for i in record.samples]
            # print GT_a
            Het_Freq = GT_a.count("0/1")/C

            for sq in record.samples:
           
                if sq['DP'] >= 10:
                    Qual_Of_Samples[sq.sample].append(record.QUAL)

            Var_Num += 1
            if record.INFO['DP'] >= TotDP:
                vcf_W.write_record(record)
            # check this position is a SNP and its AF > 0.5, and then write this line into test_output_SNP.txt
                if record.is_snp and af >= 0.5:
                    SNP_Num += 1
                    filtered_SNP.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format\
                     (record.CHROM, record.POS, record.ID, record.REF, record.ALT, record.QUAL, PMV,  1-af, record.heterozygosity, Het_Freq))
            # check this position is a SNP and its AF <= 0.5, and then write this line into test_output_SNP.txt
                elif record.is_snp and af < 0.5:
                    SNP_Num += 1
                    filtered_SNP.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format\
                     (record.CHROM, record.POS, record.ID, record.REF, record.ALT, record.QUAL, PMV,  af, record.heterozygosity, Het_Freq))

            # check INDEL and AF1 and write INDEL into test_output_IND.txt
                elif record.is_indel and af >= 0.5:
                    IND_Num += 1
                    filtered_INDEL.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format\
                     (record.CHROM, record.POS, record.ID, record.REF, record.ALT, record.QUAL, PMV,  1-af, record.heterozygosity, Het_Freq))
            # check INDEL and AF1 and write INDEL into test_output_IND.txt
                elif record.is_indel and af < 0.5:
                    IND_Num += 1
                    filtered_INDEL.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format\
                     (record.CHROM, record.POS, record.ID, record.REF, record.ALT, record.QUAL, PMV, af, record.heterozygosity, Het_Freq))

    filtered_SNP.close()
    filtered_INDEL.close()

    return (Var_Num, SNP_Num, IND_Num, PMV_count, Qual_Of_Samples)

if __name__ == "__main__":

    USAGE = "python %s.py -i vcf"%("filter_vcf")

    parser = argparse.ArgumentParser(description = USAGE)
    parser.add_argument("-i", "--input", action = "store", required = True, help = "the vcf file containing variants, including SNP and INDEL")
    parser.add_argument("-D", "--TD", type = int, default = 35, help = "the Depth in INFO field of Variants")
    parser.add_argument("-Q", "--Qst", type = float , default = 50,  help = "the QUAL standard of variant filteration for all samples")
    parser.add_argument("-o", "--Pre", action = "store", default = "Fin", help = "the prefix of output file")
    parser.add_argument("-sf", "--Samp_Qua_Suffix", action = "store", default = "QS", help = "the Quality statistics file containing each sample and qualities of its variants")
    
    parser.add_argument("-m", "--MaxDepth", type = int, default = 500, help = "the maximal Depth in INFO field of VCFfile")
    parser.add_argument("-t", "--interv", type = int , default = 50, help = "the depth interval of variant in each_sample")
    parser.add_argument("-c", "--Cov", action = "store", default = "SNP_Cov", help = "the filename including SNP depth of each sample")
    args = parser.parse_args()

    # check if the vcffile exists
    if not os.path.exists(args.input):
        sys.exit("The vcf file name is required, but it is not found.\nPlease check it and provide correct file path.")
                
    Tot_VN, SNP_N, IND_N, PMV_Tot_C, QOS = filter_SNP(args.input, args.TD, args.Qst , args.Pre)
    print ("The number of all the variants filtered through func is {0}.".format(Tot_VN))
    print ("The number of all the SNPs filtered through func is {0}.".format(SNP_N))
    print ("The number of all the INDELs filtered through func is {0}.".format(IND_N))

    for k, v in sorted(PMV_Tot_C.items(), key = lambda x: x[0]):
        print ("{0}\t{1}".format(k,v))

    QS_out = open("QUAL_Stats_"+args.Samp_Qua_Suffix +".txt", 'w')
    QS_out.write("Sample\tQual\n")
    for i in QOS:
        for j in QOS[i]:
            QS_out.write("{0}\t{1}\n".format(i, j))
    QS_out.close()


    Dep_Gradient, Samples_Depth = coverage_stats(args.input, args.MaxDepth, args.interv)
    print ("The samples in this vcf are: ")
    for each_sample in Samples_Depth:
        print ("The depth of {0}:".format(each_sample))
        Coverage_Dict = Samples_Depth[each_sample]
        for j in sorted(Coverage_Dict.items(), key = lambda x: x[0]):
            print ("{0}\t{1}".format(j[0], j[1]))
    
    output_file = open(args.Cov, 'w')

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


