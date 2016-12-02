#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from collections import defaultdict
import argparse
import os
import vcf

# Abbreviation of Samtools_results is STR
STR = "SO_var.flt.vcf"

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
                if record.is_snp and record.INFO['AF1'] >= 0.5:
                    SNP_Num += 1
                    filtered_SNP.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format\
                     (record.CHROM, record.POS, record.ID, record.REF, record.ALT, record.QUAL, PMV,  1-record.INFO['AF1'], record.heterozygosity, Het_Freq))
            # check this position is a SNP and its AF <= 0.5, and then write this line into test_output_SNP.txt
                elif record.is_snp and record.INFO['AF1'] < 0.5:
                    SNP_Num += 1
                    filtered_SNP.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format\
                     (record.CHROM, record.POS, record.ID, record.REF, record.ALT, record.QUAL, PMV,  record.INFO['AF1'], record.heterozygosity, Het_Freq))

            # check INDEL and AF1 and write INDEL into test_output_IND.txt
                elif record.is_indel and record.INFO['AF1'] >= 0.5:
                    IND_Num += 1
                    filtered_INDEL.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format\
                     (record.CHROM, record.POS, record.ID, record.REF, record.ALT, record.QUAL, PMV,  1-record.INFO['AF1'], record.heterozygosity, Het_Freq))
            # check INDEL and AF1 and write INDEL into test_output_IND.txt
                elif record.is_indel and record.INFO['AF1'] < 0.5:
                    IND_Num += 1
                    filtered_INDEL.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format\
                     (record.CHROM, record.POS, record.ID, record.REF, record.ALT, record.QUAL, PMV, record.INFO['AF1'], record.heterozygosity, Het_Freq))

    filtered_SNP.close()
    filtered_INDEL.close()

    return (Var_Num, SNP_Num, IND_Num, PMV_count, Qual_Of_Samples)

if __name__ == "__main__":

    USAGE = "python %s.py -i vcf"%("filter_vcf")

    parser = argparse.ArgumentParser(description = USAGE)
    parser.add_argument("-i", "--input", action = "store", required = True, help = "the vcffile containing variants, including SNP and INDEL")
    parser.add_argument("-D", "--TD", type = int, default = 35, help = "the Depth in INFO field of Variants")
    parser.add_argument("-Q", "--Qst", type = float , default = 50,  help = "the QUAL of variant for all samples")
    parser.add_argument("-o", "--Pre", action = "store", default = "Fin", help = "the vcffile containing variants, including SNP and INDEL")
    
    args = parser.parse_args()

    # check if the vcffile exists
    if os.path.exists(args.input):
        Tot_VN, SNP_N, IND_N, PMV_Tot_C, QOS = filter_SNP(args.input, args.TD, args.Qst , args.Pre)
        print ("The number of all the variants filtered through func is {0}.".format(Tot_VN))
        print ("The number of all the SNPs filtered through func is {0}.".format(SNP_N))
        print ("The number of all the INDELs filtered through func is {0}.".format(IND_N))
    else:
        print ("The vcf file name is required, but it is not found in providing path.")

    for k, v in sorted(PMV_Tot_C.items(), key = lambda x: x[0]):
        print ("{0}\t{1}".format(k,v))


    QS_out = open("QUAL_Stats_samtools.txt", 'w')
    QS_out.write("Sample\tQual\n")
    for i in QOS:
        for j in QOS[i]:
            QS_out.write("{0}\t{1}\n".format(i, j))

    QS_out.close()



