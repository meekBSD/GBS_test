#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from collections import defaultdict
import vcf
import gzip

# the script mainly use Quality and Depth as filteration standard and calculate MAF for vcf
# In addition, the freebayes output and samtools output is different so , the Allele Frequency field is different , please be careful 
# when using AF1 or AF

# Abbreviation of Samtools_results is STR
'''
STR = "SO_var.flt.vcf"

vcf_reader = vcf.Reader(open(STR, 'r'))
vcf_W = vcf.Writer(open("test_out_half_filtered_var.vcf", 'w'), vcf_reader)
filtered_f = open("test_output_half_SNP.txt", 'w')
filtered_INDEL = open("test_output_half_IND.txt" ,'w')
filtered_f.write("Chr\tPOS\tRec_ID\tREF\tALT\tRec_Qua\tMAF\tHET\n")
filtered_INDEL.write("Chr\tPOS\tRec_ID\tREF\tALT\tRec_Qua\tMAF\tHET\n")

Samtools_SNP_filterd = set()

Dp_more_than5_n = 0
Fil_SNP_Num = 0
Fil_IND_Num = 0
for record in vcf_reader:
    if len([s for s in record.samples if s['DP'] > 4]) > 22: 
    # if all(s['DP'] >=2 for s in record.samples): # use this standard for samtools calling results only 18 SNP sites 
    # passed this standard 3 INDEL passed. If I set this value 3, 4, or 5, no variant passed.
        Dp_more_than5_n += 1
        if record.INFO['DP'] >= 220 and record.QUAL > 50:
            #Samtools_SNP_filterd.add("{0}#{1}".format(record.CHROM, record.POS))
            vcf_W.write_record(record)
            # check this position is a SNP and its AF > 0.5, and then write this line into test_output_SNP.txt
            Allele_Freq = record.INFO['AF1']
            if record.is_snp and Allele_Freq >= 0.5:
                Fil_SNP_Num += 1
                Samtools_SNP_filterd.add("{0}#{1}".format(record.CHROM, record.POS))
                filtered_f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format\
                     (record.CHROM, record.POS, record.ID, record.REF, record.ALT,record.QUAL, 1 - Allele_Freq, record.heterozygosity))
            # check this position is a SNP and its AF <= 0.5, and then write this line into test_output_SNP.txt
            elif record.is_snp and Allele_Freq < 0.5:
                Fil_SNP_Num += 1
                Samtools_SNP_filterd.add("{0}#{1}".format(record.CHROM, record.POS))
                filtered_f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format\
                     (record.CHROM, record.POS, record.ID, record.REF, record.ALT,record.QUAL, Allele_Freq, record.heterozygosity))

            # check INDEL and AF1 and write INDEL into test_output_IND.txt
            elif record.is_indel and Allele_Freq >= 0.5:
                Fil_IND_Num += 1
                filtered_INDEL.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format\
                     (record.CHROM, record.POS, record.ID, record.REF, record.ALT,record.QUAL, 1 - Allele_Freq , record.heterozygosity))
            # check INDEL and AF1 and write INDEL into test_output_IND.txt
            elif record.is_indel and Allele_Freq < 0.5:
                Fil_IND_Num += 1
                filtered_INDEL.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format\
                     (record.CHROM, record.POS, record.ID, record.REF, record.ALT,record.QUAL, Allele_Freq, record.heterozygosity))

print ("The number of all the variants filtered is {0}.".format(Dp_more_than5_n))
print ("The number of all the SNPs filtered is {0}.".format(Fil_SNP_Num))
print ("The number of all the INDELs filtered is {0}.".format(Fil_IND_Num))
print ("Filtration of VCF file finished.")

filtered_f.close()
filtered_INDEL.close()
'''

# abbreviation of FreeBayes_results = FBR
FBR = "ZiJu.vcf.gz"

FreeBayes_V_Reader = vcf.Reader(filename= FBR)


#vcf_reader = vcf.Reader(open("five_var.flt.vcf", 'r'))
Fre_Bay_vcf_W = vcf.Writer(open("test_FreeBayes_ninty_filtered_var.vcf", 'w'), FreeBayes_V_Reader)
Fre_Bay_filtered_f = open("test_FreeBayes_ninty_SNP.txt", 'w')
Fre_Bay_filtered_INDEL = open("test_FreeBayes_ninty_IND.txt" ,'w')
Fre_Bay_filtered_f.write("Chr\tPOS\tRec_ID\tREF\tALT\tRec_Qua\tMAF\tHET\n")
Fre_Bay_filtered_INDEL.write("Chr\tPOS\tRec_ID\tREF\tALT\tRec_Qua\tMAF\tHET\n")

Qual_Of_Samples = defaultdict(list)
Free_Bayes_filtered = set()

Free_Dp_more_than5_n = 0
Free_Fil_SNP_Num = 0
Free_Fil_IND_Num = 0
for record in FreeBayes_V_Reader:
    if len([s for s in record.samples if s['DP'] > 4]) > 39:  #
    # if all(s['DP'] >= 2 for s in record.samples):  # this standand is too rigorously for SNP filtration
        
        # assign Sample quality to the Dict Qual_of_Samples

        for sq in record.samples:
            if sq['DP'] >= 10:
                Qual_Of_Samples[sq.sample].append(record.QUAL)

        Free_Dp_more_than5_n += 1
        if record.INFO['DP'] >= 220 and record.QUAL >50:
            Fre_Bay_vcf_W.write_record(record)
            #Free_Bayes_filtered.add("{0}#{1}".format(record.CHROM, record.POS))
            # check this position is a SNP and its AF > 0.5, and then write this line into test_output_SNP.txt
            Allele_Freq_thisVCF = record.INFO['AF'][0]
            if record.is_snp and Allele_Freq_thisVCF >= 0.5:
                Free_Fil_SNP_Num += 1
                Free_Bayes_filtered.add("{0}#{1}".format(record.CHROM, record.POS))
                Fre_Bay_filtered_f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format\
                     (record.CHROM, record.POS, record.ID, record.REF, record.ALT,record.QUAL, 1 - Allele_Freq_thisVCF, record.heterozygosity))
            # check this position is a SNP and its AF <= 0.5, and then write this line into test_output_SNP.txt
            elif record.is_snp and Allele_Freq_thisVCF < 0.5:
                Free_Fil_SNP_Num += 1
                Free_Bayes_filtered.add("{0}#{1}".format(record.CHROM, record.POS))
                Fre_Bay_filtered_f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format\
                     (record.CHROM, record.POS, record.ID, record.REF, record.ALT,record.QUAL, Allele_Freq_thisVCF, record.heterozygosity))

            # check INDEL and AF1 and write INDEL into test_output_IND.txt
            elif record.is_indel and Allele_Freq_thisVCF >= 0.5:
                Free_Fil_IND_Num += 1
                Fre_Bay_filtered_INDEL.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format\
                     (record.CHROM, record.POS, record.ID, record.REF, record.ALT, record.QUAL, 1 - Allele_Freq_thisVCF , record.heterozygosity))
            # check INDEL and AF1 and write INDEL into test_output_IND.txt
            elif record.is_indel and Allele_Freq_thisVCF < 0.5:
                Free_Fil_IND_Num += 1
                Fre_Bay_filtered_INDEL.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format\
                     (record.CHROM, record.POS, record.ID, record.REF, record.ALT, record.QUAL, Allele_Freq_thisVCF , record.heterozygosity))

print ("The number of all the variants filtered through FreeBayes is {0}.".format(Free_Dp_more_than5_n))
print ("The number of all the SNPs filtered through FreeBayes is {0}.".format(Free_Fil_SNP_Num))
print ("The number of all the INDELs filtered through FreeBayes is {0}.".format(Free_Fil_IND_Num))
print ("Filtration of VCF file through FreeBayes finished.")

Fre_Bay_filtered_f.close()
Fre_Bay_filtered_INDEL.close()

QS_out = open("QUAL_Stats.txt", 'w')
QS_out.write("Sample\tQual\n")
for i in Qual_Of_Samples:
    for j in Qual_Of_Samples[i]:
        QS_out.write("{0}\t{1}\n".format(i, j))

QS_out.close()

'''
print ("SNP number in samtools analysis results is {0}.".format(len(Samtools_SNP_filterd)))
print ("SNP number in FreeBayes analysis results is {0}.".format(len(Free_Bayes_filtered)))

S_NFB_number = len(Samtools_SNP_filterd - Free_Bayes_filtered) 
F_NST_number = len(Free_Bayes_filtered - Samtools_SNP_filterd)

print ("SNP number in samtools results but not in FreeBayes results is {0}".format(S_NFB_number))
print ("SNP number in FreeBayes results but not in samtools results is {0}".format(F_NST_number))
print ("The number of Common SNP in Samtools and FreeBayes is {0}".format(len(Samtools_SNP_filterd & Free_Bayes_filtered)))
'''
