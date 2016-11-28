#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import vcf
import gzip

# Abbreviation of Samtools_results is STR
STR = "SO_var.flt.vcf"

vcf_reader = vcf.Reader(open(STR, 'r'))
vcf_W = vcf.Writer(open("test_sample_filtered_var.vcf", 'w'), vcf_reader)
filtered_f = open("test_sample_SNP.txt", 'w')

Samtools_SNP_filterd = set()
n = 0

# write header to the result file.
filtered_f.write("Chr_Pos")
for s in vcf_reader.samples:
    filtered_f.write("\t{0}".format(s))
filtered_f.write("\n")

for num, record in enumerate(vcf_reader):
    if record.INFO['DP'] >= 200:
        n += 1
        Samtools_SNP_filterd.add("{0}#{1}".format(record.CHROM, record.POS))
        vcf_W.write_record(record)

        #if n <= 1:
        #    filtered_f.write("Chr_Pos")
        #    for sh in record.samples:
        #        filtered_f.write("\t{0}".format(sh))
        #    filtered_f.write("\n")

        filtered_f.write("{0}-{1}".format(record.CHROM, record.POS))
        for s in record.samples:
            filtered_f.write("\t{0}".format(s['DP']))
        filtered_f.write("\n")

print ("The number of all the variants filtered is {0}.".format(n))
print ("Filtration of VCF file finished.")

filtered_f.close()

print (len(Samtools_SNP_filterd))
