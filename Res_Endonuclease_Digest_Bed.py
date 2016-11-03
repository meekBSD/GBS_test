#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import argparse
import os


USAGE = r"""
             %prog -e CTGCAG -i genomic.fna -o out_file
             example 1: C:\Users\user\AppData\Local\Programs\Python\Python35\python.exe %prog.py -i C:\Users\user\Downloads\GCF_sample_genomic.fna -e GT^AC -o loci_in_chr.txt
                 GTAC
             example 2: C:\Users\user\AppData\Local\Programs\Python\Python35\python.exe %prog.py -i sample_genome.fasta -o rest.txt
                 CTGCAG 
          """

parser = argparse.ArgumentParser(description = USAGE)

parser.add_argument('-i', "--input", action = "store", type= str, help = "the fasta file of Genome")
parser.add_argument("-e","--restriction", action = "store", type= str, dest = "RE", default = "CTGCA^G", help = "the restriction endonuclease you choose")
parser.add_argument('-o', "--output", action = "store", type = argparse.FileType("w"))

args = parser.parse_args()


FA = args.input

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

cut_site = args.RE.find("^")
enz_seq =args.RE.replace("^","")
print (enz_seq)

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
            args.output.write("{0}\t{1}\n".format(chr_n, enz_locus + 1 + cut_site))

args.output.close()


def get_lines(file_obj):
    items = [""]*2
    #items = ("", "")
    for number, line in enumerate(file_obj):
        
        slot = number % 2
        locus = line.rstrip()
        items[slot] = locus
        if number >=1 and items[0].split("\t")[0] == items[1].split("\t")[0]:
            yield (items)

bed_file = open(args.output.name + ".bed", 'w')
#with open("rest.txt", 'r') as fobj:
with open(args.output.name, 'r') as fobj:               
    for loci in get_lines(fobj):
        chr_name = loci[0].split("\t")[0]
        loci_s = int(loci[0].split("\t")[1])
        loci_e = int(loci[1].split("\t")[1])
        #length = abs(int(loci[0]) - int(loci[1]))
        if loci_s > loci_e:
            
            length = loci_s - loci_e
            bed_file.write("{0}\t{1}\t{2}\t{3}\n".format(chr_name, loci_e, loci_s, length))
        else:
            length = loci_e - loci_s
            bed_file.write("{0}\t{1}\t{2}\t{3}\n".format(chr_name, loci_s, loci_e, length))

bed_file.close()
        
