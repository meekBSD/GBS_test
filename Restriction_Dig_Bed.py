#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import argparse
import os


USAGE = r"""
             %prog -e CTGCAG -o out_file \
             
             provide a Restriction endonuclease with cut site
             example 1: C:\Users\user\Desktop>C:\Users\user\AppData\Local\Programs\Python\Python35\python.exe 1.py -e TTCGA^A -o rest.txt \
                 TTCOT \
             example 2: C:\Users\user\Desktop>C:\Users\user\AppData\Local\Programs\Python\Python35\python.exe 1.py -o rest.txt \
                 CTGCAG \
                 CTGCA^G
          """

parser = argparse.ArgumentParser(description = USAGE)

# parser.add_argument('-e', "--restriction", action = "store", type= str, default = "CTGCA^G", help = "the restriction enzyme you choose")
parser.add_argument("--restriction", action = "store", type= str, dest = "RE", default = "CTGCA^G", help = "the restriction enzyme you choose")
parser.add_argument('-o', "--output", action = "store", type = argparse.FileType("w"))

args = parser.parse_args()


os.chdir("C:\\Users\\user\\Downloads\\")

def read_fasta(fp):
    name, seq = (None, [])
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield (name, ''.join(seq))
            name, seq = (line, [])
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))

cut_site = args.RE.find("^")
enz_seq =args.RE.replace("^","")
print (enz_seq)

with open("test_restriction.fa") as f_handle:
    for contig_n, seq in read_fasta(f_handle):
        #print (args.restriction)
        print (args.RE)
        #args.output.write (seq+"##\n")
        enz_locus = -1
        while True:
            # find enzyme sequence location from the start
            enz_locus = seq.find(enz_seq, enz_locus + 1)
            # break if not found
            if enz_locus == -1:
                break
            args.output.write("{0}\n".format(enz_locus + cut_site))

args.output.close()


os.chdir("C:\\Users\\user\\Desktop\\")

def get_lines(file_obj):
    items = [""]*2
    #items = ("", "")
    for number, line in enumerate(file_obj):
        slot=number % 2
        items[slot] = line.rstrip()
        if number >=1:
            yield (items)

bed_file = open(args.output.name + ".bed", 'w')
#with open("rest.txt", 'r') as fobj:
with open(args.output.name, 'r') as fobj:               
    for loci in get_lines(fobj):
        #print (loci)
        #length = abs(int(loci[0]) - int(loci[1]))
        if int(loci[0]) > int(loci[1]):
            
            length = int(loci[0]) - int(loci[1])
            bed_file.write("{0}\t{1}\t{2}\t{3}\n".format("test_contig", loci[1], loci[0], length))
        else:
            length = int(loci[1]) - int(loci[0])
            bed_file.write("{0}\t{1}\t{2}\t{3}\n".format("test_contig", loci[0], loci[1], length))

bed_file.close()
