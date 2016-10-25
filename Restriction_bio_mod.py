#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from Bio.Restriction import Restriction
from Bio import SeqIO

#print (Restriction.Sau3AI.site)
print (Restriction.PstI.site)

records = SeqIO.parse("C:\\Users\\user\\Downloads\\test_restriction.fa", "fasta")

for i in records:
    name = i.id
    #digest = Restriction.ApeKI.catalyse(i.seq)
    hits = Restriction.PstI.search(i.seq)

    for d in hits:
        print ("{0}\t{1}\t{2}".format(name, str(d), str(d+1)))


## Ref
##        https://github.com/daler/rdbio-scripts/blob/master/sequenceFiles/restriction-finder.py
        
##        http://coco.sam.pitt.edu/~emeneses/python/lecture8.pdf   Page 20

        ##
        ## from Bio.Restriction import Restriction
        ## from Bio import SeqIO

        ## print (Restriction.Sau3AI.site)

        ## records = SeqIO.parse(r"C:\Users\user\Downloads\hs_alt_CHM1_1.1_chr2.fa", "fasta")

        ## for i in records:
        ##     print (i.id)
        ##     digest = Restriction.ApeKI.catalyse(i.seq)
        ##     print (len(digest))

        ## for d in digest:
        ##     print (d)
