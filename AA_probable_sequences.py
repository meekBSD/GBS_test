#!/usr/bin/env python
# -*- coding: UTF-8 -*-


from itertools import product

codons= {
    'TCA' : 'S',    # Serine
    'TCC' : 'S',    # Serine
    'TCG' : 'S',    # Serine
    'TCT' : 'S',    # Serine
    'TTC' : 'F',    # Phenylalanine
    'TTT' : 'F',    # Phenylalanine
    'TTA' : 'L',    # Leucine
    'TTG' : 'L',    # Leucine
    'TAC' : 'Y',    # Tyrosine
    'TAT' : 'Y',    # Tyrosine
    'TAA' : '_',    # Stop
    'TAG' : '_',    # Stop
    'TGC' : 'C',    # Cysteine
    'TGT' : 'C',    # Cysteine
    'TGA' : '_',    # Stop
    'TGG' : 'W',    # Tryptophan
    'CTA' : 'L',    # Leucine
    'CTC' : 'L',    # Leucine
    'CTG' : 'L',    # Leucine
    'CTT' : 'L',    # Leucine
    'CCA' : 'P',    # Proline
    'CCC' : 'P',    # Proline
    'CCG' : 'P',    # Proline
    'CCT' : 'P',    # Proline
    'CAC' : 'H',    # Histidine
    'CAT' : 'H',    # Histidine
    'CAA' : 'Q',    # Glutamine
    'CAG' : 'Q',    # Glutamine
    'CGA' : 'R',    # Arginine
    'CGC' : 'R',    # Arginine
    'CGG' : 'R',    # Arginine
    'CGT' : 'R',    # Arginine
    'ATA' : 'I',    # Isoleucine
    'ATC' : 'I',    # Isoleucine
    'ATT' : 'I',    # Isoleucine
    'ATG' : 'M',    # Methionine
    'ACA' : 'T',    # Threonine
    'ACC' : 'T',    # Threonine
    'ACG' : 'T',    # Threonine
    'ACT' : 'T',    # Threonine
    'AAC' : 'N',    # Asparagine
    'AAT' : 'N',    # Asparagine
    'AAA' : 'K',    # Lysine
    'AAG' : 'K',    # Lysine
    'AGC' : 'S',    # Serine
    'AGT' : 'S',    # Serine
    'AGA' : 'R',    # Arginine
    'AGG' : 'R',    # Arginine
    'GTA' : 'V',    # Valine
    'GTC' : 'V',    # Valine
    'GTG' : 'V',    # Valine
    'GTT' : 'V',    # Valine
    'GCA' : 'A',    # Alanine
    'GCC' : 'A',    # Alanine
    'GCG' : 'A',    # Alanine
    'GCT' : 'A',    # Alanine
    'GAC' : 'D',    # Aspartic Acid
    'GAT' : 'D',    # Aspartic Acid
    'GAA' : 'E',    # Glutamic Acid
    'GAG' : 'E',    # Glutamic Acid
    'GGA' : 'G',    # Glycine
    'GGC' : 'G',    # Glycine
    'GGG' : 'G',    # Glycine
    'GGT' : 'G'  }


peptides = { "Strep-tag II": "WSHPQFEK" ,
             "c-myc": "EQKLISEEDL"
            # "S-tag":"KETAAAKFERQHMDS",
            # "HAT-tag":"KDHLIHNVHKEFHAHAHNK"
            }

result_file = open("result_p.txt","w")

for i in peptides:
    AA_Seq = peptides.get(i,"Null")
    print (i + " all the possible DNA sequences of " + AA_Seq + " is: ",file=result_file)
    AL = []
    #x = 1
    for a in AA_Seq:
        print(a, file = result_file)
        PA = []
        for c in codons:
            if codons[c] == a:
                PA.append(c)
        AL.append(PA)
    #print (AL)
    for j in list(product(*AL)):
        
        print ("".join(j),file=result_file)
    print("\n", file=result_file)

result_file.close()


"""
    ref:   http://www.bioinformatics.org/sms2/rev_trans.html
           http://matrix6ro.blog.51cto.com/1746429/1773837
           http://stackoverflow.com/questions/798854/all-combinations-of-a-list-of-lists
"""
