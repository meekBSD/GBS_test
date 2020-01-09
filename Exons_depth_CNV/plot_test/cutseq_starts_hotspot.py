#!/usr/bin/env python
# -*- coding: UTF-8 -*-


## pattern like follows:

## "CA"ACC...TG"AAC"
##
## "AAC"CGC.."CA"

## reference : Predisposition to childhood acute lymphoblastic leukemia caused by a constitutional translocation disrupting ETV6


from Bio import SeqIO

def get_acc(identifier):
    """"Given a SeqRecord identifier string, return the accession number as a string.

    e.g. "gi|2765613|emb|Z78488.1|PTZ78488" -> "Z78488.1"
    """
    parts = identifier.split(" ")
    #assert len(parts) == 5 and parts[0] == "gi" and parts[2] == "emb"
    return parts[0]


ETV6_seqs = SeqIO.index("ETV6_intron5.fa", "fasta", key_function=get_acc)
ETV6_i5 = ETV6_seqs["ETV6_intron5"].seq

NTRK3_seqs = SeqIO.index("NTRK3_intron14.fa", "fasta", key_function=get_acc)
NTRK3_i14 = NTRK3_seqs["NTRK3_intron14"].seq
n14 = NTRK3_i14.reverse_complement()

for i in range(len(ETV6_i5)-20):
    e_seq = ETV6_i5[i:i+20]
    print(e_seq)
    for j in range(7,21):
        for n in range(len(n14)-j):
            n_seq = n14[n:n+j]

            if n_seq[:5] == e_seq[-5:] and e_seq.find(n_seq[-6:]) != -1:
                print("{0}\t{1}".format(i, e_seq))
                print("{0}\t{1}\t{2}".format(j, n, n_seq))


