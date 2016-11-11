#!/usr/bin/env python
# -*- coding: utf-8 -*-

from os import path
from glob import glob
from datetime import datetime
import os
import multiprocessing

if not (path.exists("tmp_Bash_Fa_Stat") and path.isdir("tmp_Bash_Fa_Stat")):
    os.makedirs("tmp_Bash_Fa_Stat")

file1 = open("GEN160463_trimmed_Sample.files", "r")

sampleN_L = set()

## organize the Sample names in the file "trimmed_Sample.files" and store them in list
for line in file1:
    key = line.rstrip()
    if key != "":
        samn = key.split("/",2)[1]
        #sampleN_L.append(samn)
        sampleN_L.add(samn)
#print len(sampleN_L)
#print (sampleN_L)

os.chdir("tmp_Bash_Fa_Stat")
 
for i in sampleN_L:
    bash_file = open("fa_Unique_"+i+".bash", "w")
    bash_file.write("#!/bin/bash\n\n")
    bash_file.write("if [[ ! -e sample_%s.fa ]]\nthen\n    "%(i) )
    cmd_fq_to_fa1 = "/home/qzy/software/bin/fastq_to_fasta -i ../pear_merge_" + i + ".assembled.fastq -o " + "out_assembled_" + i +".fasta -Q 33"
    #print cmd_fq_to_fa1
    bash_file.write(cmd_fq_to_fa1 + "\n")
    cmd_fq_to_fa2 = "/home/qzy/software/bin/fastq_to_fasta -i ../pear_merge_" + i + ".unassembled.forward.fastq -o " + "out_unassembled.forward_" + i +".fasta -Q 33"
    #print cmd_fq_to_fa2
    bash_file.write(cmd_fq_to_fa2 + "\n")
    cmd_fq_to_fa3 = "/home/qzy/software/bin/fastq_to_fasta -i ../pear_merge_" + i + ".unassembled.reverse.fastq -o " + "out_unassembled.reverse_" + i +".fasta -Q 33"
    bash_file.write(cmd_fq_to_fa3 + "\n")
    combine_fa_cmd = "cat out_assembled_%s.fasta out_unassembled.forward_%s.fasta out_unassembled.reverse_%s.fasta > sample_%s.fa"%(i,i,i,i)
    #print cmd_fq_to_fa3
    #print combine_fa_cmd
    bash_file.write(combine_fa_cmd + "\n")
    bash_file.write("fi\n")
    bash_file.write("if [[ -e sample_%s.fa ]] && [[ ! -e Unique_Sample_%s_res.fa ]]\nthen\n    "%(i,i)) 
    CD_HIT_UniCMD = "/home/yangjy/software/cd-hit-v4.5.4-2011-03-07/cd-hit-est -M 6000 -i sample_%s.fa -o Unique_Sample_%s_res.fa -c 1.0"%(i,i)
    #print CD_HIT_UniCMD
    bash_file.write(CD_HIT_UniCMD + "\n")
    bash_file.write("fi")
    bash_file.close()    
def get_uniq_seq(f):
    Run_Shell_ample = "fa_Unique_"+f+".bash"
    log_c = os.popen("bash "+ Run_Shell_ample)
    each_log = open( f + "_CDHIT_log.txt", "w")
    each_log.write(log_c.read())
    each_log.close()

if __name__ == "__main__":
    work_pool_A = multiprocessing.Pool(processes = 16)
    for i in sampleN_L:
        work_pool_A.apply_async(get_uniq_seq, (i, ))
    
    work_pool_A.close()
    work_pool_A.join()

    print (datetime.now().strftime("%Y-%m-%d-%H:%M:%S"))
    print ("Analysis Done.")

    log_all = open("CD_HITLOG.txt" , "w")
    CDHIT_proc_Results = glob("*_CDHIT_log.txt")
    for i in CDHIT_proc_Results:
        with open(i, "r") as each_Res:
            log_all.write(each_Res.read())
    log_all.close()
