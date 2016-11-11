#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import print_function
from collections import defaultdict
from os import path
from datetime import datetime
from glob import glob
import os
import subprocess
import argparse
import multiprocessing

#OUT = subprocess.Popen("ls -l > list_Of_files.txt", shell = True, stdout = subprocess.PIPE)

#print (OUT) 

PEAR_D = "/home/yangjy/software/pear-0.9.6-bin-64/pear-0.9.6-bin-64"
bwa_bin = "/home/yangjy/software/bwa-0.7.10/bwa"
bwa_ind = "/home/qzy/GCF_000826755.1_ZizJuj_1.1_genomic.fna"

d = os.getcwd()

if not os.path.exists("trimmed_Sample.files"):
    subprocess.call('find ' +d.strip()+'/ -name "*trimed_adaptor_trimed_*.fastq" > trimmed_Sample.files', shell = True)

file1 = open("trimmed_Sample.files", "r")

sample_n = defaultdict(list)
## organize the filenames in the file "trimmed_Sample.files" and store them in the dict sample_n
for line in file1:
    key = line.rstrip()
    if key != "":
        fn = key.split("/")[-1]
        
        name = fn.split("_",1)[0]
        #print (name)
    sample_n[name].append(key)      # a dict contain sample name and its corresponding fq files

file1.close()

for k,v in sample_n.iteritems():    
    v =sorted(v,  key = lambda x: x.rstrip(".fastq").split("_")[-1])
    sample_n[k] = v

def Merge_PE_Reads(FQ_arg):
    
    Run_R = FQ_arg.split("\t")
    process_sample_name = Run_R[0].split("/")[-1].split("_")[0]
    shell_Name = "Pear_M_" + process_sample_name + ".sh"
    MerGe_Bash = open(shell_Name, "w")
    MerGe_Bash.write("#!/bin/bash\n")
    MerGe_Bash.write(PEAR_D +" -f " + Run_R[0] + " -r " + Run_R[1] +" -v 15 -n 80 -e -o pear_merge_" + process_sample_name)
    
    MerGe_Bash.write("\n")
    MerGe_Bash.close()
    return shell_Name    

def Aln_With_BWAMem(FQ_arg):
    Run_R = FQ_arg.split("\t")
    process_sample_name = Run_R[0].split("/")[-1].split("_")[0]
    
    shell_Aname = "Aln_BwaMem_"+process_sample_name+".sh"
    bash_file = open(shell_Aname, "w")
    bash_file.write("#!/bin/bash\n\n")
    # ./bwa mem -M -t 16  ref.fa read1.fq read2.fq | gzip -3 > aln-pe.sam.gz
    bash_file.write("if [[ ! -e Aln_PE_%s.sam.gz ]]\nthen\n    "%(process_sample_name) )
    cmd_fq_to_fa1 = bwa_bin + " mem -M -t 16 " + bwa_ind + " " + Run_R[0]  + " " + Run_R[1] + " | gzip -3 > Aln_PE_" + process_sample_name +".sam.gz"
    #print cmd_fq_to_fa1
    bash_file.write(cmd_fq_to_fa1 + "\n")
    bash_file.write("fi")
    bash_file.close()
    return shell_Aname


def Run_Bash(f):
    command = "bash " +f
    LOG = os.popen(command)
    each_log = open( f + f.split("_")[0] +"_log.txt", "w")
    each_log.write(LOG.read())
    each_log.close()


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()                 #simplifys the wording of using argparse as stated in the python tutorial
    parser.add_argument("-d", type=str, action='store', required= True,  dest='Bash_dir', help="input the forward read file") # a directory store bash_file scripts for bwa align commands
    parser.add_argument("-p", type = int, action='store',  dest='P_Num', default = 8, help="input the forward read file")
    parser.add_argument("-o", "--output", default = "Log_Curr.txt", help="Directs the output to a logfile of your pear assembly results")
    args = parser.parse_args()

    if not (path.exists(args.Bash_dir) and path.isdir(args.Bash_dir)):
        os.makedirs(args.Bash_dir)

    Shell_A_List =[]
    ## create bash file for processing
    os.chdir(args.Bash_dir)
    for i in sample_n:
        PE_fq = "\t".join(sample_n[i])

        shell_n = Aln_With_BWAMem(PE_fq)   
        Shell_A_List.append(shell_n)
    ## execute these bashfile with Multi Processing
    work_pool = multiprocessing.Pool(processes = args.P_Num)
    for i in Shell_A_List:         # previous use sample_n list
        #fq = "\t".join(sample_n[i])
        
        work_pool.apply_async(Run_Bash, (i, ))
        #work_pool.apply_async(Merge_PE_Reads, (fq, ))
        
    print (datetime.now().strftime("%Y-%m-%d__%H_%M_%S"))
    print ("Data analysis start.")
    work_pool.close()
    work_pool.join()
    
    print (datetime.now().strftime("%Y-%m-%d__%H_%M_%S"))
    print ("Analysis Done.")

    jobTime = datetime.now().strftime("%Y-%m-%d__%H_%M_%S")
    if os.path.exists(args.output):
        LogName = args.output + jobTime + ".txt"
    else:
        LogName = args.output

    log_all = open(LogName , "w")
    Pear_proc_Results = glob("*_Pear_merge_log.txt")
    for i in Pear_proc_Results:
        with open(i, "r") as each_Res:
            log_all.write(each_Res.read())
    log_all.close()
