#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import print_function
from collections import defaultdict
from os import path
from datetime import datetime
from glob import glob
import os, sys
import subprocess
import argparse
import multiprocessing
import shlex
import shutil
import re

# ##  define the source directory of software 

PEAR_D = "/home/yangjy/software/pear-0.9.6-bin-64/pear-0.9.6-bin-64"
bwa_bin = "/home/yangjy/software/bwa-0.7.10/bwa"
samtools_s = "/home/yangjy/software/samtools-0.1.18/samtools"
bam2bed = "/usr/bin/bamToBed"

def Pre_build(genome_fa):
    pre_name = ".".join(genome_fa.split(".")[:-1])
    genome_fa_new = pre_name + ".fasta"
    if not os.path.exists(genome_fa_new):
        os.rename(genome_fa, genome_fa_new)
    if not os.path.exists(genome_fa_new + ".bwt"):
        os.system(bwa_bin + " index -a bwtsw "+ genome_fa_new)
    if not os.path.exists(genome_fa_new + ".fai"):
        os.system(samtools_s + " faidx " + genome_fa_new)
    if not os.path.exists(pre_name + ".dict"):
        CreateDict_command = """java -jar /home/yangjy/software/picard-tools-1.86/CreateSequenceDictionary.jar REFERENCE={0} OUTPUT={1}""".format(genome_fa_new, genome_fa_new.strip("fasta")+"dict")
        os.system(CreateDict_command)
    return genome_fa_new

def get_sample(directory):
    if not os.path.exists("trimmed_fastq_files.txt"):
        subprocess.call('find ' + directory.strip()+'/ -name "*trimed_adaptor_trimed_*.fastq" > trimmed_fastq_files.txt', shell = True)

    file1 = open("trimmed_fastq_files.txt", "r")
    sample_f = defaultdict(list)
    ## organize the filenames in the file "trimmed_Sample.files" and store them in the dict sample_n
    for line in file1:
        key = line.rstrip()
        if key != "":
            fn = key.split("/")[-1]
            name = fn.split("_",1)[0]
        #print (name)
        sample_f[name].append(key)      # a dict contain sample name and its corresponding fq files
    file1.close()

    for k,v in sample_f.iteritems():    
        v =sorted(v,  key = lambda x: x.rstrip(".fastq").split("_")[-1])
        sample_f[k] = v
    return sample_f

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


def Get_Tags(i, dir_p):
    shell_n = "fa_Unique_"+i+".sh"
    bash_file = open(shell_n, "w")
    bash_file.write("#!/bin/bash\n\n")
    bash_file.write("if [[ ! -e sample_%s.fa ]]\nthen\n    "%( i) )
    cmd_fq_to_fa1 = "/home/qzy/software/bin/fastq_to_fasta -i " + dir_p +"/pear_merge_" + i + ".assembled.fastq -o " + "out_assembled_" + i +".fasta -Q 33"
    #print cmd_fq_to_fa1
    bash_file.write(cmd_fq_to_fa1 + "\n")
    cmd_fq_to_fa2 = "/home/qzy/software/bin/fastq_to_fasta -i " + dir_p +"/pear_merge_" + i + ".unassembled.forward.fastq -o " + "out_unassembled.forward_" + i +".fasta -Q 33"
    #print cmd_fq_to_fa2
    bash_file.write(cmd_fq_to_fa2 + "\n")
    cmd_fq_to_fa3 = "/home/qzy/software/bin/fastq_to_fasta -i " + dir_p +"/pear_merge_" + i + ".unassembled.reverse.fastq -o " + "out_unassembled.reverse_" + i +".fasta -Q 33"
    bash_file.write(cmd_fq_to_fa3 + "\n")
    combine_fa_cmd = "cat out_assembled_%s.fasta out_unassembled.forward_%s.fasta out_unassembled.reverse_%s.fasta > sample_%s.fa"%(i,i,i,i)
    bash_file.write(combine_fa_cmd + "\nfi\n")
    command_mkdir_mv = "mkdir %s_fa\nmv sample_%s.fa %s_fa/\ncd %s_fa/"%(i,i,i,i)
    bash_file.write(command_mkdir_mv + "\n")
    bash_file.write("if [[ -e sample_%s.fa ]] && [[ ! -e sample_%s_sorted.fasta ]]\nthen\n    "%(i,i)) 
    Uclust_UniCMD = """uclust --mergesort sample_%s.fa --output sample_%s_sorted.fasta\nfi\nif [[ ! -e results_sim_%s_tmp.uc ]]\nthen\n
uclust --input sample_%s_sorted.fasta --uc results_sim_%s_tmp.uc --id 0.98 --rev --nucleo\nfi\nif [[ ! -e ../result_%s_log.txt ]]\nthen\n
echo -ne \"Tag number count\\tresults_sim_%s_tmp.uc\\t\" > ../result_%s_log.txt\n
echo $(grep -c \"^S\" results_sim_%s_tmp.uc) >> ../result_%s_log.txt\n
#uclust --uc2fasta results_sim_%s_tmp.uc --input sample_%s.fa --output Unique_uclust_%s_res.fa --types S\n
"""%(i,i,i,i,i,i,i,i,i,i,i,i,i)
    bash_file.write(Uclust_UniCMD + "\n")
    bash_file.write("fi\ncd ..\n")
    bash_file.close()    
    return shell_n

def Aln_With_BWAMem(FQ_arg, genome):
    Run_R = FQ_arg.split("\t")
    Psn = Run_R[0].split("/")[-1].split("_")[0]    # Psn is the Abbr. of process_sample_name    
    shell_Aname = "Aln_BwaMem_"+ Psn +".sh"
    bash_file = open(shell_Aname, "w")
    bash_file.write("#!/bin/bash\n\n")
    #bash_file.write("if [[ ! -e Aln_PE_%s.sam.gz ]]\nthen\n    "%(process_sample_name) )
    # ./bwa mem -M -t 16  ref.fa read1.fq read2.fq | gzip -3 > aln-pe.sam.gz
    bash_file.write("if [[ ! -e Aln_PE_%s_sorted.bam ]]\nthen\n"%(Psn) )
    #cmd_fq_to_sam1 = bwa_bin + " mem -M -t 16 " + bwa_ind + " " + Run_R[0]  + " " + Run_R[1] + " | gzip -3 > Aln_PE_" + process_sample_name +".sam.gz"
    cmd_fq_to_sam = bwa_bin + " mem -M -t 16 " + genome + " " + Run_R[0]  + " " + Run_R[1] + " > Aln_PE_" + Psn +".sam 2> ../Aln_Err_"+Psn +"_log.txt"
    #print cmd_fq_to_fa1
    bash_file.write("    " + cmd_fq_to_sam + "\n")
    cmd_sam_to_SORTBam = "%s view -bS Aln_PE_%s.sam | %s sort - Aln_PE_%s_sorted\nrm Aln_PE_%s.sam\n"%(samtools_s, Psn, samtools_s, Psn, Psn)
    bash_file.write("    " + cmd_sam_to_SORTBam + "\nfi\n")
    bash_file.write("""if [[ ! -e Aln_PE_%s.bed ]]\nthen\n
    %s -i Aln_PE_%s_sorted.bam > Aln_PE_%s.bed\nfi\n
if [[ ! -e Aln_PE_%s_sorted.bam.bai ]]\nthen\n    %s index Aln_PE_%s_sorted.bam\n"""%(Psn, bam2bed, Psn, Psn, Psn, samtools_s, Psn))
    bash_file.write("fi\n")
    group_n = Psn + "_sorted"
    bash_file.write("if [[ ! -e %s_RG.bam ]]\nthen\n    java -jar /home/yangjy/software/picard-tools-1.86/AddOrReplaceReadGroups.jar INPUT=Aln_PE_%s_sorted.bam OUTPUT=%s_RG.bam RGID=%s RGLB=%s RGPL='Illumina' RGSM=%s RGPU=%s VALIDATION_STRINGENCY=SILENT\nfi\n"%(Psn, Psn, Psn, Psn, group_n, group_n, group_n))
    bash_file.write("if [[ ! -d \"%s_bamqc\" ]]\nthen\n    mkdir %s_bamqc\nfi\nif [[ ! -e \"%s_bamqc/qualimapReport.html\" ]]\nthen\n    qualimap bamqc -bam %s_RG.bam -gff %s.bed -outdir %s_bamqc\nfi\n"%(Psn, Psn,Psn, Psn, genome.strip(".fastq"),Psn ))
    bash_file.write("if [[ ! -e \"%s_bamqc/qualimap_stat_%s.txt\" ]]\nthen\n php /home/qzy/software/qualiHTMLParse.php %s\nfi\n"%(Psn, Psn, Psn))
    bash_file.close()
    return shell_Aname

def depth_and_MapRatio(each_bam):
    Gen_Cov = subprocess.Popen('bedtools genomecov -ibam {0} -bg'.format(each_bam), shell = True, stdout = subprocess.PIPE )
    bed_Depth = open(each_bam + "_Depth.txt", 'w')
    bed_Depth.write(each_bam + "Bed Coverage Graph.\n")
    bed_Depth.write(Gen_Cov.stdout.read())
    bed_Depth.close()
    map_sta = subprocess.Popen('samtools flagstat {0}'.format(each_bam), shell = True, stdout = subprocess.PIPE )
    Ave_cmd = "samtools depth %s | awk 'BEGIN{sum=0}{sum+=$3}END{print sum/NR}'" % (each_bam)
    Ave_dep = subprocess.Popen(Ave_cmd, shell=True, stdout = subprocess.PIPE)
    stat_res = open(each_bam+"samtools_Depth_Result.txt", 'w')
    stat_res.write(map_sta.stdout.read())
    stat_res.write(each_bam+ " Depth:\t")
    stat_res.write(Ave_dep.stdout.read())
    stat_res.close()

def Run_Bash(f):
    each_log = f + f.split("_")[0] +"_log.txt"
    command = "bash " +f + " 1> "+ each_log + " 2>&1"
    RUN = os.popen(command)

def Variant_Calling(genome, out_prefix):
    result_dir = os.getcwd()
    CS_SH = open("current_vcf_call.sh", 'w')
    CS_SH.write("#!/bin/bash\n\n")
    #CS_SH.write('java -Xmx15g -jar /home/yangjy/software/picard-tools-1.86/MergeSamFiles.jar OUTPUT=Merge.bam $(printf "INPUT=%s " RG_*.bam)\n')
    CS_SH.write("for i in `ls *RG.bam`\ndo\n     ID=$(echo $i | cut -d '_' -f1)\n     echo -e \"@RG\tID:${ID}_RG\tPL:Illumina\tPU:${ID}sorted\tLB:${ID}_sorted\tSM:${ID}_sorted\" >> tt.txt\ndone\n")
    CS_SH.write("sort tt.txt | uniq > rh.txt\n")
    CS_SH.write('samtools merge -r -h rh.txt Merged_All.bam *_RG.bam\n')
    CS_SH.write("samtools mpileup -uDgf %s %s/Merged_All.bam | /home/yangjy/software/samtools-0.1.18/bcftools/bcftools view -bvcg - > %s/%s.bcf\n\n"%(genome, result_dir, result_dir,out_prefix))
    CS_SH.write("/home/yangjy/software/samtools-0.1.18/bcftools/bcftools view %s/%s.bcf | /home/yangjy/software/samtools-0.1.18/bcftools/vcfutils.pl varFilter -D 1200 > %s/%s.vcf\n"%(result_dir, out_prefix, result_dir, out_prefix))
#
    CS_SH.write("/usr/bin/bamToBed -i Merged_All.bam > Merged_All.bed")
    CS_SH.close()
    os.popen("bash current_vcf_call.sh 1> samtools_vcf_log 2 > samtools_err.txt")

def replace_M(pre_name):
    a = glob(pre_name+"*log.txt")
    for i in a:
        with open(i, 'r')as a_s, open("new_"+i, 'w')as a_r:
            a_c = re.sub(r"[\x00-\x1f]", "\n", a_s.read())
            a_r.write(a_c)

def get_fq_list():
    file_list = defaultdict(list)
    for root,dirs,files in os.walk("."):
        for f in files:
            if f.endswith("fastq") and f.find("-") != -1 and f.find("_") != -1:
                xs = f.rindex("-")
                xe = f.index("_", xs)
                s_name = f[xs+1:xe]
                file_list[s_name].append(f)
    if list(file_list.keys()) == []:
        print("no fastq found, please check this directory.")
        return False,None
    return True,file_list

if __name__ == "__main__":

    USAGE = """Analysis stage number selection:\n\n
     1. QC and adaptor removing\n\n
     2. Pear merge\n\n
     3. Clustering with uclust\n\n
     4. bwa Mapping\n\n
     5. depth counting\n\n
     6. SNP calling\n\n"""

    parser = argparse.ArgumentParser( description = USAGE )

    parser.add_argument("-s", action='store_true', default=False, help="skip the procedure of genome library building")
    parser.add_argument("-a", type = int, help="analysis phase",default= [0], nargs="+", dest = "stage", choices = [1, 2, 3, 4, 5, 6])
    parser.add_argument("-i", type=str, action='store', dest='gen_f' , help="file name of genome") 
    parser.add_argument("-d", type=str, action='store', dest='Direct', default="." , help="input the dir of trimmed fastq file") 
    parser.add_argument("-t", type=str, action='store',  dest='Bash_dir',default="temp", help="a directory store bash_file scripts for bwa align commands")   # this directory will store pear merged results, uclust results, aln results and depth results.
    parser.add_argument("-p", type = int, action='store',  dest='P_Num', default = 16, help="input the forward read file")
    parser.add_argument("-o", "--output", default = "map_stat_COV", help="Directs the output to a logfile of your pear assembly results")
    args = parser.parse_args()

    if not (path.exists(args.Bash_dir) and path.isdir(args.Bash_dir)):
        os.makedirs(args.Bash_dir)
    
    absolute_p = os.path.abspath(args.Direct)
    sample_n = get_sample(absolute_p)

    if args.s == True:
        f_genome = os.popen("cat genome_file.txt")
        genome_name = f_genome.read().strip()
    else:
        if args.gen_f == None:
            sys.exit("Please provide species genome fasta file.")
        else:    
            genome_name = Pre_build(args.gen_f)
            gen_tmp = open("genome_file.txt",'w')
            print(genome_name, file=gen_tmp)
            gen_tmp.close()

    os.chdir(args.Bash_dir)
    proce_L = args.stage

    pear_m_dir = absolute_p +"/" + args.Bash_dir+ "/Pear_Process"

    Shell_Merge_List =[]              # a list store shell name of merge command with PEAR software
    Shell_get_tag = []                # a list store shell name of Uclust commands
    Shell_A_List =[]                  # a list store shell name of align scripts
  
    ## create bash file for processing
    if not path.exists("Pear_Process"):
        os.makedirs("Pear_Process")
        os.chdir("Pear_Process")
        os.makedirs("Uclust_unique")
        os.chdir("../")
    else:
        os.chdir("Pear_Process")
        if not path.exists("Uclust_unique"):
            os.makedirs("Uclust_unique")
        os.chdir("../")

    if not path.exists("BWA_Mapping"):
        os.makedirs("BWA_Mapping")    
    
    if not path.exists("Variants"):
        os.makedirs("Variants")

    for i in sample_n:
        PE_fq = "\t".join(sample_n[i])
        os.chdir("Pear_Process")
        shell_m = Merge_PE_Reads(PE_fq)
        os.chdir("Uclust_unique")
        Uclust_shell_n = Get_Tags(i, pear_m_dir)
        os.chdir("../../")
        os.chdir("BWA_Mapping")
        shell_n = Aln_With_BWAMem(PE_fq, genome_name)
        os.chdir("../")   
        Shell_Merge_List.append(shell_m)
        Shell_A_List.append(shell_n)
        Shell_get_tag.append(Uclust_shell_n)        # ## get Uclust shell script name and store them in a list
    
    ## check if user want to run all procedures

    if proce_L == [0]:
        print("You have select all procedures go ? Are you sure? please input Y/N\n")
        Ask = ""
        while (Ask not in ["Y", "y", "N", "n"]):
            Ask = raw_input("Please input Y or N:\n")
            if Ask == "Y" or Ask == "y":
                proce_L.extend([1,2,3,4,5,6])
                break
            if Ask == "N" or Ask == "n":
                sys.exit("You have skip all procedures.")
    else:
        pass

    ## get sample list and run QC

    if 1 in proce_L and not os.path.exists("../sample_list.txt"):
        any_fq, fqs = get_fq_list()
        if not any_fq:
            sys.exit("Invalid path.")
        fq_sample_h = open("sample_list.txt", 'w')
        for i,v in fqs.items():
            if len(v) < 2:
                print("the sample {0} does not have paired fq files, please check it.".format(i))
            else:
                v = sorted(v, key=lambda x: x.rstrip(".fastq").split("_")[-1])
                fq_sample_h.write("{0}\t{1}\t{2}\t{3}".format(i,i,v[0],v[1]))
        fq_sample_h.close()
        QC_command = "perl /home/yangjy/software/Script/RNASeq_Script/v1.0/1_RNASeq_QC.pl -i sample_list.txt -p_num 16 -num 10000 > QC.log 2>&1"
        QC_Run = subprocess.Popen( QC_command , shell=True , stdout=subprocess.PIPE )
        QC_Run.wait()
        #os.popen("perl /home/yangjy/software/Script/RNASeq_Script/v1.0/1_RNASeq_QC.pl -i sample_list.txt -p_num 16 -num 10000 > QC.log 2>&1")
    elif os.path.exists("../sample_list.txt"):
        print("The file \"sample_list.txt\" found in corresponding directory. You could check the QC_results in \"1_QC\".")
    else:
        print("You may skip the QC step for that you have finished it.")

    ## Run the merge commands for paired fastq files
    if 2 in proce_L and not os.path.exists("merge_finished.txt"):
        os.chdir("Pear_Process")
        work_pool_P = multiprocessing.Pool(processes = args.P_Num)
        for i in Shell_Merge_List:
            work_pool_P.apply_async(Run_Bash, (i, ))

        print (datetime.now().strftime("%Y-%m-%d  %H:%M:%S"))
        print ("Merging fastq Data with PEAR analysis start.")
        work_pool_P.close()
        work_pool_P.join()
        STD_Merge_L = os.popen("awk '/Assembled reads/' Pear_M*_log.txt")
        meg_log = open("merge_finished.txt",'w')
        print(STD_Merge_L.read().strip(), file=meg_log)
        meg_log.close()
        shutil.copy("merge_finished.txt", "../")
        os.chdir("../")
    
    ## run Uclust shell script in the List Shell_get_tag 
    if 3 in proce_L and (not os.path.exists("clust_finished.txt") or os.path.getsize("clust_finished.txt") ==0):
        os.chdir("Pear_Process/Uclust_unique")
        work_pool_U = multiprocessing.Pool(processes = args.P_Num)
        for i in Shell_get_tag:
            work_pool_U.apply_async(Run_Bash, (i,))
        work_pool_U.close()
        work_pool_U.join()
        #replace_M(fa_Unique)
        print (os.getcwd())    
        STD_Clust = os.popen("awk '/Tag number count/' result*log.txt")
        cls_log = open("clust_finished.txt", 'w')
        print(STD_Clust.read().strip(), file=cls_log)
        cls_log.close()
        shutil.copy("clust_finished.txt", '../../')
        os.chdir("../../")

    ## execute these bashfile wiht bwa -mem commands with Multi Processing
    if 4 in proce_L and not os.path.exists("OUT_sortbam.log"):
        os.chdir("BWA_Mapping")
        work_pool = multiprocessing.Pool(processes = args.P_Num)
        for i in Shell_A_List:         # previous use sample_n list
            #fq = "\t".join(sample_n[i])
            work_pool.apply_async(Run_Bash, (i, ))
            #work_pool.apply_async(Merge_PE_Reads, (fq, ))
        work_pool.close()
        work_pool.join()
    
        print (datetime.now().strftime("%Y-%m-%d__%H_%M_%S"))
        print ("Align fastq and sort sam Analysis Done.")
        os.chdir("../")
        log_o_bam = open("OUT_sortbam.log" , "w")

        BAM_log_Results = glob("Aln_Err_*log.txt")
        if BAM_log_Results != []:
            for i in BAM_log_Results:
                with open(i, "r") as each_Res:
                    log_o_bam.write(each_Res.read())
                if os.path.exists(i):
                    os.remove(i)
        log_o_bam.close()
        os.chdir("../")

    ## define the file name which will store the mapping ratio result and coverage result

    LogName = args.output+".txt"

    ## parallel Processing mapping ratio and coverage computation
    if 5 in proce_L and (not os.path.exists(LogName) or os.path.getsize(LogName) ==0):
        os.chdir("BWA_Mapping")
        bamfiles = glob("RG*sorted.bam")
        work_pool_M = multiprocessing.Pool(processes = args.P_Num)

        for i in bamfiles:    
            work_pool_M.apply_async(depth_and_MapRatio, (i,))
    
        work_pool_M.close()
        work_pool_M.join()

        StdCov = os.popen("awk '/0 mapped/' *Depth_Result.txt")
        os.chdir("../")
        map_results = open(LogName, "w")
        map_results.write(StdCov.read())
        map_results.close()
    
    if 6 in proce_L:
        os.chdir("BWA_Mapping")
        Variant_Calling(genome_name, "A_test" )


