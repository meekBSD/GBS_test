#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from collections import defaultdict
from glob import glob
from random import randint
import numpy as np
import argparse
import os
import sys
import vcf
import copy
import math
import subprocess
import yaml


def load_config(config_file): 
    """Load YAML config file, replacing environmental variables. 
    """ 
    with open(config_file) as in_handle: 
        config = yaml.load(in_handle) 
    #    for field, setting in config.items(): 
    #        if isinstance(config[field], dict): 
    #            for sub_field, sub_setting in config[field].items(): 
    #                config[field][sub_field] = expand_path(sub_setting) 
    #        else: 
    #            config[field] = expand_path(setting) 
    return config 

def coverage_stats(f, max, interval):     # also can set parameter max = 500 and interval = 50, but it can't be used in argparse conveniently.
    depth = {}
    for i in range(50, max, interval):
        depth[i] = 0
    depth[10] = 0

    all_sample_dep = {}
    
    AD = list(depth.keys())
    AD.sort()
    a = vcf.Reader(filename = f)
    for samp in a.samples:
        all_sample_dep[samp] = copy.deepcopy(depth)
    
    for record in a:
        for s in record.samples:
            s_n = s.sample
            Dep_curr = s['DP']
            for n,num_r in enumerate(AD[:-1]):                  # num_r is the bottom number range margin                
                if Dep_curr >= num_r and Dep_curr < AD[n+1]:    # AD[n+1] is top number margin
                    all_sample_dep[s_n][num_r] += 1
                else:
                    continue
    return (depth, all_sample_dep)

# Abbreviation of total depth is TotDP, Abbreviation of filter standard of Quality is filter_st

def filter_SNP(vcf_file, TotDP, filter_st, prefix):

    vcf_reader = open(vcf_file, 'r')
    vcf_W = open(prefix + "_out_filtered_var.vcf", 'w')
    filtered_SNP = open(prefix + "_output_SNP.txt", 'w')
    # filtered_INDEL = open(prefix + "_output_IND.txt" ,'w')
    filtered_SNP.write("Chr\tPOS\tRec_ID\tREF\tALT\tQUAL\tPMV\tMAF\tHET\tHet_Obs\n")
    # filtered_INDEL.write("Chr\tPOS\tRec_ID\tREF\tALT\tQUAL\tPMV\tMAF\tHET\tHet_Obs\n")

    Var_Num = 0
    SNP_Num = 0

    PMV_count = {
        80: 0 ,
        70: 0 ,
        60: 0 ,
        40: 0 ,
        20: 0 ,
        10: 0 ,
         5: 0 ,
         1: 0
    }

    for record_v in vcf_reader:
        if record_v.startswith("#"):
            vcf_W.write(record_v)
            
        else:
            fields = record_v.rstrip().split("\t")
            C = len(fields) - 9
            INFO = fields[7].split(";")
            v_QUAL = float(fields[5])

            #v_CHROM = fields[0]
            #v_POS = fields[1]
            #v_ID = fields[2]
            #v_REF = fields[]
            if not INFO[0].startswith("INDEL") and INFO[2].startswith("AF1"):
                af = float(INFO[2].split("=")[1])
                INFO_DP = int(INFO[0].split("=")[1])
                DP_A = [i.split(":")[2] for i in fields[9:]]
                GT_A = [i.split(":")[0] for i in fields[9:]]

                if v_QUAL > filter_st:
                    PMV_r = len([ s for s in DP_A if int(s) > 2 ]) * 100 / C
                    PMV = 100 - PMV_r            
                    for k in PMV_count:
                        if round(PMV_r) > 100 - k:
                            PMV_count[k] += 1

                    Het_Freq = GT_A.count("0/1")/C

                    Var_Num += 1
                    if INFO_DP >= TotDP and af >= 0.5:
                        vcf_W.write(record_v)
            # check this position is a SNP and its AF > 0.5, and then write this line into test_output_SNP.txt
            #    if record_v.is_snp and af >= 0.5:
                        SNP_Num += 1
                        filtered_SNP.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format\
                     ( "\t".join(fields[:6]), PMV,  1-af, 0 , Het_Freq))
            # check this position is a SNP and its AF <= 0.5, and then write this line into test_output_SNP.txt
                    elif INFO_DP >= TotDP and af < 0.5:
                        vcf_W.write(record_v)
                        SNP_Num += 1
                        filtered_SNP.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format\
                     ( "\t".join(fields[:6]), PMV,  af, 0, Het_Freq))
    vcf_W.close()
    filtered_SNP.close()
    #filtered_INDEL.close()

    #return (Var_Num, SNP_Num, IND_Num, PMV_count, Qual_Of_Samples)
    return (Var_Num, SNP_Num, PMV_count)

def data_ready(input_f):
    a = list(range(0,51))
    PMV_D = {
        80: 0 ,
        70: 0 ,
        60: 0 ,
        40: 0 ,
        20: 0 ,
        10: 0 ,
         5: 0 ,
         1: 0}

    for i in PMV_D:
        d = {}
        for j in a:
            d[j] = 0
        PMV_D[i] = d

    maf_d = {}
    for num in a:
        maf_d[num] = 0

    het_d = {}
    for i in range(0,101):
        het_d[i] = 0

    ## Reade SNP filteration result and get maf Het observation and Percent missing value , store them in corresponding dict.

    with open(input_f, 'r') as stat_handle:
        for i in stat_handle:
            if i.startswith("Chr"):
                continue
            else:
                sp_line = i.rstrip().split("\t")
                maf = int(float(sp_line[7]) * 100)   # could use math.floor in replace of int() function
                maf_d[maf] += 1
                pmv = float(sp_line[6])
                het = int( float(sp_line[9]) * 100)
                het_d[het] += 1
                for k in PMV_D:
                    if pmv <= k:
                        PMV_D[k][maf] += 1

    return (PMV_D,  maf_d,  het_d)

def data_process(SFR, window_size, step, bam_f, out_f):
    SNP_handle = open(SFR, 'r')      # SFR is abbreviation name of Records after SNP filtration, tab delimited file
    chromo = defaultdict(dict)

    for i in SNP_handle:
        if i.startswith("Chr"):
            continue
        else:
            splits = i.rstrip().split("\t")
            chr_name = splits[0]
            Position = splits[1]
            steps_in_win = window_size/step
            chr_win_start = (int(int(Position)//step) - steps_in_win +1) * step
            for win in range(chr_win_start, chr_win_start + window_size, step):
                if win >= 0 and win in chromo[chr_name]:
                    chromo[chr_name][win] += 1
                elif win >= 0:
                    chromo[chr_name][win] = 1
    SNP_handle.close()

    chr_info = {}

    chr_reading = os.popen("samtools view -H "+ bam_f)
    for i in chr_reading.read().split("\n"):
        if i.startswith("@SQ"):
            info = i.split("\t")
            ch_n = info[1].split(":")[1]
            ch_l = info[2].split(":")[1]
            chr_info[ch_n] = ch_l

    File_Handle = open(out_f, 'w')

    for i in chr_info:
        chromosome_length = int(chr_info[i])
        if chromosome_length >= window_size:
            for j in range(0, int(chr_info[i]), step):
                if j + window_size <= chromosome_length + step:
                    locus = j + window_size/2
                    File_Handle.write("{0}\t{1}-{2}\t{3}\t{4}\n".format(i, j, j+window_size ,locus,  chromo[i].get(j, 0)))

    File_Handle.close()
'''
#### the output of this function is for the Chromosome.R ####
def chrome_data_prepare(SNP_f, window_size, step, bam_file, chr_list):
    SNP_file = open(SNP_f, 'r')
    chromo = defaultdict(dict)

    for i in SNP_file:
        if i.startswith("Chr"):
            continue
        else:
            splits = i.rstrip().split("\t")
            chr_name = splits[0]
            Position = splits[1]
            steps_in_win = window_size/step
            chr_win_start = (int(int(Position)//step) - steps_in_win +1) * step
            for win in range(chr_win_start, chr_win_start + window_size, step): 
                if win >= 0 and win in chromo[chr_name]:
                    chromo[chr_name][win] += 1
                elif win >= 0:
                    chromo[chr_name][win] = 1
    SNP_file.close()

    chr_info = {}

    chr_reading = os.popen("samtools view -H "+bam_file)
    for i in chr_reading.read().split("\n"):
        if i.startswith("@SQ"):
            info = i.split("\t")
            ch_n = info[1].split(":")[1]
            ch_l = info[2].split(":")[1]
            chr_info[ch_n] = ch_l

    CHR_N = {}
    with open(chr_list, 'r') as CT:
        for i in CT:
            C = i.rstrip().split("\t")
            CHR_N[C[1]] = "chr" + C[0]+ "\t" + chr_info[C[1]]

    sort_chromosomes = sorted(CHR_N.items(), key = lambda x: int(x[1].lstrip("chr").split("\t")[0]))
    return (sort_chromosomes, chromo) 

def prepare_DataFor_Circ(SNP_filter, window, List_F, Spe_Abbr, bam_F, bed_F ):
    chromo = defaultdict(dict)

    CHR_N = {}
    with open(List_F, 'r') as CT:
        for i in CT:
            C = i.rstrip().split("\t")
            CHR_N[C[1]] = Spe_Abbr + C[0]

    maf_c = open("maf_norm.txt", 'w')
    pmv_c = open("pmv_norm.txt", 'w')

    maf_d = defaultdict(list)
    pmv_d = defaultdict(list)

    with open(SNP_filter , 'r') as SNP_summary:
        for i in SNP_summary:
            if i.startswith("Chr"):
                continue
            else:
                s = i.rstrip().split("\t")
                pos = int(s[1])
                chr_r = int( pos / window) * window
                maf_tmp = 1/(float(s[7]) + 0.001)
                mp_key = int(pos/(window/10)) * window/10
                if s[0] in CHR_N:
                    #maf_c.write("{0} {1} {2} {3}\n".format(CHR_N[s[0]], pos, pos , math.log(maf_tmp, 2)))
                    #pmv_c.write("{0} {1} {2} {3}\n".format(CHR_N[s[0]], pos, pos , 100/(float(s[6]) + 1)))
                    cn = CHR_N[s[0]]
                    maf_d[cn+"_"+str(mp_key)].append(math.log(maf_tmp,2))
                    pmv_d[cn+"_"+str(mp_key)].append(100/(float(s[6]) + 1))
                if chr_r in chromo[s[0]]:
                    chromo[s[0]][chr_r] += 1
                else:
                    chromo[s[0]][chr_r] = 1

    for i in sorted(maf_d.items(), key = lambda x: x[0].split("_")[0]):
        cs = i[0].split("_")
        chr = cs[0]
        start = cs[1]
        end = int(start) + window/10 -1
        aver = sum(i[1])/len(i[1])
        maf_c.write("{0} {1} {2} {3}\n".format(chr, start, end , aver))

    for i in sorted(pmv_d.items(), key = lambda x: x[0].split("_")[0]):
        cs = i[0].split("_")
        chr = cs[0]
        start = cs[1]
        end = int(start) + window/10 -1
        aver = sum(i[1])/len(i[1])
        pmv_c.write("{0} {1} {2} {3}\n".format(chr, start, end , aver))

    maf_c.close()
    pmv_c.close()
    chr_info = {}

    #  get chromosome info and create genome file
    chr_reading = os.popen("samtools view -H "+ bam_F)
    for i in chr_reading.read().split("\n"):
        if i.startswith("@SQ"):
            info = i.split("\t")
            ch_n = info[1].split(":")[1]
            ch_l = info[2].split(":")[1]
            chr_info[ch_n] = ch_l

    colors = [["vvlred","vlred","lred","red","dred","vdred","vvdred"],
          ["vvlpred", "vlpred", "lpred", "pred", "dpred", "vdpred", "vvdpred"],
          ["vvlgreen", "vlgreen" ,"lgreen", "green", "dgreen", "vdgreen", "vvdgreen"],
          ["vvlpgreen", "vlpgreen", "lpgreen", "pgreen", "dpgreen", "vdpgreen", "vvdpgreen"],
          ["vvlblue", "vlblue", "lblue", "blue", "dblue","vdblue","vvdblue"],
          ["vvlpblue", "vlpblue", "lpblue", "pblue", "dpblue", "vdpblue", "vvdpblue"],
          ["vvlpurple", "vlpurple", "lpurple", "purple", "dpurple", "vdpurple", "vvdpurple"],
          ["vvlppurple", "vlppurple", "lppurple", "ppurple", "dppurple", "vdppurple", "vvdppurple"],
          ["vvlorange", "vlorange", "lorange", "orange", "dorange", "vdorange", "vvdorange"],
          ["vvlporange", "vlporange", "lporange", "porange", "dporange", "vdporange", "vvdporange"],
          ["vvlyellow", "vlyellow", "lyellow", "yellow", "dyellow", "vdyellow", "vvdyellow"],
          ["vvlpyellow", "vlpyellow", "lpyellow", "pyellow", "dpyellow", "vdpyellow", "vvdpyellow"]]

    genome_C = open("SPE_GENOME.txt", 'w')
    sort_chromosome = sorted(CHR_N.items(), key = lambda x: int(x[1].lstrip(Spe_Abbr)))
    for n, i in enumerate(sort_chromosome):
        genome_C.write("chr - {0} {1} 0 {2} {3}\n".format(i[1], n+1, chr_info[i[0]], colors[randint(0,11)][randint(0,6)]))
    genome_C.close()

    tag_c = defaultdict(dict)
    for i in chr_info:
        CL = int(chr_info[i])   # CL is chromosome Length
        if CL >= window:
            for j in range(0, CL + 1, window):
                tag_c[i][j] = 0

    with open(bed_F , 'r') as Bed_tags:
        for i in Bed_tags:
            sp = i.rstrip().split("\t", 3)
            if len(sp) > 3:
                key_pos = (int(sp[1])//window) * window
                site = key_pos + window
                start = int(sp[1])
                end = int(sp[2])
                if 2 * site >=  start + end and key_pos in tag_c[sp[0]]:
                    tag_c[sp[0]][key_pos] += 1
                elif site in tag_c[sp[0]] and key_pos+window in tag_c[sp[0]]:
                    tag_c[sp[0]][key_pos + window] += 1

    print ("Tag number analysis starts:")

    tag_OUT = open("tag_density.txt", 'w')

    for i in CHR_N:
        for j in sorted(tag_c[i]):
        #tag_norm = math.log(tag_c[i][j]/1000 + 1, 10)
            tag_norm = (tag_c[i][j]/1000 + 1)/ 10**1.5
            tag_OUT.write("{0} {1} {2} {3}\n".format(CHR_N[i], j, j+window-1, math.log(tag_norm +1, 2)))

    tag_OUT.close()

    SNP_OUT = open("snp_density.txt", 'w')
    for i in CHR_N:
        CL = int(chr_info[i])   # CL is chromosome Length
        if CL >= window:
            for j in range(0, CL, window):
                if j in chromo[i]:
                    end = j + window
                    snp_norm = math.log(chromo[i][j] + 1, 10)
                    SNP_OUT.write("{0}\t{1}\t{2}\t{3}\n".format(CHR_N[i], j, end-1 , snp_norm))

    SNP_OUT.close()

def LenStat_sample(uc_dir):
    if not os.path.isdir(uc_dir):
        sys.exit("Directory does not exist. Please check it.")
    results_UC = {}
    sample = []
    for root, dirs, files in os.walk(uc_dir):
        for fr in files:
            if fr.endswith("uc"):
                SI = fr.rstrip("_tmp.uc").lstrip("results_sim_")
                results_UC[SI] = os.path.join(root, fr)
            elif fr.startswith("sample") and fr.endswith("fa"):
                sample.append(fr)
    Sample_N = list(results_UC.keys())
    for i in sample:
        sample_ID = i.lstrip("sample_").rstrip(".fa")
        if sample_ID not in Sample_N:
            print("uclust analysis of {0} may not been finished.".format(i))

    return results_UC

def get_length_stat(f, stat_f):
    the_Stat_File = open( f , "r")
    sample_name = f.lstrip("Unique_").rstrip("_ucstat.txt")
    len_Of_Seq = []
    Tot_Len = 0

    for record in the_Stat_File:
        rsp = record.rstrip().split("\t")
        length = int(rsp[1])
        len_Of_Seq.append(length)
        Tot_Len += length
    arr = np.array(len_Of_Seq)

    mx_len = np.max(arr)
    mean_len = np.mean(arr)
    General_stat = open(stat_f, "a")
    General_stat.write("The Total number of Tags in sample {0} is: {1}\n".format(sample_name, arr.size))
    General_stat.write("The Min length of Tags in sample {0} is: {1}\n".format(sample_name, np.min(arr)))
    General_stat.write("The Max length of Tags in sample {0} is: {1}\n".format(sample_name, mx_len))
    General_stat.write("The Mean length of Tags in sample {0} is: {1}\n".format(sample_name, mean_len))
    General_stat.close()

    hist, bin_edges = np.histogram(arr, bins = np.arange(0, mx_len + 20, 20), density = False)

    Matri_Len = {}

    for n,i in enumerate(bin_edges):
        try:
            Matri_Len[str(i) +"-" + str(i+20)] = hist[n]
            #print ("{0}-{1}\t{2}".format(i,i+20, hist[n]))
        except IndexError:
            pass
            #print ("{0}-{1}\texceed Number range".format(i, i+20))
    stat_result = open(sample_name + "_SD.txt", 'w')
    stat_result.write("The distribution of length is:\n")
    Sort_Mat = sorted(Matri_Len.iteritems(), key = lambda x: int(x[0].split("-")[0]))
    for n,v in Sort_Mat:
        stat_result.write("{0}\t{1}\n".format(n, v))
    stat_result.close()

def Tag_Num(summary_file):
    LEN_dis = defaultdict(list)
    try:
        File_Handle = open(summary_file, "r" )
        for line in File_Handle:
            lc = line.rstrip()
            if line.startswith("The Total number of Tags "):
                start_i = lc.index("in sample ")
                end_i = lc.index(" is:")
                S_Name_Tag = lc[start_i+10:end_i]
                LEN_dis[S_Name_Tag].append("Tot\t"+lc.split()[-1])
            elif line.startswith("The Mean length of Tags "):
                start_i = lc.index("in sample ")
                end_i = lc.index(" is:")
                S_Name_Tag = lc[start_i+10:end_i]
                LEN_dis[S_Name_Tag].append("Ave\t"+lc.split()[-1])
    except IOError:
        print ("Summary file not Found.")
    else:
        File_Handle.close()
    return LEN_dis
'''
if __name__ == "__main__":

    Args = load_config("options.yaml")

    R_handle = open("functions_file.R", 'w')

    R_handle.write("#!/usr/bin/env Rscript\n\n")
    R_handle.write("library(RColorBrewer)\nlibrary(ggplot2)\nlibrary(reshape2)\n\nlibrary(gdsfmt)\nlibrary(SNPRelate)\n\n")
    R_handle.write("trim.suf <- function(x)sub(\"\\..*$\",\"\",x)\n\n")

    for j in Args["R_functions"]:
        R_handle.write(j)
        R_handle.write(Args["R_functions"][j].lstrip() + "\n\n")
    R_handle.close()

    Script_h = open ("Invoke_R_Fun.R", 'w')
    Script_h.write('#!/usr/bin/env Rscript\n\nsource("functions_file.R")\n\n')
    Script_h.write(Args["Plot_script"]+"\n")
    Script_h.close()


    USAGE = "Define input file and Total Sequencing Depth of each variant site"

    parser = argparse.ArgumentParser(description = USAGE)
    parser.add_argument("-i", "--input", action = "store", required = True, help = "the vcf file containing variants, including SNP and INDEL")
    parser.add_argument("-D", "--TD", type = int, default = 35, help = "the Depth in INFO field of Variants")
    args = parser.parse_args()

    # Vcf_Input = Args["receive_VCF"]
    # TD = Args["filter_Op"]["FO_TD"]
    Vcf_Input = args.input
    print Vcf_Input

    TD = args.TD
    print TD
    
    Qst = Args["filter_Op"]["FO_Qst"]
    print Qst

    Pre = Args["filter_Op"]["output_Prefix"]
    print Pre
    SQ_Suffix = Args["filter_Op"]["FO_sf"]
    print SQ_Suffix
    MaxDepth = Args["filter_Op"]["FO_MaxDepth"]
    print MaxDepth
    Interv = Args["filter_Op"]["FO_interv"]
    print Interv
    Cov =Args["filter_Op"]["FO_Cov"]
    print Cov
    
    Filtered_SNP = Pre + "_output_SNP.txt"
    print Filtered_SNP    

    # check if the vcffile exists
    if not os.path.exists(Vcf_Input):
        sys.exit("The vcf file name is required, but it is not found.\nPlease check it and provide correct file path.")
                
    Tot_VN, SNP_N, PMV_Tot_C = filter_SNP(Vcf_Input, TD, Qst , Pre)
    print ("The number of all the variants filtered through func is {0}.".format(Tot_VN))
    print ("The number of all the SNPs filtered through func is {0}.".format(SNP_N))
    #print ("The number of all the INDELs filtered through func is {0}.".format(IND_N))

    for k, v in sorted(PMV_Tot_C.items(), key = lambda x: x[0]):
        print ("{0}\t{1}".format(k,v))

    QS_out = open("QUAL_Stats_"+SQ_Suffix +".txt", 'w')
    QS_out.write("Sample\tQual\n")
    # for i in QOS:
    #    for j in QOS[i]:
    #        QS_out.write("{0}\t{1}\n".format(i, j))
    QS_out.close()


    Dep_Gradient, Samples_Depth = coverage_stats(Vcf_Input, MaxDepth, Interv)
    print ("The samples in this vcf are: ")
    for each_sample in Samples_Depth:
        print ("The depth of {0}:".format(each_sample))
        Coverage_Dict = Samples_Depth[each_sample]
        for j in sorted(Coverage_Dict.items(), key = lambda x: x[0]):
            print ("{0}\t{1}".format(j[0], j[1]))
    
    output_file = open(Cov, 'w')

    sample_list = []
    output_file.write("Minimal_Cov")
    for i in Samples_Depth:
        output_file.write("\t{0}".format(i))
        sample_list.append(i)
    output_file.write("\n")

    Depths = [x[0] for x in sorted(Dep_Gradient.items(), key = lambda x: x[0])]
    for Gra in Depths:
        output_file.write("{0}".format(Gra))
        for s in sample_list:
            output_file.write("\t{0}".format(Samples_Depth[s][Gra]))
        output_file.write("\n")
    output_file.close()
    print ("The sample number in vcf is {0}.\n Analysis finished.".format(len(Samples_Depth)))


############################################################################################
### the output files is ready for Plot_e.R ###

    
    PMV_SNP = Pre + "_PMV_SNP.txt"

    stat_1, stat_2, stat_3 = data_ready(Filtered_SNP)
    ## write SNP statistic results according PMV and MAF
    PMV_out = open(PMV_SNP, 'w')
    PMV_out.write("MAF_value")

    sort_v = sorted(list(stat_1.keys()))
    for v in sort_v:
        PMV_out.write("\tP{0}".format(v))
    PMV_out.write("\n")

    for i in sorted(stat_1[1]):
        PMV_out.write("{0}".format(i * 0.01))
        for j in sort_v:
            PMV_out.write("\t{0}".format( stat_1[j][i] ))
        PMV_out.write("\n")

    PMV_out.close()

    ## write SNP number of different maf thresholds
    handle_maf = open("maf_stats.txt", 'w')
    handle_maf.write("MAF\tSNP_NUM\n")
    for i in sorted(stat_2.items(), key = lambda x: x[0]):
        handle_maf.write("{0}\t{1}\n".format(i[0]*0.01, i[1]))
    handle_maf.close()

    ## write SNP number of different Herterozygosity
    handle_het = open("het_stats_v.txt", 'w')
    handle_het.write("HET\tSNP_NUM\n")

    for i in sorted(list(stat_3.keys())):
        if i % 2 == 0 and i < 98:      
            handle_het.write("{0}\t{1}\n".format( i * 0.01, stat_3[i] + stat_3[i+1]))
        elif i == 98:
            handle_het.write("{0}\t{1}\n".format( i * 0.01, stat_3[i] + stat_3[i+1] + stat_3[i+2]))
        elif i == 100:
            handle_het.write("{0}\t{1}\n".format( i * 0.01, 0))
    handle_het.close()


###########################################################################################################

### the output file will be ready for plotting with Plot_SNP.R ###

    
    Win_S = Args["Plot_Chr_Op"]["PCO_WS"]
    Arg_Step = Args["Plot_Chr_Op"]["PCO_Step"]
    Chr_name = Args["Plot_Chr_Op"]["List_F"]
    bamfile = Args["Plot_Chr_Op"]["bam_fY"]
    Res_endo_bed = Args["Plot_Chr_Op"]["PCO_Res_endo"]
    SD_output = Pre + "_SNPnum_distribution.txt"
    
    if not os.path.exists(Chr_name):
        print ("""the list file for chromosome does not exist. please provide this file. It should seem like this:  
          1       NC_029679.1
          2       NC_029680.1
          3       NC_029681.1
          4       NC_029682.1
          5       NC_029683.1
          6       NC_029684.1
""")
        sys.exit("Process aborted abnormally.")
    
    # invoke two functions to process data
    data_process(Filtered_SNP, Win_S, Arg_Step, bamfile, SD_output)
    SC, CHO = chrome_data_prepare(Filtered_SNP, Win_S, Arg_Step, bamfile, Chr_name)

    # write the data to files
    SC_rev = sorted(SC, reverse=True)
    sim_bed = open("simple_BED.txt" ,  "w")  # this is a simple bed file , or simulated bed for chromosomes
    sim_bed.write("Chr\tLength\n")
    for i in SC_rev:
        print ("The name of chromosome in this specie is: {0}\t{1}".format(i[1],i[0]))
        sim_bed.write("{0}\n".format(i[1]))
    sim_bed.close()

    sim_region = open("simple_region.txt", 'w')
    sim_region.write("Chr\tstart\tend\tvalue\tcolor\n")

    for i in SC:
        chr_length = int(i[1].split("\t")[1])
        chr_n = i[1].split("\t")[0]
        if chr_length >= Win_S:
            for j in range(0, chr_length, Arg_Step):
                value = CHO[i[0]].get(j, 0)
                if j + Win_S <= chr_length + Arg_Step:
                    locus = j + Win_S/2
                    sim_region.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(chr_n, j, j+Win_S-1, value, 0))   
                #if j + Win_S <= chr_length + Arg_Step and value > 80:
                #    locus = j + Win_S/2
                #    sim_region.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(chr_n, j, j+Win_S-1, value * 200, 0))
    sim_region.close()
    
    KI = [i[1] for i in SC]
    KC = [i[0] for i in SC]

    if os.path.exists(Res_endo_bed):
        o_b = open("Res_bed.txt", 'w')
        d_c = defaultdict(list)
        b_f = open(Res_endo_bed, 'r')
        for i in b_f:
            fi = i.rstrip().split("\t",3)
            if fi[0] in KC:
                ci = KC.index(fi[0])
                ch_name_Bed = KI[ci].split("\t")[0]
                o_b.write(ch_name_Bed + "\t" + "\t".join(fi[1:]) + "\n")
                r_w = int(((int(fi[1]) + int(fi[2]))/2)/Win_S)*Win_S
                d_c[ch_name_Bed+"\t"+str(r_w)].append(int(fi[2]) - int(fi[1]))
        o_b.close()
        b_f.close()
        Res_bed_s = open("Size_bed_win.txt", 'w')
        Res_bed_s.write("Chr\tstart\tend\tvalue\tcolor\n")
        for i in d_c:
            ir = int(i.split("\t")[1])+Win_S
            Res_bed_s.write("{0}\t{1}\t{2}\t0\n".format(i, ir, sum(d_c[i])))
        Res_bed_s.close()
 
    else:
        print("No digest bedfile found. Please check it.")


############################################################################################
    
    Circ_win = Args["Circos_Op"]["CO_win"]
    List_C = Args["Circos_Op"]["List_F"]
    C_bamfile = Args["Circos_Op"]["bam_fY"]
    Abbr = Args["Circos_Op"]["CO_Abbr"]
    bedfile = Args["bed_fY"]
    
    if (not os.path.exists(C_bamfile)) or (not os.path.exists(bedfile)):
        print ("""bam and bed file required, please the path of these two file, 
    if you have bam file but not bed file, use the following command to generated it.
         /usr/bin/bamToBed -i some.bam > some.bed""")
        sys.exit("Quit for some file missing.")
    prepare_DataFor_Circ(Filtered_SNP, Circ_win, List_C, Abbr, C_bamfile, bedfile )
    #CHR_Tran = "/home/yangjy/16T/GENOME/GEN160781/GEN160781/Ref_Genome/List.txt"

################################################################################################
    
    Tag_OutDir = Args["SLDO"]["SO_outDir"]     # A directory containg fasta files
    Tag_Stat_sum = Args["SLDO"]["SO_Summ"]
    final_summary = Args["SLDO"]["SO_final_sum"]

    ucs = LenStat_sample(Tag_OutDir)
    for i in list(ucs.keys()):
        sn = "Unique_"+ i + "_ucstat.txt"
        if not os.path.exists(sn):
            cmd = "awk -F\'\\t\' -v OFS=\'\\t\' \'/^S/{print $1,$3}\' " + ucs[i] + " > " + sn
            print (cmd)
            Unique_UC = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            Unique_UC.wait()

        #if Unique_UC.poll() == 0 and os.path.getsize(sn) !=0:
        if os.path.getsize(sn) != 0:
            get_length_stat(sn, Tag_Stat_sum)
        else:
            sys.exit("Sample_{0} Tag information cannot be retrieved.".format(sn.lstrip("Unique_").rstrip("_ucstat.txt")))
    
    D = Tag_Num(Tag_Stat_sum)

    if os.path.exists(final_summary):
        sys.exit("the summary file exists, please check and use another name.")

    STs = open(final_summary, 'w')   # STs is abbreviation of Sample_Tag_stat

    STs.write("Sample_Name\tTotal_Num\tAverage_Length\n")
    for k,i in D.items():
        STs.write("%s\t"%(k))
        num_L = ["", ""]
        for j in i:
            if j.startswith("Tot"):
                num_L[0] = str(round(float(j.split("\t")[1])))
            elif j.startswith("Ave"):
                num_L[1] = str(round(float(j.split("\t")[1])))
        STs.write("\t".join(num_L)+"\n")
    STs.close()

    for i in glob("Unique_*ucstat.txt"):
        os.remove(i)



