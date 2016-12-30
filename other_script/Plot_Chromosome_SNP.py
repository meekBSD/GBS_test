#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from collections import defaultdict
import os, sys
import argparse

### the output file will be ready for plotting with Plot_SNP.R ###

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

if __name__ == "__main__":

    USAGE = "python %s.py -i SNP_filter_result"%("Plot_Ch.py")

    parser = argparse.ArgumentParser(description = USAGE)
    parser.add_argument("-i", "--input", action = "store", required = True, help = "the SNP tab delimited file containing variants, pos and chromosome")
    parser.add_argument("-w", "--WS", type = int, default = 100000, help = "the window size of region")
    parser.add_argument("-s", "--Step", type = int , default = 10000,  help = "the step ")
    parser.add_argument("-rd", "--Res_endo_bed",  help = "the bedfile of Restriction endonuclease digestion ")
    parser.add_argument("-b", "--bamfile", action = "store", required = True, help = "the bam_file containing chromosome length")
    parser.add_argument("-c", "--Chr_name", action = "store", default = "List.txt", help = "the filename containing alias name of chromosome for NCBI deposited ID")
    parser.add_argument("-o", "--output", action = "store", default = "SNP_num_distribution.txt", help = "the Output file")

    args = parser.parse_args()
    Win_S = args.WS
    Arg_Step = args.Step
    
    if not os.path.exists(args.Chr_name):
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
    data_process(args.input, args.WS, args.Step, args.bamfile, args.output)
    SC, CHO = chrome_data_prepare(args.input, Win_S, Arg_Step, args.bamfile, args.Chr_name)

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

    if os.path.exists(args.Res_endo_bed):
        o_b = open("Res_bed.txt", 'w')
        d_c = defaultdict(list)
        b_f = open(args.Res_endo_bed, 'r')
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


