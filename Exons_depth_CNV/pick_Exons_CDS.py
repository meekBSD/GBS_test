

import sys

transID_file = sys.argv[1]
refGene_file = sys.argv[2]


def get_RegionInfo(transcript_ID, chromosome, strand, exStart, exStop, cStart, cStop):
    exonNum = len(exStart)
    exon_total   = 0
    coding_total = 0
    if strand == "-":
        cStart, cStop = cStop, cStart
        ExonList = []
        CDSList  = []

        for j in range(exonNum-1, -1 , -1):
            
            exon_Len = abs(int(exStart[exonNum-1-j]) - int(exStop[exonNum-1-j]))
            if j == exonNum -1:
                CDS_Len = abs(int(cStart) - int(exStop[0]))
            elif j == 0:
                CDS_Len = abs(int(exStart[exonNum -1]) - int(cStop))
            else:
                CDS_Len = exon_Len
            
            ExonList.append(exon_Len)
            CDSList.append(CDS_Len)
        for j in range(exonNum-1, -1 , -1):    
            n = exonNum-1-j
            exon_total   = sum(ExonList[n:])
            coding_total = sum(CDSList[n:])
            print("{0}\texon{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}".format( transcript_ID,j+1, x[2] ,exStart[n], exStop[n], strand, ExonList[n], CDSList[n], exon_total, coding_total ))
    else:
        #for n, s in enumerate(xsta):
        
        for j in range(0, exonNum):
            exon_Len = abs(int(exStart[j]) - int(exStop[j]))
            if j == 0:
                CDS_Len = abs(int(cStart) - int(exStop[0]))
            elif j == exonNum -1:
                CDS_Len = abs(int(exStart[j]) - int(cStop))
            else:
                CDS_Len = exon_Len
                
            exon_total   += exon_Len
            coding_total += CDS_Len
            print("{0}\texon{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}".format( transcript_ID ,j+1, x[2],exStart[j], exStop[j], strand, exon_Len, CDS_Len, exon_total, coding_total ))



transID_list = []

a = open(transID_file, "r")
for line in a:
    x = line.rstrip().split("\t")[0].split('.')
    transID_list.append(x[0])

a.close()

refHandle = open(refGene_file, 'r')
for line in refHandle:
    x = line.rstrip().split('\t')
    if x[1] in transID_list:
        strand = x[3]
        
        cd_s = x[6]
        cd_e = x[7]
        estar  = x[9]
        estop  = x[10]

        xsta = estar.split(",")
        xsta.remove("")
        xsto = estop.split(",")
        xsto.remove("")
        get_RegionInfo(x[1], x[2], strand, xsta, xsto, cd_s, cd_e)
        #print(line.rstrip())


refHandle.close()


