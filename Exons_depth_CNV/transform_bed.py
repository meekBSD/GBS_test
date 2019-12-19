import sys

transID_file = sys.argv[1]
transID_dict = {}

a = open(transID_file, "r")
for line in a:
    sp = line.rstrip().split('\t')
    x = line.rstrip().split("\t")[0].split('.')
    transID_dict[x[0]] = sp[1]

a.close()

b = open("Test_genes_exons.xls", 'r')
for line in b:
    k = line.rstrip().split('\t')
    print('{0}\t{1}\t{2}\t{3}\t0\t{4}'.format(k[2], k[3], k[4], transID_dict[k[0]]+'_'+k[1]+'.' + k[0], k[5] ))
b.close()
