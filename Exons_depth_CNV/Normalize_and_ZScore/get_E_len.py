

## get bed region Length

## exon_file

r_file = open("exons_length.txt", 'w')
a = open("66genes_exons.txt", 'r')

for line in a:
    x = line.rstrip().split('\t')

    length = int(x[2]) + 1 - int(x[1])

    r_file.write("{0}\t{1}\n".format(x[3],length))

r_file.close()
a.close()

