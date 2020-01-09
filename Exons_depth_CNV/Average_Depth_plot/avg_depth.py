import numpy as np
from glob import glob


png_files = glob("out*png")

Exonfiles = glob("Exons*xls")

outDep= open("all_exons_dep.xls", "w")
for i in png_files:
    ki = i[:-4].split("_")
    geneName = ki[1] + '_' + ki[2]

    depth_A = []
    for ef in Exonfiles:
        a = open(ef, 'r')
        depthL = []
        sample = ef.split('_')[1]

        for line in a:
            x = line.rstrip().split('\t')
            if x[-1] == geneName:
                
                depthL.append(int(x[2]))
        a.close()

        arrDep = np.array(depthL)
        print("{0}\t{1}\t{2}".format(geneName, sample, np.mean(arrDep) ))
        depth_A.append(np.mean(arrDep))
    
    geneDep = np.array(depth_A)
    outDep.write("{0}\t{1}\n".format(geneName, np.mean(geneDep)))

outDep.close()
