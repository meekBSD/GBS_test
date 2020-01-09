#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

e_start = 12022903

## store Y values

X = []
Y = []

n = 0
h1 = open("a1.txt", "r")
for line in h1:
    s = line.rstrip().split("\t")
    if len(s) == 2 :
        Y.append(e_start + int(s[0]))

        X.append(n)
        n += 1

h1.close()
   
plt.rcParams['figure.figsize'] = (9, 3)

#x2=range(0,10) 
#y2=[5,8,0,30,20,40,50,10,40,15] 
plt.plot(X,Y,label='Frist line',linewidth=0.2,color='r',marker='o', 
markerfacecolor='blue',markersize=3.2) 
#plt.plot(x2,y2,label='second line') 
plt.xlabel('Plot Number') 
plt.ylabel('Genome Coordinate') 
plt.title('Points Graph\nCheck it out') 
plt.legend() 
#plt.show()
plt.savefig("out1.png", bbox_inches="tight")
plt.close() 


