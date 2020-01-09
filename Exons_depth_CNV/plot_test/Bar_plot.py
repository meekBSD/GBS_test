#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import matplotlib
import numpy as np
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import Counter

e_start = 12022903

## store Y values

Y = []

n = 0
h1 = open("a1.txt", "r")
for line in h1:
    s = line.rstrip().split("\t")
    if len(s) == 2 :
        Y.append(e_start + int(s[0]))

h1.close()
   
c = Counter(Y)
Top_c = c.most_common(15)

xa = list(np.arange(15))

xt_labels = []
data = []

for n, f in Top_c:
    xt_labels.append(str(n))
    data.append(f)

plt.rcParams['figure.figsize'] = (9, 3)

fig=plt.figure(1)
ax1=plt.subplot(111)

rect=ax1.bar( x =xa, height=data,width=0.4,color="lightblue")
for rec in rect:
    x=rec.get_x()
    height=rec.get_height()
    ax1.text(x+0.1,1.02*height,str(height))

ax1.set_xticks( xa )
ax1.set_xticklabels(tuple(xt_labels), rotation=90)
ax1.set_ylabel("Frequency")
ax1.set_title("BreakPoint Frequency in Intron Region")
ax1.grid(True)
#ax1.set_ylim(0,28)

#plt.legend() 
#plt.show()
plt.savefig("frequency_1.png", bbox_inches="tight")
plt.close() 


