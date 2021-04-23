#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import pandas as pd

a = pd.read_csv("low_cov_previous.xls", encoding="utf-8", sep="\t", header=None)
a.index = a.iloc[:,1]
b = pd.read_csv("TS00001.base.coverage.newcov.txt", encoding="utf-8", sep="\t", header=None)
b.index = b.iloc[:,0]

b= b.iloc[:,3]
result = pd.merge(a, b, left_index=True, right_index=True, how="inner")

c = pd.read_csv("TS00002.base.coverage.newcov.txt", encoding="utf-8", sep="\t", header=None)
c.index = c.iloc[:,0]
c = c.iloc[:,3]
result = pd.merge(result, c, left_index=True, right_index=True, how="inner")

d = pd.read_csv("TS00003.base.coverage.newcov.txt", encoding="utf-8", sep="\t", header=None)
d.index = d.iloc[:,0]
d = d.iloc[:,3]
result = pd.merge(result, d, left_index=True, right_index=True, how="inner")

result.to_csv('Result1.csv',index=0,  header=0, na_rep='NA', sep="\t")

