import os
import sys
import pandas as pd
from scipy.io import loadmat
import numpy as np


# Retrieve current working directory (`cwd`)
cwd = os.getcwd()
print (cwd)
print(sys.version)

LU_cat_file = '../Local/LanduseTable_2017_0515.xlsx'
LU_cat_xl = pd.ExcelFile(LU_cat_file)

# https://pandas.pydata.org/pandas-docs/version/0.23/generated/pandas.read_excel.html
df = LU_cat_xl.parse(LU_cat_xl.sheet_names[0],1,None,None,'A:B')

print(LU_cat_xl.sheet_names[0])
print(df.iloc[10,1])

print('Loading Data...')
Data = loadmat('../Local/data4python.mat')
Data1 = loadmat('../Local/URFdata.mat')

LUmaps = {}
for i in range(1945, 2005, 15):
    LUmaps[i] = Data['LU' + str(i)]

Ngw = {}
for i in range(1945,2050,15):
    Ngw[i] = Data['Ngw' + str(i)]

LUcat = Data['LUcat']
print('# Categories: ' + str(len(LUcat)))

Sxyv = Data1['Sxyv']
Sid = Data1['Sid']
Sij = Data1['SIJ']
urfs = Data1['urfs']
urfV = Data1['urfV']

print(len(urfs))

del Data
del Data1

LF = np.empty((len(urfs), 105))

# Build Loading Functions
print('Building Loadinf Functions...')
for i in range(0,len(urfs),1):
    for k in range(1945,2049,15):
        I = Sij[i][0]
        J = Sij[i][1]
        val_start = Ngw[k][I][J]
        val_end = Ngw[k+15][I][J]

        klu_s = k
        klu_e = k + 15
        if k > 1990:
            klu_s = 2005
            klu_e = 2005









print(Ngw[1945][6200,2876])

