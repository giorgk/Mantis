import os
import sys
import pandas as pd
from scipy.io import loadmat
import numpy as np
import random


# Retrieve current working directory (`cwd`)
cwd = os.getcwd()
print (cwd)
print(sys.version)

LU_cat_file = '../Local/LanduseTable_2017_0515.xlsx'
LU_cat_xl = pd.ExcelFile(LU_cat_file)

# https://pandas.pydata.org/pandas-docs/version/0.23/generated/pandas.read_excel.html
df = LU_cat_xl.parse(LU_cat_xl.sheet_names[0],1,None,None,'A:B')

print(LU_cat_xl.sheet_names[0])
print(df.iloc[0,1])
print(df.iloc[0,0])
print(df.size)
# Assign dummy loading reductions
Perc_reduc = []
for i in range(0,df.size,1):
    Perc_reduc.append(random.uniform(50,100))
    print(Perc_reduc[i])


exit()

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
    # For each loading function we want to build
    # We will loop through the time snapshots of land use and Ngw ( every 15 years)
    # We will extract a value at the start of the period and a value at the end of the period
    # We will also extract the land use id at the start of the period and at the end
    for k in range(1945,2049,15):
        I = Sij[i][0]
        J = Sij[i][1]
        val_start = Ngw[k][I][J]
        val_end = Ngw[k+15][I][J]

        # klu is the index for the land use. After 2005 we use the 2005 land use map
        klu_s = k
        klu_e = k + 15
        if k > 1990:
            klu_s = 2005
            klu_e = 2005

        lu_s = LUmaps[klu_s][I][J]
        lu_e = LUmaps[klu_e][I][J]

        # find out the user defined reduction for the land use categories at the start and end of the current
        # 15 year period (lu_s and lu_e)











print(Ngw[1945][6200,2876])

