import os
import sys
import pandas as pd
from scipy.io import loadmat
import numpy as np
import random
import matplotlib.pyplot as plt


# Retrieve current working directory (`cwd`)
cwd = os.getcwd()
print (cwd)
print(sys.version)

LU_cat_file = '../Local/LanduseTable_2017_0515.xlsx'
LU_cat_xl = pd.ExcelFile(LU_cat_file)

# https://pandas.pydata.org/pandas-docs/version/0.23/generated/pandas.read_excel.html
df = LU_cat_xl.parse(LU_cat_xl.sheet_names[0],0,None,None,'A:B')

print(LU_cat_xl.sheet_names[0])
print(df.iloc[15,1])
print(df.iloc[15,0])
print(df.size)

# Assign dummy loading reductions
Perc_reduc = []
for i in range(0,df.size,1):
    Perc_reduc.append(random.uniform(50,100))
    #print(Perc_reduc[i])

# # find column names in the dataframe
# print(df.columns)
# # search in a column for a particular value/code
# ind = df.loc[df['DWR/CAML Code'] == 20008]
# print('# Rows ' + str(len(ind)))
# # get the row index of the first row of the search.
# #print(ind.index[0])
#
# #print(Perc_reduc[ind.index[0]])

print('Loading Data...')
Data = loadmat('../Local/data4python.mat')
Data1 = loadmat('../Local/URFdata.mat')

LUmaps = {}
for i in range(1945, 2006, 15):
    LUmaps[i] = Data['LU' + str(i)]

Ngw = {}
for i in range(1945,2051,15):
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

LF_base = np.empty((len(urfs), 105))
LF_red = np.empty((len(urfs), 105))

# Build Loading Functions
print('Building Loadinf Functions...')
for i in range(0,1000,1):  #len(urfs)
    #print(i)
    # For each loading function we want to build
    # We will loop through the time snapshots of land use and Ngw ( every 15 years)
    # We will extract a value at the start of the period and a value at the end of the period
    # We will also extract the land use id at the start of the period and at the end

    # loading function counter
    lfi = 0

    for k in range(1945,2050,15): # This will loop through the 2035
        #print('--------' + str(k) + '------------')
        I = Sij[i][0]
        J = Sij[i][1]
        val_start = Ngw[k][I][J]
        val_end = Ngw[k+15][I][J]

        #print('V start ' + str(val_start) + ' ' + ' V end ' + str(val_end))


        # klu is the index for the land use. After 2005 we use the 2005 land use map
        klu_s = k
        klu_e = k + 15
        if k > 1990:
            klu_s = 2005
            klu_e = 2005

        lu_s = LUmaps[klu_s][I][J]
        lu_e = LUmaps[klu_e][I][J]
        #print('klu start ' + str(klu_s) + ' ' + ' klu end ' + str(klu_e))
        #print('LU start ' + str(lu_s) + ' ' + ' LU end ' + str(lu_e))

        # find out the user defined reduction for the land use categories at the start and end of the current
        # 15 year period (lu_s and lu_e)
        ind = df.loc[df['DWR/CAML Code'] == lu_s]
        # initialize reduction values to no reduction
        red_s = 1
        red_e = 1
        if len(ind) > 0:
            red_s = Perc_reduc[ind.index[0]]/100

        ind = df.loc[df['DWR/CAML Code'] == lu_e]
        if len(ind) > 0:
            red_e = Perc_reduc[ind.index[0]]/100

        #print('red start ' + str(red_s) + ' ' + ' red end ' + str(red_e))

        # distribute the loading function
        lf_base_temp = np.linspace(val_start, val_end, num=15)


        lf_red_temp = np.linspace(val_start*red_s, val_end*red_e, num=15)
        #print(lf_red_temp)
        for t in range(0,15,1):
            LF_base[i,lfi] = lf_base_temp[t]
            LF_red[i,lfi] = lf_red_temp[t]
            lfi += 1
        #print(lfi)

    #print(LF_base[i,:])


# Plot and convolute URF with LF

for i in range(0,10,1):
    r = random.randint(0,1000) #len(LF_base)
    print(r)
    plt.plot(range(1945,2050,1),LF_base[r,:], 'b', range(1945,2050,1),LF_red[r,:], 'r')
    plt.legend(('Base LF', 'Reduced LF'), loc='upper left', shadow=True)
    plt.ylabel('Loading function')
    plt.ylabel('Time')
    plt.show()

    # convolute URFs with LFs
    btc_base = np.zeros(105)
    btc_red = np.zeros(105)
    shift = 0
    for ii in range(0,105,1):
        for jj in range(shift,105,1):
            btc_base[jj] = btc_base[jj] + urfs[r, jj-shift] * LF_base[r, ii]
            btc_red[jj] = btc_red[jj] + urfs[r, jj - shift] * LF_red[r, ii]
        shift += 1

    plt.plot(range(1945, 2050, 1), btc_base, 'b', range(1945, 2050, 1), btc_red, 'r')
    plt.show()



# How to index this dictionary
print(Ngw[1945][6200,2876])










