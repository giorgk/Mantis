import os
import sys
import pandas as pd
#import matplotlib.pyplot as plt
from PIL import Image

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

#LU_1945 = plt.imread('../Local/model_input_LU1945.tif')
im = Image.open('../Local/Ngw_1945.tif')
#im.show()
#print(LU_1945)
