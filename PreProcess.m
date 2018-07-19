%% Make a unique list of land uses
LU = imread('Local/model_input_LU2005.tif');
LU_cat = unique(LU);
%% Find names for each LU category
% temp1 has all the numerical values and temp2 the text values
[temp1, temp2] = xlsread('Local/LanduseTable_2017_0515.xlsx', 'FINAL Landuse Table','A2:E208');
