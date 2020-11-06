top_level_path = fullfile('..','..');
swat_gis_path = fullfile(top_level_path,'SWAT_Model','gis_data');
%% Read Shapefiles
SAC_cvhm = shaperead(fullfile(swat_gis_path,'CVHM_SAC_3310'));
SJV_cvhm = shaperead(fullfile(swat_gis_path,'CVHM_SJV_3310'));
SSJV_cvhm = shaperead(fullfile(swat_gis_path,'CVHM_SSJV_3310'));
%% Unique Codes in each area
SAC_HRU_codes = unique([SAC_cvhm.HRU_GIS1]');
SJV_HRU_codes = unique([SJV_cvhm.HRU_GIS1]');
SSJV_HRU_codes = unique([SSJV_cvhm.HRU_GIS1]');
save('CVHM_HRU_GIS_CODES_UNIQUE','SAC_HRU_codes','SJV_HRU_codes','SSJV_HRU_codes');
%% Or load the unique codes
load('CVHM_HRU_GIS_CODES_UNIQUE');
%% Load N and Rch data
scen = 4;
SAC_TAB = readtable(fullfile(swat_gis_path,'..','NLoading','SWAT',['SAC_scen' num2str(scen) '.dat']));
SJV_TAB = readtable(fullfile(swat_gis_path,'..','NLoading','SWAT',['SJV_scen' num2str(scen) '.dat']));
SSJV_TAB = readtable(fullfile(swat_gis_path,'..','NLoading','SWAT',['SSJV_scen' num2str(scen) '.dat']));
%% Make a unique list of SWAT LAND USES
load(fullfile(swat_gis_path, '..','SWAT_LULC_NAMES'))
for ii = 1:length(SWATLULCNAMES)
    CODE{ii,1} = strip(SWATLULCNAMES{ii,1}{1});
    NAME{ii,1} = strip(SWATLULCNAMES{ii,2}{1});
end
SWAT_LULC = table(CODE, NAME);
SWAT_LULC.CODE = categorical(SWAT_LULC.CODE);
SWAT_LULC = sortrows(SWAT_LULC);
SWAT_LULC.ID = [1:57]';
% writetable(SWAT_LULC,fullfile('..','MantisData','SWAT_LULC.csv'));
%% Keep the yearly values from the dataset
SAC_TAB_yr = SAC_TAB(SAC_TAB.MON>1000,:);
SJV_TAB_yr = SJV_TAB(SJV_TAB.MON>1000,:);
SSJV_TAB_yr = SSJV_TAB(SSJV_TAB.MON>1000,:);
%% Reshape the data according to years 
clear SAC SJV SSJV
% SAC
SAC.hru = reshape(SAC_TAB_yr.HRUGIS, size(SAC_TAB_yr,1)/25, 25);
SAC.hru = SAC.hru(:,1);
SAC.lulc = reshape(SAC_TAB_yr.LULC, size(SAC_TAB_yr,1)/25, 25);
SAC.lulc = categorical(SAC.lulc(:,1));
[Lia, Locb] = ismember(SAC.lulc,SWAT_LULC.CODE);
SAC.luid = Locb;
SAC.NO3 = reshape(SAC_TAB_yr.NO3Lkg_ha, size(SAC_TAB_yr,1)/25, 25);
SAC.Rch = reshape(SAC_TAB_yr.GW_RCHGmm, size(SAC_TAB_yr,1)/25, 25);
% SJV
SJV.hru = reshape(SJV_TAB_yr.HRUGIS, size(SJV_TAB_yr,1)/25, 25);
SJV.hru = SJV.hru(:,1);
SJV.lulc = reshape(SJV_TAB_yr.LULC, size(SJV_TAB_yr,1)/25, 25);
SJV.lulc = categorical(SJV.lulc(:,1));
[Lia, Locb] = ismember(SJV.lulc,SWAT_LULC.CODE);
SJV.luid = Locb;
SJV.NO3 = reshape(SJV_TAB_yr.NO3Lkg_ha, size(SJV_TAB_yr,1)/25, 25);
SJV.Rch = reshape(SJV_TAB_yr.GW_RCHGmm, size(SJV_TAB_yr,1)/25, 25);
% SSJV
SSJV.hru = reshape(SSJV_TAB_yr.HRUGIS, size(SSJV_TAB_yr,1)/25, 25);
SSJV.hru = SSJV.hru(:,1);
SSJV.lulc = reshape(SSJV_TAB_yr.LULC, size(SSJV_TAB_yr,1)/25, 25);
SSJV.lulc = categorical(SSJV.lulc(:,1));
[Lia, Locb] = ismember(SSJV.lulc,SWAT_LULC.CODE);
SSJV.luid = Locb;
SSJV.NO3 = reshape(SSJV_TAB_yr.NO3Lkg_ha, size(SSJV_TAB_yr,1)/25, 25);
SSJV.Rch = reshape(SSJV_TAB_yr.GW_RCHGmm, size(SSJV_TAB_yr,1)/25, 25);
%% Convert to concentration
% SAC
% The units of NO3 are Kg/ha/year
NO3_tmp = SAC.NO3 * 1000000; %mg/ha/year
NO3_tmp = NO3_tmp / 10000; %mg/m^2/year
NO3_tmp = NO3_tmp / 1000000; %mg/mm^2/year
% The units of recharge are mm/year
SAC.Con = NO3_tmp ./ SAC.Rch;% mg/mm^3
SAC.Con = SAC.Con * 1000000;% mg/L
SAC.Con(isnan(SAC.Con)) = 0;
% SJV
% The units of NO3 are Kg/ha/year
NO3_tmp = SJV.NO3 * 1000000; %mg/ha/year
NO3_tmp = NO3_tmp / 10000; %mg/m^2/year
NO3_tmp = NO3_tmp / 1000000; %mg/mm^2/year
% The units of recharge are mm/year
SJV.Con = NO3_tmp ./ SJV.Rch;% mg/mm^3
SJV.Con = SJV.Con * 1000000;% mg/L
SJV.Con(isnan(SJV.Con)) = 0;
% SSJV
% The units of NO3 are Kg/ha/year
NO3_tmp = SSJV.NO3 * 1000000; %mg/ha/year
NO3_tmp = NO3_tmp / 10000; %mg/m^2/year
NO3_tmp = NO3_tmp / 1000000; %mg/mm^2/year
% The units of recharge are mm/year
SSJV.Con = NO3_tmp ./ SSJV.Rch;% mg/mm^3
SSJV.Con = SSJV.Con * 1000000;% mg/L
SSJV.Con(isnan(SSJV.Con)) = 0;
%% Isolate the HRUs within the Central Valley
% SAC
[C, ia, ib] = intersect(SAC_HRU_codes, SAC.hru);
SAC.CVhru = C;
SAC.CVcon = SAC.Con(ib,:);
SAC.CVluid = SAC.luid(ib,1);
% SJV
[C, ia, ib] = intersect(SJV_HRU_codes, SJV.hru);
SJV.CVhru = C;
SJV.CVcon = SJV.Con(ib,:);
SJV.CVluid = SJV.luid(ib,1);
% SSJV
[C, ia, ib] = intersect(SSJV_HRU_codes, SSJV.hru);
SSJV.CVhru = C;
SSJV.CVcon = SSJV.Con(ib,:);
SSJV.CVluid = SSJV.luid(ib,1);
%% Concatenate the three regions
pre = 10000000;
SWAT_conc = [pre + SAC.CVhru SAC.CVluid SAC.CVcon; ...
    pre + SJV.CVhru SJV.CVluid SJV.CVcon; ...
    pre + SSJV.CVhru SSJV.CVluid SSJV.CVcon];
SWAT_conc = [[1:size(SWAT_conc,1)]' SWAT_conc];
%% Write loading
frmt = '%d %d';
for ii = 1:25
    frmt = [frmt ' %.5f'];
end
frmt = [frmt '\n'];
fid = fopen(fullfile('..','MantisData',['SWAT_LOADING_SCEN_' num2str(scen) '.dat']),'w');
fprintf(fid, frmt,SWAT_conc(:,2:end)');
fclose(fid);
%% Identify the swat_index 
