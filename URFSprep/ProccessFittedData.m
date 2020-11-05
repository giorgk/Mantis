top_level_path = fullfile('..','..');
%% Load data for CVHM BUD0
scen_name = 'CVHM_92_03_BUD0';
data_path = fullfile(top_level_path,'CVHM_NPSAT', 'Simulations',...
    'Sim_1992m10_2003m9','BudCoef_0','Hstd_406','output','final1');
load(fullfile(data_path,'FittedUrfs.mat'));
% remove the ones that do not exit from the top or side
id_rmv = find([wells.Ntop]' + [wells.Nside]' == 0);
wells(id_rmv,:) = [];
% Load the groundwater recharge function
cvhm_rch = read_Scattered(fullfile(top_level_path, 'CVHM_NPSAT','Simulations',...
    'Sim_1992m10_2003m9','BudCoef_0','cvhm_1992m10_2003m9_bf0_RCH.npsat'),2);
%% Load data for CVHM BUD1
scen_name = 'CVHM_92_03_BUD1';
data_path = fullfile(top_level_path,'CVHM_NPSAT', 'Simulations',...
    'Sim_1992m10_2003m9','BudCoef_1','Hstd_406','output');
load(fullfile(data_path,'FittedUrfs.mat'));
% remove the ones that do not exit from the top or side
id_rmv = find([wells.Ntop]' + [wells.Nside]' == 0);
wells(id_rmv,:) = [];
% Load the groundwater recharge function
cvhm_rch = read_Scattered(fullfile(top_level_path, 'CVHM_NPSAT','Simulations',...
    'Sim_1992m10_2003m9','BudCoef_1','cvhm_1992m10_2003m9_bf1_RCH.npsat'),2);
%% make a modified interpolation recharge with values all positive
irch_neg = find(cvhm_rch.v*365*1000 < 30);
irch_pos = cvhm_rch.v*365*1000 >= 30;
neg_xy = cvhm_rch.p(irch_neg,:);
pos_xy = cvhm_rch.p(irch_pos,:);
urf_rch_pos = cvhm_rch.v(irch_pos,:);
dst = zeros(size(neg_xy,1),1);
for ii = 1:size(neg_xy,1)
    [ccc, ddd] = min(sqrt((neg_xy(ii,1) - pos_xy(:,1)).^2 + (neg_xy(ii,2) - pos_xy(:,2)).^2));
    cvhm_rch.v(irch_neg(ii),1) = urf_rch_pos(ddd,1);
    dst(ii,1) = ccc;
end
Fcvhm_rch = scatteredInterpolant(cvhm_rch.p(:,1), cvhm_rch.p(:,2), cvhm_rch.v,'linear','nearest');
%% Create a matrix with the streamline URFs
% Fields of URFS 
% Eid Sid gnlm_index swat_index mu std wgh rch
URF_DATA = nan(100*length(wells),5); %Eid Sid mu std wgh
urf_xy = nan(100*length(wells),2);
cnt = 1;
for ii = 1:length(wells)
    for jj = 1:size(wells(ii).urfs,1)
        if ~isempty(wells(ii).urfs(jj).p_lnd)
            URF_DATA(cnt,:) = [ii jj wells(ii).urfs(jj).m wells(ii).urfs(jj).s wells(ii).urfs(jj).v_cds];
            urf_xy(cnt,:) = wells(ii).urfs(jj).p_lnd;
            cnt = cnt + 1;
        end 
    end
end
URF_DATA(cnt:end,:) = [];
urf_xy(cnt:end,:) = [];
urf_rch = Fcvhm_rch(urf_xy(:,1),urf_xy(:,2));
[urfI, urfJ] = calc_IJ_Mantis(urf_xy);
urf_linear_index = sub2ind([12863 7046],urfI, urfJ); 
%% Load the GNLM Loading
fid = fopen(fullfile('..','MantisData', 'GNLM_LU_NGW.dat'),'r');
GNLM_LOAD = textscan(fid, '%d %d %d %d %d %d %f %f %f %f %f %f %f %f');
fclose(fid);
%% Find the rows of the linear indices in the GNLM LOAD matrix
[Lia, urf_gnlm_index] = ismember(urf_linear_index, GNLM_LOAD{1,1});
%% create a shape file from the coordinates
S(size(urf_xy,1),1).Geometry = 'Point';
for ii = 1:size(urf_xy,1)
    S(ii,1).Geometry = 'Point';
    S(ii,1).X = urf_xy(ii,1);
    S(ii,1).Y = urf_xy(ii,2);
    S(ii,1).Eid = URF_DATA(ii,1);
    S(ii,1).Sid = URF_DATA(ii,2);
    S(ii,1).mu = URF_DATA(ii,3);
    S(ii,1).std = URF_DATA(ii,4);
    S(ii,1).wgh = URF_DATA(ii,5);
end
shapewrite(S, [scen_name '_URFdata'])
%% load the URF data with the HRUGIS fields for BUD0
URF_SAC = readtable(fullfile('CVHM_BUD0','UrfPntsSAC.csv'));
URF_SJV = readtable(fullfile('CVHM_BUD0','UrfPntsSJV.csv'));
URF_SSJV = readtable(fullfile('CVHM_BUD0','UrfPntsSSJV.csv'));
%% load the URF data with the HRUGIS fields for BUD1
URF_SAC = readtable(fullfile('CVHM_BUD1','UrfPntsSAC.csv'));
URF_SJV = readtable(fullfile('CVHM_BUD1','UrfPntsSJV.csv'));
URF_SSJV = readtable(fullfile('CVHM_BUD1','UrfPntsSSJV.csv'));
%% Load the SWAT loading.
% any scenario would do as the indices are identical
frmt = '%f %f %f';
for ii = 1:25
    frmt = [frmt ' %f'];
end
fid = fopen(fullfile('..','MantisData',['SWAT_LOADING_SCEN_' num2str(1) '.dat']),'r');
SWAT_LOAD = textscan(fid, frmt);
fclose(fid);
%% Identify the swat index in the shapefiles
[Lia, Locb] = ismember([URF_SAC.HRU_GIS1], SWAT_LOAD{1,2}-10000000);
SAC_swat_index = [[URF_SAC.Eid] [URF_SAC.Sid] Locb];
[Lia, Locb] = ismember([URF_SJV.HRU_GIS1], SWAT_LOAD{1,2}-10000000);
SJV_swat_index = [[URF_SJV.Eid] [URF_SJV.Sid] Locb];
[Lia, Locb] = ismember([URF_SSJV.HRU_GIS1], SWAT_LOAD{1,2}-10000000);
SSJV_swat_index = [[URF_SSJV.Eid] [URF_SSJV.Sid] Locb];
%% Find the location of the swat_index in the URF_DATA
urf_swat_index = nan(size(URF_DATA,1),1);
[Lia, Locb] = ismember(SAC_swat_index(:,1:2), URF_DATA(:,1:2), 'rows');
urf_swat_index(Locb,1) = SAC_swat_index(:,3);
[Lia, Locb] = ismember(SJV_swat_index(:,1:2), URF_DATA(:,1:2), 'rows');
urf_swat_index(Locb,1) = SJV_swat_index(:,3);
[Lia, Locb] = ismember(SSJV_swat_index(:,1:2), URF_DATA(:,1:2), 'rows');
urf_swat_index(Locb,1) = SSJV_swat_index(:,3);
urf_swat_index(isnan(urf_swat_index),1) = -9;
%% Assmble and print the URF file
% Eid Sid gnlm_index swat_index mu std wgh rch
URF_FINAL_MAT = [URF_DATA(:,1:2) urf_gnlm_index urf_swat_index URF_DATA(:,3:5) urf_rch*365*1000];
%% Write URF data
fid = fopen(fullfile('..','MantisData',['URFS_' scen_name '.dat']), 'w');
fprintf(fid, '%d %s LGNRM\n',size(URF_FINAL_MAT,1), scen_name);
fprintf(fid,'%d %d %d %d %.5f %.5f %.5f %.3f\n', URF_FINAL_MAT');
fclose(fid);
%% code for Mantis debug
ind = 16489831+1;
aa = [];
for ii = 1:14
    aa = [aa;GNLM_LOAD{1,ii}(ind)];
end