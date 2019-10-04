%% Create MAPS.mat file
Names = {'Outline', 'Basins', 'Counties', 'B118', 'CVHM Farms'};
shps = {'CV_outline_simple', 'Basins_simple', 'counties_simple', 'B118_simple', 'CVHM_FarmsTA'};
for ii = 1:length(shps)
    CVmap(ii,1).name = Names{ii};
    CVmap(ii,1).data = shaperead(['gis_data/' shps{ii}]);
end
%% Make a unique list of land uses
%
LU = imread('Local/model_input_LU2005.tif');
LU_cat = unique(LU);
%% Find names for each LU category
% temp1 has all the numerical values and temp2 the text values
[temp1, temp2] = xlsread('Local/LanduseTable_2017_0515.xlsx', 'FINAL Landuse Table','A2:E208');
for ii = 1:length(LU_cat)
   id = find( temp1(:,1) == LU_cat(ii));
   LU_name{ii,1} = temp2{id,1};
end
%% save this to a mat file
save('LU_data', 'LU_cat','LU_name');
%% Create an Raster Ascii for GIS. This is just for test purposes
% from the Local/Ngw_2005.tif.xml it appears that the 
% left lower corner of the raster is at -223300, -344600
% and that the cell size is 50 m
% The coordinate system is the EPSG: 3310
WriteAscii4Raster('Local/LU_2005_ascii',LU, -223300, -344600, 50, 0);
%% Load URF data and make sure they are in the same coordinate system
% First make one variable with all streamline points
URFS = [];
for ii = 1:6 % This is the number of processors used in the simulation
    % w = load(['Local/Tule/wellURFS_000' num2str(ii-1) '.mat']);
    w = load(['Local/Tule/TuleRiverURFs_' num2str(ii-1) '.mat']);
    URFS = [URFS; w.WellURF];
end
%% Create a shapefile with the streamlines points at the land side.
% This shape file will be overlaid onto raster and convert the coordinates
% The coordinates of this shapefile are in EPSG:26911 
clear S
S = [];
S(size(URFS,1), 1).Geometry = [];
S(size(URFS,1), 1).X = [];
S(size(URFS,1), 1).Y = [];
S(size(URFS,1), 1).Eid = [];
S(size(URFS,1), 1).Sid = [];
for ii = 1:size(URFS,1)
   S(ii,1).Geometry = 'Point';
   S(ii,1).X = URFS(ii,1).p_lnd(1);
   S(ii,1).Y = URFS(ii,1).p_lnd(2);
   S(ii,1).Eid = double(URFS(ii,1).Eid);
   S(ii,1).Sid = double(URFS(ii,1).Sid);
   S(ii,1).Vland = URFS(ii,1).v_lnd;
end
shapewrite(S,'Local/Tule/TuleStrmlnPoints');
%% load the converted shapefile
% The converted shapefile has coordinates on EPSG:3310
% S = shaperead('gis_data/TuleStrmlnPoints');
S = shaperead('gis_data/TuleStrmlnPointsHome');
%% Save all data into one file for loading from python
years = 1945:15:2050;
for ii = 1:8
    if ii <= 5
        eval(['LU' num2str(years(ii)) ' = LUmaps{' num2str(ii) ',1};']);
    end
    eval(['Ngw' num2str(years(ii)) ' = Ngw{' num2str(ii) ',1};']);
end
%%
LUcat = LU.LU_cat;
save('data4python.mat', 'LUcat', 'LU1945','Ngw1945','-v7');
for ii = 2:8
    if ii <= 5
        save('data4python.mat',['LU' num2str(years(ii))], ['Ngw' num2str(years(ii))], '-append');
    else
        save('data4python.mat', ['Ngw' num2str(years(ii))], '-append');
    end
    
end
%% 
Spnts = shaperead('gis_data/TuleStrmlnPointsHome');
Sxyv = [[Spnts.X]' [Spnts.Y]' [Spnts.Vland]'];
Sid =  [[Spnts.Eid]' [Spnts.Sid]'];
urfs = zeros(size(URFS.URFS,1),200);
for ii = 1:size(URFS.URFS,1)
    urfs(ii,:) = URFS.URFS(ii,1).URF;
end
urfV = [URFS.URFS.v_lnd]';
save('URFdata.mat', 'Sxyv', 'Sid', 'urfs', 'urfV', '-v7');
%% Write dummy land categories
load('LU_data');
LU_name{1,2} = 'Trees';
LU_name{2,2} = 'Bushes';
LU_name{3,2} = 'Dairy lands';
LU_name{4,2} = 'Urban';
LU_name{5,2} = 'Native vegetation';
LU_name{6,2} = 'Rainforests';
LU_name{1,3} = 'Agricultural land';
LU_name{2,3} = 'Urban';
LU_name{3,3} = 'Native vegetation';
LU_name{1,4} = 'All land uses';
LU_groups(:,1) = [1:100]';
LU_groups(:,2) = ceil(6*rand(100,1));
LU_groups(:,3) = ceil(3*rand(100,1));
LU_groups(:,4) = 1;
Ncat = [100;6;3;1];
groupNames = {'Individual land use', 'Sub Groups', 'Super Groups',  'All land uses'};
save('LU_data_test', 'LU_cat', 'LU_name', 'LU_groups', 'groupNames','Ncat');
%% Land use groups based on the excel
% import the excel land use file as column vectors
for ii = 1:size(LanduseTable20170515,1)
    LUType{ii,1} = LanduseTable20170515{ii,1}{1};
    DWRCAMLCode(ii,1) = LanduseTable20170515{ii,2};
    CropLandGroup{ii,1} = LanduseTable20170515{ii,6}{1};
end
%}
%%
clear class
class.GroupNames{1,1} = 'Individual Crop/Land';
class.GroupNames{2,1} = 'Group Crop/Land';
class.GroupNames{3,1} = 'All Crop/Land';
class.classes(1,1).Names =  unique(LUType);
class.classes(2,1).Names =  unique(CropLandGroup);
class.classes(3,1).Names =  {'All Crops/Lands'};

%%
for ii = 1:length(class.classes(1,1).Names)
    class.classes(1,1).CAMLcodes{ii,1} = [];
    for jj = 1:length(LUType)
        if strcmp(LUType{jj,1}, class.classes(1,1).Names{ii,1})
            class.classes(1,1).CAMLcodes{ii,1} = ...
                [class.classes(1,1).CAMLcodes{ii,1}; jj];
        end
        
    end
end

for ii = 1:length(class.classes(2,1).Names)
    class.classes(2,1).CAMLcodes{ii,1} = [];
    for jj = 1:length(CropLandGroup)
        if strcmp(CropLandGroup{jj,1}, class.classes(2,1).Names{ii,1})
            class.classes(2,1).CAMLcodes{ii,1} = ...
                [class.classes(2,1).CAMLcodes{ii,1}; jj];
        end
        
    end
end
class.classes(3,1).CAMLcodes{1,1} = [1:length(DWRCAMLCode)]';
%%
save('LU_data_v2','LUType', 'DWRCAMLCode', 'class');
%% Export land use categories as javascript variable
% print to javascript map variable the cropMap. 
% An association between CAMLcode and group name
fid = fopen('tmp.tmp','w');
fprintf(fid, 'var cropMap = {};\n');
for ii = 1:length(LU.DWRCAMLCode)
    fprintf(fid, 'cropMap[%d] = {percentage: %0.2f, name: ''%s''};\n', LU.DWRCAMLCode(ii,1), 1, LU.LUType{ii,1});  
end
fclose(fid);
%% export as ee.Dictionary
fid = fopen('tmp.tmp','w');
fprintf(fid, 'var cropDict = ee.Dictionary.fromLists([\n');
for ii = 1:length(LU.DWRCAMLCode)
    fprintf(fid, '''%d'',', LU.DWRCAMLCode(ii,1));
    if mod(ii,20) == 0
        fprintf(fid, '\n');
    end
end
fprintf(fid, '],[\n');
for ii = 1:length(LU.DWRCAMLCode)
    fprintf(fid, '{percentage: %0.2f, name: ''%s''},\n', 1, LU.LUType{ii,1});
end
fprintf(fid, ']);\n');

fclose(fid);

%% Create a list of groups
fid = fopen('tmp1.tmp','w');
fprintf(fid, 'var cropGroups = {};\n');
for ii = [3 2 1]
   fprintf(fid, 'cropGroups[%d] = {className: ''%s'', groupNames: [\n', ii, LU.class.GroupNames{ii,1});
   for jj = 1:length(LU.class.classes(ii,1).Names)
       fprintf(fid, '''%s'', ', LU.class.classes(ii,1).Names{jj,1});
       if mod(jj, 3) == 0
           fprintf(fid, '\n');
       end
   end
   fprintf(fid, '], codes : [\n');
   for jj = 1:length(LU.class.classes(ii,1).Names)
       fprintf(fid, '[');
       for kk = 1:length(LU.class.classes(ii, 1).CAMLcodes{jj, 1})
          fprintf(fid, '%d,', LU.DWRCAMLCode(LU.class.classes(ii, 1).CAMLcodes{jj, 1}(kk)));
          if mod(kk,20) == 0
              fprintf(fid, '\n');
          end
       end
       fprintf(fid,'],\n');
   end
   fprintf(fid, ']};\n');
end

fclose(fid);
%% As dictionary
fid = fopen('tmp1.tmp','w');
fprintf(fid, 'var cropGroups = ee.Dictionary.fromLists(\n');
fprintf(fid, '[''1'', ''2'', ''3''],[\n');
for ii = [3 2 1]
    fprintf(fid, '{names:[\n');
    for jj = 1:length(LU.class.classes(ii,1).Names)
        fprintf(fid, '''%s'',\n', LU.class.classes(ii,1).Names{jj});
    end
    fprintf(fid, '],\n');
    fprintf(fid, 'groups:[\n');
    for jj = 1:length(LU.class.classes(ii,1).CAMLcodes)
        fprintf(fid, '[');
        for kk = 1:length(LU.class.classes(ii,1).CAMLcodes{jj,1})
            fprintf(fid, '%d,', LU.DWRCAMLCode(LU.class.classes(ii,1).CAMLcodes{jj,1}(kk)));
            if mod(kk, 20) == 0
                fprintf(fid, '\n');
            end
        end
        fprintf(fid, '],\n');
    end
    fprintf(fid, ']},\n');
end
fprintf(fid, ']\n');
fprintf(fid, ');\n');
fclose(fid);
%% load the image maps and save them as ascii files
yrs = 1945:15:2050;
for jj = 1:length(yrs)
    Ngw{jj,1} = imread(['Local/Ngw_' num2str(yrs(jj)) '.tif']);
    Ngw{jj,1}(Ngw{jj,1} == Ngw{jj,1}(1,1)) = 0;
end
%% find pixels with non zero values for at least one year
IJ = [];
for jj = 1:length(yrs)
    jj
    [I, J] = find(Ngw{jj,1} ~=0 );
    IJ = unique([IJ;[ I J]], 'rows');
end
%% Write them as one matrix
NGWs = [];
ind = sub2ind(size(Ngw{1,1}), IJ(:,1), IJ(:,2));
for jj = 1:length(yrs)
    NGWs = [NGWs Ngw{jj,1}(ind)];
end
%% Write NGws as ascii file
fid = fopen('Local/NGWs.dat','w');
fprintf(fid,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n', NGWs');
fclose(fid);
%% load Land uses
for jj = 1:5
    LUmaps{jj,1} = imread(['Local/model_input_LU' num2str(yrs(jj)) '.tif']);
end
%% Isolate only the non zero NGWs pixel for the landuses
LUs = [];
for jj = 1:5
    LUs = [LUs LUmaps{jj,1}(ind)];
end
%% Write LUs as ascii file
fid = fopen('Local/LUs.dat','w');
fprintf(fid,'%d %d %d %d %d %d %d\n', [IJ LUs]');
fclose(fid);
%% write them as binary
fileID = fopen('Local/LUs.bin','w');
fwrite(fileID,[IJ LUs],'uint');
fclose(fileID);

fileID = fopen('Local/NGWs.bin','w');
fwrite(fileID, NGWs,'double');
fclose(fileID);
%% Write MAPS for C++ server input
load('MAPS.mat')
fid = fopen('Local/MantisMaps.dat','w');
fprintf(fid, '%d\n', length(CVmap));
for ii = 1:length(CVmap)
    fprintf(fid, '%d\n', length(CVmap(ii,1).data));
    for jj = 1:length(CVmap(ii,1).data)
        [Xs, Ys] = msim_polysplit(CVmap(ii,1).data(jj,1).X, CVmap(ii,1).data(jj,1).Y);
        fprintf(fid, '%d\n', length(Xs) );
        for kk = 1:length(Xs)
            xx = Xs{kk}(1); yy = Ys{kk,1}(1);
            % remove duplicates
            for nn = 2:length(Xs{kk,1})
                dst = min(sqrt((xx - Xs{kk,1}(nn)).^2 + (yy - Ys{kk,1}(nn)).^2));
                if dst > 0.1
                    xx = [xx; Xs{kk,1}(nn)];
                    yy = [yy; Ys{kk,1}(nn)];
                end
            end
            
            fprintf(fid, '%d\n', length(xx) );
            for nn = 1:length(xx)
                fprintf(fid,'%f %f\n', [xx(nn) yy(nn)]);
            end
        end
    end
end

fclose(fid);



