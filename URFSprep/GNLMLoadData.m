top_level_path = fullfile('..','..');
yrs = 1945:15:2050;
%% Load GNLM loading and Land uses
clear LUs
for ii = 1:5
    num2str(yrs(ii))
    LUs{ii,1} = imread(fullfile(top_level_path,'SWAT_Model','NLoading','GNLM',['model_input_LU' num2str(yrs(ii)) '.tif']));  
end

for ii = 1:length(yrs)
    num2str(yrs(ii))
    GNLM{ii,1} = imread(fullfile(top_level_path,'SWAT_Model','NLoading','GNLM',['Ngw_' num2str(yrs(ii)) '.tif'])); 
    %GNLM{ii,1} = imread(['GNLM/Ngw_nondirect_' num2str(yrs(ii)) '.tif']);
    GNLM{ii,1}(GNLM{ii,1} == GNLM{ii,1}(1,1)) = 0;
end
%% Find the Central valley cells
[II, JJ] = find(LUs{1,1} ~=  0);
% linear index of central valley pixels
lind = sub2ind(size(LUs{1,1}),II,JJ);
LUlinear = nan(length(lind),7);
LUlinear(:,1) = (1:length(lind))';
LUlinear(:,2) = lind;
for ii = 1:5
    LUlinear(:,ii+2) = LUs{ii,1}(lind);
end
%% Write land uses
fid = fopen(fullfile('..','MantisData','GNLM_LUs.dat'), 'w');
fprintf(fid, '%d %d %d %d %d %d %d\n', LUlinear');
fclose(fid);
%% find the loading
NGW = nan(length(lind),10);
NGW(:,1) = (1:length(lind))';
NGW(:,2) = lind;
for ii = 1:8
    NGW(:, ii+2) = GNLM{ii,1}(lind);
end
%% Write N loading
fid = fopen(fullfile('..','MantisData','GNLM_NGW.dat'), 'w');
fprintf(fid, '%d %d %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n', NGW');
fclose(fid);
%% Combine LU and NGW
LU_NGW = [LUlinear NGW(:,3:10)];
LU_NGW(sum(LU_NGW(:,8:15),2) == 0,:) = [];
%% Write LU and loading combined
fid = fopen(fullfile('..','MantisData','GNLM_LU_NGW.dat'), 'w');
fprintf(fid, '%d %d %d %d %d %d %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n', LU_NGW(:,2:end)');
fclose(fid);