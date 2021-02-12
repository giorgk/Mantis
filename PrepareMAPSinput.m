%% Read background map data
bckgrMaps = {'CentralValley','CV_outline_simple';...
             'Basins','Basins_simple';...
             'Counties', 'counties_simple';...
             'B118', 'B118_simple';...
             'Townships', 'CVHM_Townships_3310_simplified';...
             'CVHMfarms', 'CVHM_FarmsTA'; ...
             'C2VsimSubregions', 'C2Vsim_Subregions_3310'};
for ii = 1:length(bckgrMaps)
    CVmap(ii,1).name = bckgrMaps{ii,1};
    CVmap(ii,1).data = shaperead(fullfile('gis_data', bckgrMaps{ii,2}));
end
%% 
fid = fopen(fullfile('MantisData','BackgroundMapsv1.dat'),'w');
fprintf(fid, '%d\n', length(CVmap)); % Number of Background maps
for ii = 1:length(CVmap)
    % Name and number of subregions this bacnkground maps is devided
    fprintf(fid, '%s %d\n', CVmap(ii,1).name, length(CVmap(ii,1).data));
    for jj = 1:length(CVmap(ii,1).data)
        [Xs, Ys] = polysplit(CVmap(ii,1).data(jj,1).X, CVmap(ii,1).data(jj,1).Y);
        is_not_hole = false(length(Xs),1);
        for k = 1:length(Xs)
            if ispolycw(Xs{k,1}, Ys{k,1})
                is_not_hole(k) = true;
            end
        end
        if ii == 1
            subregion_name = 'CentralValley';
        elseif ii == 2
            subregion_name = CVmap(ii,1).data(jj,1).CVHM_Basin;
        elseif ii == 3
            subregion_name = CVmap(ii,1).data(jj,1).name;
        elseif ii == 4
            subregion_name =  CVmap(ii,1).data(jj,1).Basin_Subb;
        elseif ii == 5
            subregion_name = CVmap(ii,1).data(jj,1).CO_MTR;
        elseif ii == 6
            subregion_name = ['Farm ' num2str(CVmap(ii,1).data(jj,1).dwr_sbrgns)];
        elseif ii == 7
            subregion_name = ['Subregion ' num2str(CVmap(ii,1).data(jj,1).IRGE)];
        end
        subregion_name(:,isspace(subregion_name)) = [];
        subregion_name = replace(subregion_name,{'-','.'},'_');
        
        % name of subregion and number of polygons
        fprintf(fid, '%s %d\n', subregion_name, sum(is_not_hole));
        for k = 1:length(is_not_hole)
            if ~is_not_hole(k)
                continue;
            end
            xx = Xs{k}(1); yy = Ys{k,1}(1);
            % remove duplicates
            for nn = 2:length(Xs{k,1})
                dst = min(sqrt((xx - Xs{k,1}(nn)).^2 + (yy - Ys{k,1}(nn)).^2));
                if dst > 0.1
                    xx = [xx; Xs{k,1}(nn)];
                    yy = [yy; Ys{k,1}(nn)];
                end
            end
            
            fprintf(fid, '%d\n', length(xx));
            fprintf(fid,'%.2f %.2f\n', [xx yy]');
        end
    end
end
fclose(fid);
