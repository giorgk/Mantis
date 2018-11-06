function out = MainRun(MAPS, LUinfo, RunTag, accuracy)
%MainRun runs the convolution for a given Land use reduction
% MAPS : a structure that containts information about the area to compute
%       the statistics
% LUinfo: This is a N x 2 matrix. 
%         The first column is the land use id and 
%         the second column is the percentage of loading for the given category
%         100 % means no reduction
% RunTag: A tag for a given run

out = [];
yrs = 1945:15:2050;
sim_yrs = 1945:2100;
Nsim_yrs = length(sim_yrs);
%ax = findobj('Tag','MainPlot');
hstat = findobj('Tag','Stats');
hcppmode = findobj('Tag', 'cpp_mode');
% ====================LOAD DATA AREA=======================================
set(hstat,'String', 'Boiling data...');
drawnow
URFS = evalin('base','URFS');
%Spnts = evalin('base','Spnts');
Ngw = evalin('base','Ngw');
LUmaps = evalin('base','LUmaps');

%Spnts_Eid = [Spnts.Eid]';
%Spnts_X = [Spnts.X]';
%Spnts_Y = [Spnts.Y]';
%Spnts_V = [Spnts.Vland]';
Spnts_Eid = [URFS.URFS.Eid]';
Spnts_X = [URFS.URFS.X]';
Spnts_Y = [URFS.URFS.Y]';
Spnts_V = [URFS.URFS.Vland]';

Wellids = unique(Spnts_Eid);

% choose a number of wells based on the accuracy
Nwells = ceil(max(10, length(Wellids)*accuracy));
well_rand = randperm(length(Wellids));
well_sim_id = Wellids(well_rand(1:Nwells));
well_sim_id = sort(well_sim_id);

% find the streamlines associated with the selected wells
Spnt_sim_id = findElemAinB(well_sim_id, Spnts_Eid);
Spnts_Eid = Spnts_Eid(Spnt_sim_id);
Spnts_X = Spnts_X(Spnt_sim_id);
Spnts_Y = Spnts_Y(Spnt_sim_id);
Spnts_V = Spnts_V(Spnt_sim_id);
tempurf = URFS.URFS(Spnt_sim_id,:);

% find the IJ in a vectorized fashion
IJ = findIJ(Spnts_X, Spnts_Y);

LFNC = zeros(length(Spnt_sim_id), Nsim_yrs);
ALLURFS = zeros(length(Spnt_sim_id),length(URFS.URFS(1).URF));
set(hstat,'String', 'Building Loading functions...');
drawnow
tic
for ii = 1:length(Spnt_sim_id)
    % create the loading function
    LF = nan(1,Nsim_yrs);
    %LF_base = nan(1,Nsim_yrs);
    % find pixel landuse
    for k = 1:length(yrs)-1
        Val_start = Ngw{k,1}(IJ(ii,1),IJ(ii,2));
        Val_end = Ngw{k+1,1}(IJ(ii,1),IJ(ii,2));
        
        klu_s = k; % index for land use map. it should go up to 5
        klu_e = k+1;
        if k >= 5
            klu_s = 5;
            klu_e = 5;
        end

        lu_s = double(LUmaps{klu_s,1}(IJ(ii,1),IJ(ii,2)));
        lu_e = double(LUmaps{klu_e,1}(IJ(ii,1),IJ(ii,2)));
        
        % find the reduction that coresponds to those land uses
        id_s = find(LUinfo == lu_s, 1);
        id_e = find(LUinfo == lu_e, 1);
        if isempty(id_s)
            red_s = 100;
        else
            red_s = LUinfo(id_s,2);
        end
        
        if isempty(id_e) 
            red_e = 100;
        else
            red_e = LUinfo(id_e,2);
        end
        
        % distribute the reduced loading after 2020
        if yrs(k) >= 2020
            LF((k-1)*15+1:k*15) = linspace(Val_start*(red_s/100), ...
                                             Val_end*(red_e/100), 15);
        else
            LF((k-1)*15+1:k*15) = linspace(Val_start, Val_end, 15);
        end
        %LF_base((k-1)*15+1:k*15) = linspace(Val_start, Val_end, 15);
    end
    
    % after 2050 assume constant loading
    LF(k*15+1:end) = LF(k*15);
    %LF_base(k*15+1:end) = LF_base(k*15);
    
    LFNC(ii,:) = LF;
    %LFNC_base(ii,:) = LF_base;
    %ALLURFS(ii,:) = tempurf(ii).URF;
    ALLURFS(ii,:) = reBuildURF(tempurf(ii).URF.urf.x, tempurf(ii).URF.urf.y);
end
time_lf = toc;
set(hstat,'String', 'Calculating BTC...');
drawnow
tic
if get(hcppmode,'Value') == 1
    BTC = ConvoluteURF(ALLURFS, LFNC, 'cpp');
else
    BTC = ConvoluteURF(ALLURFS, LFNC, 'vect');
end
%BTC_base = ConvoluteURF(ALLURFS, LFNC_base, 'vect');
time_bct = toc;

tic
wells_btc = zeros(length(well_sim_id), size(LFNC,2));
wgh = [tempurf.v_cds]';
for ii = 1:length(well_sim_id)
    % find the streamlines of well ii
    id = find(Spnts_Eid == well_sim_id(ii));
    if isempty(id)
        continue;
    end
    btc_temp = BTC(id,:);
    weight = wgh(id,1);%[URFS.URFS(id,1).v_cds]';
    weight = weight/sum(weight);
    btc_temp = bsxfun(@times,btc_temp,weight);
    wells_btc(ii,:) = sum(btc_temp,1);
end
time_stat = toc;
stat_str{1,1} = ['Stats: Lfnc : ' num2str(time_lf) ' sec'];
stat_str{1,2} = ['          BTC  : ' num2str(time_bct) ' sec'];
stat_str{1,3} = ['          Stat : ' num2str(time_stat) ' sec'];
set(hstat,'String', stat_str);

out.Tag = RunTag;
out.WellBTC = wells_btc;
out.wellid = well_sim_id;
return;

perc = prctile(wells_btc,10:10:90,1);
perc_base = prctile(wells_btc_base,10:10:90,1);
time_stat = toc;

plot(sim_yrs, perc','r', 'linewidth', 1.5);
hold on
plot(sim_yrs, perc_base','--k', 'linewidth', 1.5);
xlabel('Time[years]');
ylabel('Concentration [mg/l]');
xticks([1950:20:2100]);
xticklabels(datestr(datenum(1950:20:2100,1,1),'YY'))
xlim([sim_yrs(1) sim_yrs(end)]);
grid on
stat_str{1,1} = ['Stats: Lfnc : ' num2str(time_lf) ' sec'];
stat_str{1,2} = ['           BTC  : ' num2str(time_bct) ' sec'];
stat_str{1,3} = ['           Stat : ' num2str(time_stat) ' sec'];
set(hstat,'String', stat_str);

end

% ==================Sub functions==================================
function IJ = findIJ(x, y)
    % Returns the IJ of the cell for each coordinate x,y
    Xmin = -223300;
    Ymin = -344600;
    csz = 50; % cell size 
    Nr = 12863;
    Nc = 7046;
    Xgrid = Xmin:csz:Xmin+csz*Nc;
    Ygrid = Ymin:csz:Ymin+csz*Nr;
    IJ = nan(length(x),2);
    
    Nx = length(Xgrid);
    Ny = length(Ygrid);
    mx_dim = max(Nx-1, Ny-1);
    
    for it = 1:mx_dim
        if it < Nx
            ids = x >= Xgrid(it) & x <= Xgrid(it+1);
            IJ(ids,2) = it;
        end
        
        if it < Ny
            ids = y >= Ygrid(it) & y <= Ygrid(it+1);
            IJ(ids,1) = Nr - it+1;
        end
    end
end

function id = findElemAinB(A,B)
    id = [];
    for ii = 1:length(A)
        id = [id; find(B == A(ii))];
    end
end

function urf = reBuildURF(x,y)
    urf = interp1(x,y,1:200);
end


