function out = MainRun( LUinfo, ax )
%MainRun runs the convolution for a given Land use reduction
%   The input of the function is a N x 2 matrix. 
%   The first column is the land use id and 
%   the second column is the percentage of loading for the given category
%   100 % means no reduction

out = [];
yrs = 1945:15:2050;
ax = findobj('Tag','MainPlot');
% ====================LOAD DATA AREA=======================================
URFS = evalin('base','URFS');
Spnts = evalin('base','Spnts');
Ngw = evalin('base','Ngw');
LUmaps = evalin('base','LUmaps');

Wellids = unique([Spnts.Eid]');

% find the IJ in a vectorized fashion
IJ = findIJ([Spnts.X]', [Spnts.Y]');

LFNC = zeros(size(Spnts, 1),105);
LFNC_base = zeros(size(Spnts, 1),105);
ALLURFS = zeros(size(Spnts, 1),length(URFS.URFS(1).URF));

for ii = 1:size(Spnts, 1)
    % create the loading function
    LF = nan(1,105);
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
        
        % distribute the reduced loading 
        LF((k-1)*15+1:(k)*15) = linspace(Val_start*(red_s/100), ...
                                         Val_end*(red_e/100), 15);
        LF_base((k-1)*15+1:(k)*15) = linspace(Val_start, Val_end, 15);
    end
    LFNC(ii,:) = LF;
    LFNC_base(ii,:) = LF_base;
    ALLURFS(ii,:) = URFS.URFS(ii).URF;
end

BTC = ConvoluteURF(ALLURFS, LFNC, 'vect');
BTC_base = ConvoluteURF(ALLURFS, LFNC_base, 'vect');

wells_btc = zeros(length(Wellids), size(LFNC,2));
wells_btc_base = zeros(length(Wellids), size(LFNC,2));
for ii = 1:length(Wellids)
    ii
    % find the streamlines of well ii
    id = find([Spnts.Eid]' == Wellids(ii));
    if isempty(id)
        continue;
    end
    btc_temp = BTC(id,:);
    btc_temp_base = BTC_base(id,:);
    weight = [URFS.URFS(id,1).v_cds]';
    weight = weight/sum(weight);
    btc_temp = bsxfun(@times,btc_temp,weight);
    btc_temp_base = bsxfun(@times,btc_temp_base,weight);
    wells_btc(ii,:) = sum(btc_temp,1);
    wells_btc_base(ii,:) = sum(btc_temp_base,1);
end

perc = prctile(wells_btc,10:10:90,1);
perc_base = prctile(wells_btc_base,10:10:90,1);
plot(ax,perc','g');
hold on
plot(ax,perc_base','k');

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

