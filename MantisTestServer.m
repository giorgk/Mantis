function MantisTestServer()
disp('The main folder is:');
disp(pwd);
folder = pwd;
prefix = '/MantisServer.';
inputfile = [folder prefix 'inp'];
lockfile = [folder prefix 'lock'];
outfile = [folder prefix 'out'];
quitfile = [folder prefix 'quit'];

icounter = 0;
yrs = 1945:15:2050;
sim_yrs = 1945:2100;
Nsim_yrs = length(sim_yrs);

% ====== Load data ==================
disp('loading data...')
URFS = load('Local/CVHM/CVHM_ALLURFS_TA');
for jj = 1:length(yrs)
    Ngw{jj,1} = imread(['Local/Ngw_' num2str(yrs(jj)) '.tif']);
    Ngw{jj,1}(Ngw{jj,1} == Ngw{jj,1}(1,1)) = 0;
end

for jj = 1:5
    LUmaps{jj,1} = imread(['Local/model_input_LU' num2str(yrs(jj)) '.tif']);
end
load('Local/CVHM/CVHMWells');
load('MAPS');
%====================================



while true
    pause(1);
    icounter = icounter+1
    if exist(inputfile, 'file') == 2
        % create lock server file
        fid = fopen(lockfile,'w');
        fclose(fid);
        % read input file
        [MapID, RegIDs, LUinfo] = readInputfile(inputfile);
        
        % Pick wells based on the map selection
        wx = [Wells.X]';
        wy = [Wells.Y]';
        well_ids = [];
        for ii = 1:length(RegIDs)
            in = inpolygon(wx, wy, CVmap(MapID,1).data(RegIDs(ii),1).X,...
                                   CVmap(MapID,1).data(RegIDs(ii),1).Y);
            well_ids = [well_ids; [Wells(in,1).Eid]'];
        end

        Spnts_Eid = [URFS.URFS.Eid]';
        Spnts_X = [URFS.URFS.X]';
        Spnts_Y = [URFS.URFS.Y]';
        Spnts_V = [URFS.URFS.Vland]';
        
        Wellids = unique(Spnts_Eid);
        Wellids = intersect(well_ids, Wellids);
        well_sim_id = sort(Wellids);
        
        Spnt_sim_id = findElemAinB(well_sim_id, Spnts_Eid);
        Spnts_Eid = Spnts_Eid(Spnt_sim_id);
        Spnts_X = Spnts_X(Spnt_sim_id);
        Spnts_Y = Spnts_Y(Spnt_sim_id);
        Spnts_V = Spnts_V(Spnt_sim_id);
        tempurf = URFS.URFS(Spnt_sim_id,:);
        
        IJ = findIJ(Spnts_X, Spnts_Y);
        
        LFNC = zeros(length(Spnt_sim_id), Nsim_yrs);
        ALLURFS = zeros(length(Spnt_sim_id),200);
        
        for ii = 1:length(Spnt_sim_id)
            LF = nan(1,Nsim_yrs);
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
            end
            
            LF(k*15+1:end) = LF(k*15);
            LFNC(ii,:) = LF;
            
            if isempty(tempurf(ii).urf)
                ALLURFS(ii,:) = zeros(1,200);
            else
                %ALLURFS(ii,:) = tempurf(ii).urf; 
                ALLURFS(ii,:) = reBuildURF(tempurf(ii).urf.x, tempurf(ii).urf.y);
            end
        end
        
        BTC = ConvoluteURF(ALLURFS, LFNC, 'cpp');
        
        wells_btc = zeros(length(well_sim_id), size(LFNC,2));
        wgh = [tempurf.Vcds]';
        for ii = 1:length(well_sim_id)
            id = find(Spnts_Eid == well_sim_id(ii));
            if isempty(id)
                continue;
            end
            btc_temp = BTC(id,:);
            weight = wgh(id,1);
            weight = weight/sum(weight);
            btc_temp = bsxfun(@times,btc_temp,weight);
            wells_btc(ii,:) = sum(btc_temp,1);
        end
        
        writeoutfile(outfile, wells_btc);
        
        delete(inputfile);
        delete(lockfile);
    end
    
    if exist(quitfile, 'file') == 2
        delete(quitfile);
        break;
    end
end

end % END MAIN FUNCTION

function [MapID, RegIDs, LUinfo] = readInputfile(filename)
    fid = fopen(filename,'r');
    temp = fgetl(fid);
    C = textscan(temp, '%f');
    MapID = C{1,1}(1);
    if MapID == 1
        RegIDs = 1;
    else
        RegIDs = zeros(length(C{1,1})-1,1);
        for jj = 2:length(C{1,1})
            RegIDs(jj-1,1) = C{1,1}(jj);
        end
    end
    Ncrops = fscanf(fid, '%d',1);
    LUinfo = zeros(Ncrops, 2);
    for ii = 1:Ncrops
        lu = fscanf(fid, '%f',2)';
        LUinfo(ii,:) = lu;
    end
    fclose(fid);
end

function id = findElemAinB(A,B)
    id = [];
    for ii = 1:length(A)
        id = [id; find(B == A(ii))];
    end
end

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

function writeoutfile(outfile, BTC)
    fid = fopen(outfile,'w');
    fprintf(fid,'%d\n', size(BTC,1));
    frmt = '%f';
    for ii = 2:size(BTC,2)
        frmt = [frmt ' %f'];
    end
    fprintf(fid, [frmt '\n'], BTC');
    fclose(fid);
end

function urf = reBuildURF(x,y)
    urf = interp1(x,y,1:200);
end