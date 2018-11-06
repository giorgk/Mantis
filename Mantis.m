function varargout = Mantis(varargin)
% MANTIS MATLAB code for Mantis.fig
%      MANTIS, by itself, creates a new MANTIS or raises the existing
%      singleton*.
%
%      H = MANTIS returns the handle to a new MANTIS or the handle to
%      the existing singleton*.
%
%      MANTIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MANTIS.M with the given input arguments.
%
%      MANTIS('Property','Value',...) creates a new MANTIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Mantis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Mantis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Mantis

% Last Modified by GUIDE v2.5 24-Sep-2018 00:50:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Mantis_OpeningFcn, ...
                   'gui_OutputFcn',  @Mantis_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Mantis is made visible.
function Mantis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Mantis (see VARARGIN)

% Choose default command line output for Mantis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Mantis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Mantis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
ha = axes('units','normalized', 'Tag','logo',...
            'position',[0.95 0.91 0.09*0.7 0.09]);
uistack(ha,'bottom');
I=imread('mantis1.jpeg');
hi = imagesc(I);
set(ha,'handlevisibility','off', ...
            'visible','off')
axis(ha,'equal')



% --- Executes on slider movement.
function LoadingSlider_Callback(hObject, eventdata, handles)
% hObject    handle to LoadingSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set_percentage(get(hObject,'Value')*100);






% --- Executes during object creation, after setting all properties.
function LoadingSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LoadingSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in LandUseList.
function LandUseList_Callback(hObject, eventdata, handles)
% hObject    handle to LandUseList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns LandUseList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LandUseList

%display(get(hObject,'Value'));
LU = evalin('base','LU');

hSelectionMode = findobj('Tag', 'CropSelectionMenu'); % get the crop selection object
IselMode = get(hSelectionMode,'Value'); % get the id of the selected grouping

ilu = get(hObject,'Value'); % get the id of the selected entry in the crop list

hlulabel = findobj('Tag','SelectedLU'); % get the uicontrol that displays the name of the crop/land


set(hlulabel, 'String', LU.class.classes(IselMode).Names{ilu}); % set the name

hluslider = findobj('Tag', 'LoadingSlider'); % get the slider object

% find the percentage that corresponds to the selected crop.
% if there are more that two crops in the selected ilu set the slider to
% the percentage of the first
ids_perc = LU.class.classes(IselMode).CAMLcodes{ilu}(1);

LU.selected = ilu;% maybe this not needed


set(hluslider, 'Value', LU.perc(ids_perc,1)/100);
percdisp = findobj('Tag','PercDisplay');
set(percdisp, 'String', num2str(LU.perc(ids_perc,1), '%.2f'));
assignin('base', 'LU', LU)



% --- Executes during object creation, after setting all properties.
function LandUseList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LandUseList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%load('LU_data');
%set(hObject, 'Items',LU_name);
%f = uifigure;
%l = uilistbox(f);

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%LU = load('LU_data');
LU = load('LU_data_v2');
% set up random percentages
LU.perc = 100*rand(length(LU.DWRCAMLCode), 1);
RUNS = [];
opt = setMantisoptions;
% Append percentage loading
LU.selected = -9;
assignin('base', 'LU', LU);
assignin('base', 'RUNS', RUNS);
assignin('base', 'opt', opt);

hSelectionMode = findobj('Tag', 'CropSelectionMenu');
set(hSelectionMode,'String', LU.class.GroupNames);
set(hSelectionMode, 'Value', 1);
set(hObject, 'String', LU.class.classes(1).Names);



function SelectedLU_Callback(hObject, eventdata, handles)
% hObject    handle to SelectedLU (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SelectedLU as text
%        str2double(get(hObject,'String')) returns contents of SelectedLU as a double


% --- Executes during object creation, after setting all properties.
function SelectedLU_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelectedLU (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PercDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to PercDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PercDisplay as text
%        str2double(get(hObject,'String')) returns contents of PercDisplay as a double
set_percentage(get(hObject,'String'));



% --- Executes during object creation, after setting all properties.
function PercDisplay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PercDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RunButton.
function RunButton_Callback(hObject, eventdata, handles)
% hObject    handle to RunButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

LU = evalin('base','LU');
opt = evalin('base', 'opt');
RUNS = evalin('base','RUNS');
MAPS = evalin('base','MAPS');

% set a tag
htagbox = findobj('Tag', 'NickNameRun');
hTagList = findobj('Tag', 'ListofRuns');
runtag = get(htagbox,'String');
if isempty(runtag) || strcmp('Set a tag for each run',runtag)
    runtag = 'Run 1';
end
% append the accuracy
runtag = [runtag '_{ac' num2str(round(opt.sim_accuracy*100)) '}'];

ListofRuns = cellstr(get(hTagList,'String'));
if length(ListofRuns) == 1
    if isempty(ListofRuns{1,1})
        ListofRuns{1,1} = runtag;
    else
        ListofRuns{length(ListofRuns)+1,1} = runtag;
    end
else
    ListofRuns{length(ListofRuns)+1,1} = runtag;
end




%ax = findobj('Tag','MainPlot');
out = MainRun(MAPS, [LU.DWRCAMLCode, LU.perc], runtag, opt.sim_accuracy );
if ~isempty(out)
    RUNS = [RUNS;out];
    assignin('base', 'RUNS', RUNS);
    set(hTagList, 'String', ListofRuns);
    set(hTagList, 'Value', length(ListofRuns));
    drawnow
end



function Stats_Callback(hObject, eventdata, handles)
% hObject    handle to Stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Stats as text
%        str2double(get(hObject,'String')) returns contents of Stats as a double


% --- Executes during object creation, after setting all properties.
function Stats_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadDataButton.
function LoadDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


hstat = findobj('Tag','Stats');
set(hstat,'String', 'Loading URFS...');
drawnow
yrs = 1945:15:2050;

% Unit Response Function data
%URFS = load('Local/Tule/ALLURFS'); 
%Spnts = shaperead('gis_data/TuleStrmlnPointsHome'); % URF end Points at the land side
URFS = load('Local/CVHM/CVHM_ALLURFS_TA');

assignin('base', 'URFS', URFS);
%assignin('base', 'Spnts', Spnts);


% Load Ngw
set(hstat,'String', 'Loading NGW maps...');
drawnow
for jj = 1:length(yrs)
    Ngw{jj,1} = imread(['Local/Ngw_' num2str(yrs(jj)) '.tif']);
    Ngw{jj,1}(Ngw{jj,1} == Ngw{jj,1}(1,1)) = 0;
end
assignin('base', 'Ngw', Ngw);
% Load Land use historic maps
set(hstat,'String', 'Loading LU maps...');
drawnow
for jj = 1:5
    LUmaps{jj,1} = imread(['Local/model_input_LU' num2str(yrs(jj)) '.tif']);
end

assignin('base', 'LUmaps', LUmaps);


% Load wells
%load('Local/Tule/TuleWells');
load('Local/CVHM/CVHMWells');

% load wells
set(hstat,'String', 'Loading Wells...');
drawnow
%Wells = shaperead('Local/Tule/TuleWells3310');
%load('Local/Tule/TuleWells');
assignin('base', 'Wells', Wells);


set(hstat,'String', 'Loading Done');
drawnow


% --- Executes on selection change in CropSelectionMenu.
function CropSelectionMenu_Callback(hObject, eventdata, handles)
% hObject    handle to CropSelectionMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CropSelectionMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CropSelectionMenu
SelectModeOptions = cellstr(get(hObject,'String'));
LU = evalin('base', 'LU');
hCropList = findobj('Tag', 'LandUseList');

set(hCropList, 'Value', 1);
set(hCropList, 'String', LU.class.classes(get(hObject, 'Value')).Names)

LU.selected = -9;
assignin('base', 'LU', LU);



% --- Executes during object creation, after setting all properties.
function CropSelectionMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CropSelectionMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function set_percentage(val)
% get data and uicontrols
LU = evalin('base','LU');
hstat = findobj('Tag','Stats');
hslider = findobj('Tag','LoadingSlider');
hpercdisp = findobj('Tag','PercDisplay');
hSelectmode = findobj('Tag', 'CropSelectionMenu');
hCropList = findobj('Tag', 'LandUseList');
%SelectModeOptions = cellstr(get(hSelectmode,'String'));

if LU.selected < 0
    set(hstat,'String', 'Select a first a category');
else
    if ~isnumeric(val)
        val = str2double(val);
    end
    if isnan(val)
        set(hstat,'String', 'Input value must be numeric only (dont use (,)!');
    else
        if val >= 0 && val <= 100
            iSelMode = get(hSelectmode, 'Value');
            ilu = get(hCropList, 'Value');
            
            % find the ids that correspond to the selection
            ids = LU.class.classes(iSelMode).CAMLcodes{ilu};  
            
            LU.perc(ids, 1) = val;
            set(hslider, 'Value', val/100);
            set(hpercdisp, 'String', num2str(val, '%.2f'));
            set(hstat,'String', '');
            assignin('base', 'LU', LU)
            
            %%% propagate the change to the finer categories
            %id = LU.LU_groups(:,iSelMode) == LU.selected;
            %if iSelMode > 1
            %    for ii = iSelMode-1:-1:1
            %        LU.perc(id, ii) = val;
            %    end
            %end
            
        else
            set(hstat,'String', ['[' num2str(val) '] is not between (0,100)']);
        end
    end
end



function NickNameRun_Callback(hObject, eventdata, handles)
% hObject    handle to NickNameRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NickNameRun as text
%        str2double(get(hObject,'String')) returns contents of NickNameRun as a double


% --- Executes during object creation, after setting all properties.
function NickNameRun_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NickNameRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ListofRuns.
function ListofRuns_Callback(hObject, eventdata, handles)
% hObject    handle to ListofRuns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ListofRuns contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ListofRuns


% --- Executes during object creation, after setting all properties.
function ListofRuns_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ListofRuns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PlotingStylesMenu.
function PlotingStylesMenu_Callback(hObject, eventdata, handles)
% hObject    handle to PlotingStylesMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PlotingStylesMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PlotingStylesMenu


% --- Executes during object creation, after setting all properties.
function PlotingStylesMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlotingStylesMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

set(hObject,'String', {'10:10:90','5 50 95','Mean'})
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PlotButton.
function PlotButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get the data
RUNS = evalin('base','RUNS');
% get the GUI elements
hScenarioList = findobj('Tag', 'ListofRuns');
hPlotType = findobj('Tag', 'PlotingStylesMenu');
hPlot = findobj('Tag', 'MainPlot');

iscenario = get(hScenarioList,'Value');
ScenarioNames = cellstr(get(hScenarioList,'String'));
SelectedScenario = ScenarioNames{iscenario};

opt = evalin('base', 'opt');

for ii = 1:size(RUNS,1)
    if strcmp(SelectedScenario, RUNS(ii,1).Tag)
        PlotTypes = cellstr(get(hPlotType,'String'));
        selectedPlotType = PlotTypes{get(hPlotType, 'Value')};
        
        if strcmp(selectedPlotType, '10:10:90')
            perc = prctile(RUNS(ii,1).WellBTC, 10:10:90, 1);
            
        elseif strcmp(selectedPlotType, '5 50 95')
            perc = prctile(RUNS(ii,1).WellBTC, [5 50 95], 1);
            
        elseif strcmp(selectedPlotType, 'Mean')
            perc = mean(RUNS(ii,1).WellBTC,1);
            
        end
        
        sim_yrs = 1945:2100;
        hold(hPlot, 'on');
        linesPlot = plot(hPlot, sim_yrs, perc', 'color', opt.clr_list(opt.clr_id,:), ...
                        'linewidth', 1.5); %, 'DisplayName', RUNS(ii,1).Tag
        linesGroup = hggroup(hPlot,'DisplayName', RUNS(ii,1).Tag);
        set(linesPlot,'Parent',linesGroup);
        set(get(get(linesGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
        
        opt = updateColorIndex(opt);
        set(get(hPlot,'XLabel'), 'String', 'Time[years]');
        set(get(hPlot,'YLabel'), 'String', 'Concentration [mg/l]');
        set(hPlot, 'XTick', [1950:20:2100]);
        set(hPlot, 'XTickLabel', datestr(datenum(1950:20:2100,1,1),'YY'));
        set(hPlot, 'XLim', [sim_yrs(1) sim_yrs(end)]);
        legend(hPlot,'Location', 'northwest');
        grid(hPlot, 'on');
        %xticks([1950:20:2100]);
        %xticklabels(datestr(datenum(1950:20:2100,1,1),'YY'))
        %xlim([sim_yrs(1) sim_yrs(end)]);
        %grid on
    end
end
assignin('base', 'opt', opt);


% --- Executes on button press in ClearScenarioListButton.
function ClearScenarioListButton_Callback(hObject, eventdata, handles)
% hObject    handle to ClearScenarioListButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hListofRuns = findobj('Tag', 'ListofRuns');
set(hListofRuns, 'String', {});
set(hListofRuns, 'Value', []);
RUNS = evalin('base', 'RUNS');
RUNS = [];
assignin('base', 'RUNS', RUNS);

function opt = setMantisoptions
hPlot = findobj('Tag', 'MainPlot');
opt.clr_list = get(hPlot, 'colororder');
opt.clr_id = 1;

hAccSlider = findobj('Tag', 'AccuracySlider');
accuracy = get(hAccSlider,'Value');
opt.sim_accuracy = accuracy;

function opt = updateColorIndex(opt)
opt.clr_id = opt.clr_id + 1;
if opt.clr_id > size(opt.clr_list,1)
    opt.clr_id = 1;
end

% --- Executes on button press in ClearPlotButton.
function ClearPlotButton_Callback(hObject, eventdata, handles)
% hObject    handle to ClearPlotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hPlot = findobj('Tag', 'MainPlot');
cla(hPlot);
opt = evalin('base', 'opt');
opt.clr_id = 1;
assignin('base', 'opt', opt);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over ClearScenarioListButton.
function ClearScenarioListButton_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ClearScenarioListButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on slider movement.
function AccuracySlider_Callback(hObject, eventdata, handles)
% hObject    handle to AccuracySlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
accuracy = get(hObject,'Value');
opt = evalin('base', 'opt');
opt.sim_accuracy = accuracy;
assignin('base', 'opt', opt);




% --- Executes during object creation, after setting all properties.
function AccuracySlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AccuracySlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in cpp_mode.
function cpp_mode_Callback(hObject, eventdata, handles)
% hObject    handle to cpp_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cpp_mode


% --- Executes on selection change in SpatialSelection.
function SpatialSelection_Callback(hObject, eventdata, handles)
% hObject    handle to SpatialSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SpatialSelection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SpatialSelection
imap = get(hObject,'Value');
SetMap(imap);



% --- Executes during object creation, after setting all properties.
function SpatialSelection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SpatialSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

MAPS = load('MAPS.mat');

set(hObject,'String', {MAPS.CVmap.name});
set(hObject,'Value', 1);
MAPS.imap = 1;
MAPS.SelPoly = false(1,1);
assignin('base', 'MAPS', MAPS);
SetMap(1);

%hPlotMap = findobj('Tag', 'MapPlot');
%plot(hPlotMap, MAPS.CVmap(1,1).data.X, MAPS.CVmap(1,1).data.Y, 'linewidth',2);
%axis(hPlotMap,'off','equal');

function SetMap(imap)
MAPS = evalin('base','MAPS');
if imap > length(MAPS.CVmap)
    warndlg(['You tried to set map ' num2str(imap) ' but there are only ' num2str(length(MAPS.CVmap)) ' maps in the list']);
else
    hPlotMap = findobj('Tag','MapPlot');
    
    cla(hPlotMap);
    hold(hPlotMap, 'on');
    for ii = 1:length(MAPS.CVmap(imap).data)
        plot(hPlotMap, MAPS.CVmap(imap).data(ii,1).X, MAPS.CVmap(imap).data(ii,1).Y, 'color', [0 0.4470 0.7410], 'linewidth', 1.5);
    end
    axis(hPlotMap, 'equal','off')
    
    hold(hPlotMap, 'off');
end

MAPS.imap = imap;
MAPS.SelPoly = false(length(MAPS.CVmap(imap).data),1);
assignin('base', 'MAPS', MAPS);


% --- Executes on button press in PointSelect.
function PointSelect_Callback(hObject, eventdata, handles)
% hObject    handle to PointSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hPlotMap = findobj('Tag','MapPlot');
axes(hPlotMap);

MAPS = evalin('base','MAPS');

while 1
    [x, y, button] = ginput(1);
    if button ~= 1
        break;
    end
    p = [x y];

    %get(gcf,'SelectionType')

    
    imap = MAPS.imap;
    for ii = 1:length(MAPS.CVmap(imap).data)
        in = inpolygon(p(1), p(2), MAPS.CVmap(imap).data(ii,1).X, MAPS.CVmap(imap).data(ii,1).Y);
        if in
            if MAPS.SelPoly(ii)
                MAPS.SelPoly(ii)= false;
            else
                MAPS.SelPoly(ii)= true;
            end
           cla(hPlotMap);
           hold(hPlotMap, 'on');
           for jj = 1:length(MAPS.CVmap(imap).data)
               if ~MAPS.SelPoly(jj)
                    plot(hPlotMap, MAPS.CVmap(imap).data(jj,1).X, MAPS.CVmap(imap).data(jj,1).Y, 'color', [0 0.4470 0.7410], 'linewidth', 1.5);
               else
                   [Xs, Ys] = polysplit(MAPS.CVmap(imap).data(jj,1).X, MAPS.CVmap(imap).data(jj,1).Y);
                   for kk = 1:length(Xs)
                       if ispolycw(Xs{kk,1}, Ys{kk,1})
                            patch(Xs{kk,1}, Ys{kk,1},[0.8500 0.3250 0.0980], 'FaceAlpha',0.5'); %, 'FaceColor',[0.8500 0.3250 0.0980]
                       end
                   end
               end
           end
           hold(hPlotMap, 'off');
           break;
        end
    end
end
assignin('base', 'MAPS', MAPS);
