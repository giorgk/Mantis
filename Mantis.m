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

% Last Modified by GUIDE v2.5 05-Aug-2018 22:54:17

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
            'position',[0.94 0.91 0.09*0.7 0.09]);
uistack(ha,'bottom');
I=imread('mantis1.jpeg');
hi = imagesc(I);
set(ha,'handlevisibility','off', ...
            'visible','off')


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
hSelectionMode = findobj('Tag', 'CropSelectionMenu');
IselMode = get(hSelectionMode,'Value');

ilu = get(hObject,'Value');
lulabel = findobj('Tag','SelectedLU');
LU = evalin('base','LU');
LU.selected = ilu;
set(lulabel, 'String', LU.LU_name{ilu, IselMode});
luslider = findobj('Tag', 'LoadingSlider');
set(luslider, 'Value', LU.perc(ilu,IselMode)/100);
percdisp = findobj('Tag','PercDisplay');
set(percdisp, 'String', num2str(LU.perc(ilu,IselMode), '%.2f'));
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
LU = load('LU_data_test');
RUNS = [];
opt = setMantisoptions;
% Append percentage loading
LU.perc = 100*rand(size(LU.LU_name));
LU.selected = -9;
assignin('base', 'LU', LU);
assignin('base', 'RUNS', RUNS);
assignin('base', 'opt', opt);

hSelectionMode = findobj('Tag', 'CropSelectionMenu');
set(hSelectionMode,'String', LU.groupNames);
set(hObject, 'String', {LU.LU_name{1:LU.Ncat(1),1}}');
set(hObject, 'Value', 1);
%h=findobj('Tag','LandUseList');
%h.Items = LU_name';
% =============== LOAD data==============================
if false
    yrs = 1945:15:2050;
    URFS = load('Local/Tule/ALLURFS'); % Unit Response Function data
    Spnts = shaperead('gis_data/TuleStrmlnPointsHome'); % URF end Points at the land side
    assignin('base', 'URFS', URFS);
    assignin('base', 'Spnts', Spnts);
    % Load Ngw
    for jj = 1:length(yrs)
        Ngw{jj,1} = imread(['Local/Ngw_' num2str(yrs(jj)) '.tif']);
        Ngw{jj,1}(Ngw{jj,1} == Ngw{jj,1}(1,1)) = 0;
    end
    assignin('base', 'Ngw', Ngw);
    % Load Land use historic maps
    for jj = 1:5
        LUmaps{jj,1} = imread(['Local/model_input_LU' num2str(yrs(jj)) '.tif']);
    end
    assignin('base', 'LUmaps', LUmaps);
end


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
out = MainRun( [LU.LU_cat, LU.perc], runtag, opt.sim_accuracy );
RUNS = [RUNS;out];
assignin('base', 'RUNS', RUNS);
set(hTagList, 'String', ListofRuns);
set(hTagList, 'Value', length(ListofRuns));
drawnow



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
set(hstat,'String', 'Loading Data...');
drawnow
yrs = 1945:15:2050;
URFS = load('Local/Tule/ALLURFS'); % Unit Response Function data
Spnts = shaperead('gis_data/TuleStrmlnPointsHome'); % URF end Points at the land side
assignin('base', 'URFS', URFS);
assignin('base', 'Spnts', Spnts);
% Load Ngw
for jj = 1:length(yrs)
    Ngw{jj,1} = imread(['Local/Ngw_' num2str(yrs(jj)) '.tif']);
    Ngw{jj,1}(Ngw{jj,1} == Ngw{jj,1}(1,1)) = 0;
end
assignin('base', 'Ngw', Ngw);
% Load Land use historic maps
for jj = 1:5
    LUmaps{jj,1} = imread(['Local/model_input_LU' num2str(yrs(jj)) '.tif']);
end
assignin('base', 'LUmaps', LUmaps);
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

for ii = 1:length(LU.groupNames)
    if strcmp(SelectModeOptions{get(hObject, 'Value')}, LU.groupNames{ii})
        set(hCropList, 'String', {LU.LU_name{1:LU.Ncat(ii),ii}}')
        set(hCropList, 'Value', 1);
        break;
    end
end


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
            LU.perc(LU.selected, iSelMode) = val;
            set(hslider, 'Value', LU.perc(LU.selected, iSelMode)/100);
            set(hpercdisp, 'String', num2str(LU.perc(LU.selected, iSelMode), '%.2f'));
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
        linesGroup = hggroup('DisplayName', RUNS(ii,1).Tag);
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
