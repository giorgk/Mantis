function varargout = TestGUI(varargin)
% TESTGUI MATLAB code for TestGUI.fig
%      TESTGUI, by itself, creates a new TESTGUI or raises the existing
%      singleton*.
%
%      H = TESTGUI returns the handle to a new TESTGUI or the handle to
%      the existing singleton*.
%
%      TESTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TESTGUI.M with the given input arguments.
%
%      TESTGUI('Property','Value',...) creates a new TESTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TestGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TestGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TestGUI

% Last Modified by GUIDE v2.5 21-Jul-2018 04:18:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TestGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @TestGUI_OutputFcn, ...
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


% --- Executes just before TestGUI is made visible.
function TestGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TestGUI (see VARARGIN)

% Choose default command line output for TestGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TestGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TestGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
ha = axes('units','normalized', 'Tag','logo',...
            'position',[0.92 0.89 0.11*0.7 0.11]);
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
newVal = get(hObject,'Value');
LU = evalin('base','LU');
LU.perc(LU.selected,1) = newVal*100;
percdisp = findobj('Tag','PercDisplay');
set(percdisp, 'String', num2str(LU.perc(LU.selected), '%.2f'));
assignin('base', 'LU', LU);





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
ilu = get(hObject,'Value');
lulabel = findobj('Tag','SelectedLU');
LU = evalin('base','LU');
LU.selected = ilu;
set(lulabel, 'String', LU.LU_name{ilu});
luslider = findobj('Tag', 'LoadingSlider');
set(luslider, 'Value', LU.perc(ilu)/100);
percdisp = findobj('Tag','PercDisplay');
set(percdisp, 'String', num2str(LU.perc(ilu), '%.2f'));
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

LU = load('LU_data');
% Append percentage loading
LU.perc = 100*rand(length(LU.LU_cat),1);
LU.selected = -9;
assignin('base', 'LU', LU)

set(hObject, 'String', LU.LU_name);
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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LU = evalin('base','LU');
%ax = findobj('Tag','MainPlot');
out = MainRun( [LU.LU_cat, LU.perc] );



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


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
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