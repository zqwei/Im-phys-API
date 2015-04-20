function varargout = Load_DataSets(varargin)
% LOAD_DATASETS MATLAB code for Load_DataSets.fig
%      LOAD_DATASETS, by itself, creates a new LOAD_DATASETS or raises the existing
%      singleton*.
%
%      H = LOAD_DATASETS returns the handle to a new LOAD_DATASETS or the handle to
%      the existing singleton*.
%
%      LOAD_DATASETS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOAD_DATASETS.M with the given input arguments.
%
%      LOAD_DATASETS('Property','Value',...) creates a new LOAD_DATASETS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Load_DataSets_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Load_DataSets_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Load_DataSets

% Last Modified by GUIDE v2.5 16-Jan-2015 12:06:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Load_DataSets_OpeningFcn, ...
                   'gui_OutputFcn',  @Load_DataSets_OutputFcn, ...
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


% --- Executes just before Load_DataSets is made visible.
function Load_DataSets_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Load_DataSets (see VARARGIN)

% Choose default command line output for Load_DataSets
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% Comparison_V1_0

% UIWAIT makes Load_DataSets wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Load_DataSets_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

if strcmp(handles.pushbuttonLoadData.Enable,'off') && ...
        strcmp(handles.listboxDataSets.Enable,'off')
    disp('Loading all datasets -----');
end


% --- Executes on selection change in listboxDataSets.
function listboxDataSets_Callback(hObject, eventdata, handles)
% hObject    handle to listboxDataSets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxDataSets contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxDataSets


% --- Executes during object creation, after setting all properties.
function listboxDataSets_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxDataSets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonLoadData.
function pushbuttonLoadData_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    disp(handles.listboxDataSets.Value)



% --- Executes on selection change in popupmenuAnalysisTypes.
function popupmenuAnalysisTypes_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuAnalysisTypes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuAnalysisTypes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuAnalysisTypes
    switch hObject.Value
        case 1
            Raw_Activity();
        case 3
            disp('3');
        otherwise
            disp('no');
    end

% --- Executes during object creation, after setting all properties.
function popupmenuAnalysisTypes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuAnalysisTypes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonAddData.
function pushbuttonAddData_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAddData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonRemoveData.
function pushbuttonRemoveData_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRemoveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
