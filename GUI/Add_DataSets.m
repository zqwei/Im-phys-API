function varargout = Add_DataSets(varargin)
% ADD_DATASETS MATLAB code for Add_DataSets.fig
%      ADD_DATASETS, by itself, creates a new ADD_DATASETS or raises the existing
%      singleton*.
%
%      H = ADD_DATASETS returns the handle to a new ADD_DATASETS or the handle to
%      the existing singleton*.
%
%      ADD_DATASETS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADD_DATASETS.M with the given input arguments.
%
%      ADD_DATASETS('Property','Value',...) creates a new ADD_DATASETS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Add_DataSets_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Add_DataSets_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Add_DataSets

% Last Modified by GUIDE v2.5 16-Jan-2015 17:20:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Add_DataSets_OpeningFcn, ...
                   'gui_OutputFcn',  @Add_DataSets_OutputFcn, ...
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


% --- Executes just before Add_DataSets is made visible.
function Add_DataSets_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Add_DataSets (see VARARGIN)

% Choose default command line output for Add_DataSets
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Add_DataSets wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Add_DataSets_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editName_Callback(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editName as text
%        str2double(get(hObject,'String')) returns contents of editName as a double


% --- Executes during object creation, after setting all properties.
function editName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editDir_Callback(hObject, eventdata, handles)
% hObject    handle to editDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDir as text
%        str2double(get(hObject,'String')) returns contents of editDir as a double
    pushbuttonDir_Callback(handles.pushbuttonDir, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function editDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonDir.
function pushbuttonDir_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    dialog_title            = handles.popupmenuDataType.String{handles.popupmenuDataType.Value};
    start_path              = handles.editDir.String;
    folder_name             = uigetdir(start_path,dialog_title);
    handles.editDir.String  = folder_name;

% --- Executes on button press in radiobuttonIsSpike.
function radiobuttonIsSpike_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonIsSpike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonIsSpike


% --- Executes on button press in radiobuttonIsCa.
function radiobuttonIsCa_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonIsCa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonIsCa


% --- Executes on button press in pushbuttonAdd.
function pushbuttonAdd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    if ~exist(handles.editDir.String, 'dir')
        msgbox('Choose directory location first please.')
        return;
    end

    popupmenuDataType_Callback(handles.popupmenuDataType, eventdata, handles);
    params.frameRate       =  str2double(handles.editFrameRate.String);
    params.binsize         =  str2double(handles.editBinSize.String)/1000;
    params.polein          =  str2double(handles.editPoleIn.String);
    params.poleout         =  str2double(handles.editPoleOut.String);
    minTimeToAnalysis      =  round(str2double(handles.editMinTime.String) * params.frameRate);
    maxTimeToAnalysis      =  round(str2double(handles.editMaxTime.String) * params.frameRate);
    params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
    params.timeSeries      = params.timeWindowIndexRange * params.binsize;
    params.minNumTrialToAnalysis =  str2double(handles.editMinTrial.String);
    params.Fm              =  str2double(handles.editFm.String);
    params.Kd              =  str2double(handles.editKd.String);
    params.n               =  str2double(handles.editn.String);
    params.tau_decay       =  str2double(handles.editTd.String)/1000;
    params.tau_rise        =  str2double(handles.editTr.String)/1000;
    
    if exist('TempDat/DataList.mat','file')
        load('TempDat/DataList.mat')
    else
        DataSetList        = [];
    end
    
    nList                  = length(DataSetList);

    % DataSetList{1}     = handles.editName.String;
    
    switch handles.popupmenuDataType.Value
        case 1 % Nuo Li's data
            nDataDir                       = [handles.editDir.String '/'];
            nFileList                      = dir([nDataDir '*.mat']);
            nDataSet                       = getSpikeData(nDataDir, nFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize); %#ok<NASGU>
            handles.editName.String        = ['Nuo_Li_' datestr(datetime('now'), 'yyyy_mm_dd_HH_MM')];
            DataSetList(nList + 1).name    = ['Spikes_' handles.editName.String];
            DataSetList(nList + 1).params  = params; 
            save(['TempDat/Spikes_' handles.editName.String '.mat'], 'nDataSet');
            nList                          = nList + 1;
            nDataSet                       = getFakeCaImagingData(nDataDir, nFileList, params.minNumTrialToAnalysis, params); %#ok<NASGU>
            DataSetList(nList + 1).name    = ['Fake_Ca_' handles.editName.String];
            DataSetList(nList + 1).params  = params; 
            save(['TempDat/Fake_Ca_' handles.editName.String '.mat'], 'nDataSet');            
        case 2 % Tsai-wen Chen's data
            nDataDir                       = [handles.editDir.String '/'];
            nFileList                      = dir([nDataDir '*.mat']);
            nDataSet                       = getCaImagingData(nDataDir, nFileList, params.minNumTrialToAnalysis, params); %#ok<NASGU>
            handles.editName.String        = ['Tsaiwen_Chen_' datestr(datetime('now'), 'yyyy_mm_dd_HH_MM')];
            DataSetList(nList + 1).name    = ['Ca_' handles.editName.String];
            DataSetList(nList + 1).params  = params; 
            save(['TempDat/Ca_' handles.editName.String '.mat'], 'nDataSet');            
        case 3 % Yinan Wan's data
            nDataDir                       = [handles.editDir.String '/'];
            nFileList                      = dir([nDataDir '*.mat']);
            nDataSet                       = [];
            handles.editName.String        = ['Yinan_Wan_' datestr(datetime('now'), 'yyyy_mm_dd_HH_MM')];
            DataSetList(nList + 1).name    = ['Ca_' handles.editName.String];
            DataSetList(nList + 1).params  = params; 
            save(['TempDat/Ca_' handles.editName.String '.mat'], 'nDataSet');
    end
    
    save('TempDat/DataList.mat', 'DataSetList');
    
    
    

function editFrameRate_Callback(hObject, eventdata, handles)
% hObject    handle to editFrameRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFrameRate as text
%        str2double(get(hObject,'String')) returns contents of editFrameRate as a double
    handles.editBinSize.String    = num2str(1000/str2double(hObject.String));

% --- Executes during object creation, after setting all properties.
function editFrameRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFrameRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editBinSize_Callback(hObject, eventdata, handles)
% hObject    handle to editBinSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBinSize as text
%        str2double(get(hObject,'String')) returns contents of editBinSize as a double
    handles.editFrameRate.String    = num2str(1000/str2double(hObject.String));

% --- Executes during object creation, after setting all properties.
function editBinSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBinSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPoleIn_Callback(hObject, eventdata, handles)
% hObject    handle to editPoleIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPoleIn as text
%        str2double(get(hObject,'String')) returns contents of editPoleIn as a double


% --- Executes during object creation, after setting all properties.
function editPoleIn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPoleIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPoleOut_Callback(hObject, eventdata, handles)
% hObject    handle to editPoleOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPoleOut as text
%        str2double(get(hObject,'String')) returns contents of editPoleOut as a double


% --- Executes during object creation, after setting all properties.
function editPoleOut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPoleOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMinTime_Callback(hObject, eventdata, handles)
% hObject    handle to editMinTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinTime as text
%        str2double(get(hObject,'String')) returns contents of editMinTime as a double


% --- Executes during object creation, after setting all properties.
function editMinTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMaxTime_Callback(hObject, eventdata, handles)
% hObject    handle to editMaxTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMaxTime as text
%        str2double(get(hObject,'String')) returns contents of editMaxTime as a double


% --- Executes during object creation, after setting all properties.
function editMaxTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFm_Callback(hObject, eventdata, handles)
% hObject    handle to editFm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFm as text
%        str2double(get(hObject,'String')) returns contents of editFm as a double


% --- Executes during object creation, after setting all properties.
function editFm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editKd_Callback(hObject, eventdata, handles)
% hObject    handle to editKd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editKd as text
%        str2double(get(hObject,'String')) returns contents of editKd as a double


% --- Executes during object creation, after setting all properties.
function editKd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editKd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuDataType.
function popupmenuDataType_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuDataType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuDataType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuDataType
    switch get(hObject,'Value')
        case 1 % Nuo Li (Spikes)
            handles.editFm.Enable  = 'on';
            handles.editKd.Enable  = 'on';
            handles.editn.Enable   = 'on';
            handles.editTd.Enable  = 'on';
            handles.editTr.Enable  = 'on';
        otherwise
            handles.editFm.Enable  = 'off';
            handles.editKd.Enable  = 'off';
            handles.editn.Enable   = 'off';
            handles.editTd.Enable  = 'off';
            handles.editTr.Enable  = 'off';
    end
    
    % handles.editName.String        = [hObject.String{hObject.Value} ': ' datestr(datetime('now'), 'yyyy_mm_dd_HH_MM')];
            


% --- Executes during object creation, after setting all properties.
function popupmenuDataType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuDataType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editn_Callback(hObject, eventdata, handles)
% hObject    handle to editn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editn as text
%        str2double(get(hObject,'String')) returns contents of editn as a double


% --- Executes during object creation, after setting all properties.
function editn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editTr_Callback(hObject, eventdata, handles)
% hObject    handle to editTr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTr as text
%        str2double(get(hObject,'String')) returns contents of editTr as a double


% --- Executes during object creation, after setting all properties.
function editTr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function editTd_Callback(hObject, eventdata, handles)
% hObject    handle to editTd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTd as text
%        str2double(get(hObject,'String')) returns contents of editTd as a double


% --- Executes during object creation, after setting all properties.
function editTd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMinTrial_Callback(hObject, eventdata, handles)
% hObject    handle to editMinTrial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinTrial as text
%        str2double(get(hObject,'String')) returns contents of editMinTrial as a double


% --- Executes during object creation, after setting all properties.
function editMinTrial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinTrial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
