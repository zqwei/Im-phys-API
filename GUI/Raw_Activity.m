function varargout = Raw_Activity(varargin)
% RAW_ACTIVITY MATLAB code for Raw_Activity.fig
%      RAW_ACTIVITY, by itself, creates a new RAW_ACTIVITY or raises the existing
%      singleton*.
%
%      H = RAW_ACTIVITY returns the handle to a new RAW_ACTIVITY or the handle to
%      the existing singleton*.
%
%      RAW_ACTIVITY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RAW_ACTIVITY.M with the given input arguments.
%
%      RAW_ACTIVITY('Property','Value',...) creates a new RAW_ACTIVITY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Raw_Activity_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Raw_Activity_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Raw_Activity

% Last Modified by GUIDE v2.5 14-Jan-2015 16:51:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Raw_Activity_OpeningFcn, ...
                   'gui_OutputFcn',  @Raw_Activity_OutputFcn, ...
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


% --- Executes just before Raw_Activity is made visible.
function Raw_Activity_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Raw_Activity (see VARARGIN)

% Choose default command line output for Raw_Activity
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Raw_Activity wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Raw_Activity_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenuListDataSet.
function popupmenuListDataSet_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuListDataSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuListDataSet contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuListDataSet


% --- Executes during object creation, after setting all properties.
function popupmenuListDataSet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuListDataSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxExampleNeuons.
function checkboxExampleNeuons_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxExampleNeuons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxExampleNeuons
     if get(hObject,'Value') == true
         handles.editNeuronSeq.Enable = 'on';
         handles.radiobuttonRandom.Enable = 'on';
     elseif get(hObject,'Value') == false
         handles.editNeuronSeq.Enable = 'off';
         handles.radiobuttonRandom.Enable = 'off';
     end


function editNeuronSeq_Callback(hObject, eventdata, handles)
% hObject    handle to editNeuronSeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNeuronSeq as text
%        str2double(get(hObject,'String')) returns contents of editNeuronSeq as a double
    NeuroSeq = str2num(get(hObject,'String')); %#ok<NASGU,ST2NM>
    disp(NeuroSeq)

% --- Executes during object creation, after setting all properties.
function editNeuronSeq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNeuronSeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function editRandomSize_Callback(hObject, eventdata, handles)
% hObject    handle to editRandomSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editRandomSize as text
%        str2double(get(hObject,'String')) returns contents of editRandomSize as a double
    SampleSize = str2double(get(hObject,'String'));
    NeuroSeq   = randn(1, SampleSize);
    handles.editNeuronSeq.String = num2str(NeuroSeq);
    editNeuronSeq_Callback(handles.editNeuronSeq, eventdata, handles)
    


% --- Executes during object creation, after setting all properties.
function editRandomSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editRandomSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobuttonRandom.
function radiobuttonRandom_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonRandom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonRandom
    if get(hObject,'Value')
        handles.editRandomSize.Enable = 'on';
    else
        handles.editRandomSize.Enable = 'off';
    end
