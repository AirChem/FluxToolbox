function varargout = FluxLegPicker(varargin)
%FLUXLEGPICKER M-file for fluxLegPicker.fig
% This GUI is designed to help you pick out individual flight legs for fluxification.
%      FLUXLEGPICKER, by itself, creates a new FLUXLEGPICKER or raises the existing
%      singleton*.
%
%      H = FLUXLEGPICKER returns the handle to a new FLUXLEGPICKER or the handle to
%      the existing singleton*.
%
%      FLUXLEGPICKER('Property','Value',...) creates a new FLUXLEGPICKER using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to FluxLegPicker_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      FLUXLEGPICKER('CALLBACK') and FLUXLEGPICKER('CALLBACK',hObject,...) call the
%      local function named CALLBACK in FLUXLEGPICKER.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% INPUTS (all required)
% data: a structure containing the following fields:
%   time:   UTC seconds
%   alt:    altitude, m
%   lat:    latitude
%   lon:    longitude
%   hdg:    aircraft heading
%   roll:   aircraft roll angle
%   pitch:  aircraft pitch angle
%   w:      vertical wind speed, m/s
% legTimes: a 2-column matrix of start and stop times for individual legs.

% Edit the above text to modify the response to help FluxLegPicker

% Last Modified by GUIDE v2.5 10-May-2017 08:33:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FluxLegPicker_OpeningFcn, ...
                   'gui_OutputFcn',  @FluxLegPicker_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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

% --- Executes just before FluxLegPicker is made visible.
function FluxLegPicker_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for FluxLegPicker
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FluxLegPicker wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% break out varargin
handles.d        = varargin{1};
handles.legTimes = varargin{2};
guidata(hObject,handles);

% INITIALIZE PLOTS
hold(handles.axes1,'on');
hold(handles.axes2,'on');
hold(handles.axes3,'on');

% time series
plot(handles.axes1,handles.d.time,handles.d.alt,'c-')
plot(handles.axes1,handles.d.time,handles.d.hdg,'-','color',[1 0.6 0])
set(gcf,'CurrentAxes',handles.axes1)
text(0.05,0.95,'Alt','color','c','Units','normalized')
text(0.10,0.95,'Hdg','color',[1 0.6 0],'Units','normalized')
% ylabel(handles.axes1,'Alt (m)')
set(handles.axes1,'xticklabel',[])

plot(handles.axes2,handles.d.time,handles.d.w,'c-')
plot(handles.axes2,handles.d.time,handles.d.pitch,'-','color',[1 0.6 0])
plot(handles.axes2,handles.d.time,handles.d.roll,'g-')
% ylabel(handles.axes2,'Hdg (deg)')
set(gcf,'CurrentAxes',handles.axes2)
text(0.05,0.95,'W','color','c','Units','normalized')
text(0.10,0.95,'Pitch','color',[1 0.6 0],'Units','normalized')
text(0.20,0.95,'Roll','color','g','Units','normalized')
xlabel(handles.axes2,'UTC (s)')

set([handles.axes1 handles.axes2],'XGrid','on','YGrid','on')
linkaxes([handles.axes1 handles.axes2],'x')

%map
mlat = IBWread('USEast_outline_lat.ibw'); mlat=mlat.y;
mlon = IBWread('USEast_outline_lon.ibw'); mlon=mlon.y;
plot(handles.axes3,mlon,mlat,'-','color',[0.6 0.6 0.6])
plot(handles.axes3,handles.d.lon,handles.d.lat,'k-')
set(handles.axes3,'xlim',[min(handles.d.lon) max(handles.d.lon)],...
    'ylim',[min(handles.d.lat) max(handles.d.lat)])
set(handles.axes3,'XGrid','on','YGrid','on')

% plot all legs
N = size(handles.legTimes,1);
set(handles.edit1,'String',['1:' num2str(N)]);
i=[];
for j=1:N
    i=[i;find(handles.d.time>=handles.legTimes(j,1) & handles.d.time<=handles.legTimes(j,2))];
end
handles.p1 = plot(handles.axes1,handles.d.time(i),handles.d.alt(i),'b.');
handles.p2 = plot(handles.axes1,handles.d.time(i),handles.d.hdg(i),'r.');
handles.p3 = plot(handles.axes2,handles.d.time(i),handles.d.w(i),'b.');
handles.p4 = plot(handles.axes2,handles.d.time(i),handles.d.pitch(i),'r.');
handles.p5 = plot(handles.axes2,handles.d.time(i),handles.d.roll(i),'.','color',[0 0.5 0]);
handles.p6 = plot(handles.axes3,handles.d.lon(i),handles.d.lat(i),'r.');
guidata(hObject,handles)

% populate table
tableStart = num2cell(zeros(20,2));
tableStart(1:N,:) = num2cell(handles.legTimes);
set(handles.uitable1,'Data',tableStart)


% --- Outputs from this function are returned to the command line.
function varargout = FluxLegPicker_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1 (PLOT LEGS).
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

legs = str2num(get(handles.edit1,'String'));
i=[];
for j=legs
    legTime = handles.legTimes(j,:);
    i=[i;find(handles.d.time>=legTime(1) & handles.d.time<=legTime(2))];
end

set(handles.p1,'XData',handles.d.time(i),'YData',handles.d.alt(i))
set(handles.p2,'XData',handles.d.time(i),'YData',handles.d.hdg(i))
set(handles.p3,'XData',handles.d.time(i),'YData',handles.d.w(i))
set(handles.p4,'XData',handles.d.time(i),'YData',handles.d.pitch(i))
set(handles.p5,'XData',handles.d.time(i),'YData',handles.d.roll(i))
set(handles.p6,'XData',handles.d.lon(i),'YData',handles.d.lat(i))


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

% check if ok
legs = str2num(get(hObject,'String'));
oklegs = 1:size(handles.legTimes,1);
if any(~ismember(legs,oklegs))
    set(hObject,'String',['1:' num2str(size(handles.legTimes,1))])
end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

% update index
handles.legTimes(eventdata.Indices(1),eventdata.Indices(2)) = eventdata.NewData;
guidata(hObject,handles)


% --- Executes on button press in pushbutton2 (DUMP TO CW).
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% dump indices to command window for easy cut 'n paste
outtie = sprintf('%1.2f, %1.2f;\n',handles.legTimes');
disp('New Leg Times:')
disp(outtie)
% assignin('base','NewLegTimesString',outtie);
