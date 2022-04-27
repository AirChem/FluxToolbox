function varargout = FluxInspector(varargin)
%FLUXINSPECTOR M-file for fluxinspector.fig
% This GUI is designed to examine airborne flux results.
% Included plots:
% 1) Time Series (selectable)
% 2) vertical profile
% 3) error summary
% 4) Lags (selectable)
% 5) Cospectra (selectable)
% 6) Spectra (selectable)
% 7) Map (selectable)
%
% INPUTS (all required)
% F: a big-ass structure contain all the stuff you need.
% name: name of flux species to plot.
%
% OUTPUT (optional) is a structure of handles to the various GUI objects.
%
% 20170701 GMW
% 20170710 GMW  Changed time series errors to std of mean
%               Fixed error on window close
% 20170801 GMW  Modified to handle new REwave flux error.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FluxInspector_OpeningFcn, ...
                   'gui_OutputFcn',  @FluxInspector_OutputFcn, ...
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

% --- Executes just before FluxInspector is made visible.
function FluxInspector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for FluxInspector
handles.output = hObject;

% updateplots handles structure
guidata(hObject, handles);

% GET INPUT VARIABLES
F  = varargin{1};
n  = varargin{2};
nL = sum(strncmp(fieldnames(F),'L',1));
c  = distColors(nL); %line colors

% if ~isempty(F.Div.divLegs)
%     divLegs = cell2mat(F.Div.divLegs(:,1)');
% else
%     divLegs=[];
% end

set(hObject,'name',[n ' - FluxInspector'])

%% SCALAR TIME SERIES
hTX = nan(nL,1); %handles
hTW = nan(nL,1);
set(gcf,'CurrentAxes',handles.ScalarTimeSeries)
hold on
for j=1:nL
    Ln = ['L' num2str(j)];
    t = F.(Ln).(n).wave.data.t;
    w = F.(Ln).(n).wave.data.w; w = (w - mean(w,'omitnan'))./std(w,'omitnan');
    x = F.(Ln).(n).wave.data.x; x = (x - mean(x,'omitnan'))./std(x,'omitnan');
    hTX(j) = plot(t,x,'k-');
    hTW(j) = plot(t,w,'m-');
end
ylabel('(x - <x>)/\sigma_x')
text(0.01,0.98,'w''','color','m')
text(0.08,0.98,'x''','color','k')
set(gca,'XTickLabel',[])

%% FLUX TIME SERIES
hTF = nan(nL,1); %time series handles
hTE = hTF; %errors
hTT = hTF; %text
hTC = hTF; %coi
hTD = hTF; %divergence points
set(gcf,'CurrentAxes',handles.FluxTimeSeries)
hold on
Ferr = sqrt((F.Avg1s.(n).*F.Avg1s.([n '_SEtot'])).^2 + F.Avg1s.([n '_REwave']).^2);
for j=1:nL
    k = F.Avg1s.Leg==j;
    h = shadedErrorBar(F.Avg1s.Time(k),F.Avg1s.(n)(k),Ferr(k),{'color',c(j,:)});
    delete(h.edge)
    hTF(j) = h.mainLine;
    hTE(j) = h.patch;
    hTT(j) = text(mean(F.Avg1s.Time(k),'omitnan'),min(F.Avg1s.(n)(k)),['L' num2str(j)],'color',c(j,:),'Units','data');
    
    %high coi points
    k2 = k & abs(F.Avg1s.([n '_qcoi']))>0.5;
    if sum(k2)
        hTC(j) = plot(F.Avg1s.Time(k2),F.Avg1s.(n)(k2),'kx'); 
    else
        hTC(j) = hTF(j); %dummy handle
    end
    
    %used points for divergence
    k2 = k & F.Div.([n '_usedpts']);
    if sum(k2)
        hTD(j) = plot(F.Avg1s.Time(k2),F.Avg1s.(n)(k2),'o','color',c(j,:),'MarkerSize',6);
    else
        hTD(j) = hTF(j);
    end
end
xlabel('Time')
ylabel('Flux')
text(0.01,0.15,'x: qcoi>0.5, o: div','color','k','Units','normalized')
set(gca,'xticklabelmode','auto')

linkaxes([handles.FluxTimeSeries handles.ScalarTimeSeries],'x')

%% MAP
set(gcf,'CurrentAxes',handles.Map)
hold on
mlat = IBWread('USEast_outline_lat.ibw'); mlat=mlat.y;
mlon = IBWread('USEast_outline_lon.ibw'); mlon=mlon.y;
plot(mlon,mlat,'-','color',[0.6 0.6 0.6]) %boundaries
plot(F.Avg1s.Lon,F.Avg1s.Lat,'k-') %whole flight track
hMA = nan(nL,1); %map handles
hMD = hMA; %divergence points
for j=1:nL %plot each leg
    k = F.Avg1s.Leg==j;
    hMA(j) = plot(handles.Map,F.Avg1s.Lon(k),F.Avg1s.Lat(k),'-','color',c(j,:),'lineWidth',3);
    
    %used points for divergence
    k2 = k & F.Div.([n '_usedpts']);
    if sum(k2)
        hMD(j) = plot(handles.Map,F.Avg1s.Lon(k2),F.Avg1s.Lat(k2),'o','color',c(j,:),'MarkerSize',6);
    else
        hMD(j) = hMA(j);
    end
end
set(handles.Map,'xlim',[min(F.Avg1s.Lon) max(F.Avg1s.Lon)],...
    'ylim',[min(F.Avg1s.Lat) max(F.Avg1s.Lat)],...
    'XTick',[],'YTick',[])
text(0.8,0.1,'o: div','color','k','Units','normalized','backgroundcolor','w')

%% VERTICAL PROFILE
set(gcf,'CurrentAxes',handles.VertProf);
hold on
hVP = nan(nL,1);
hVT = hVP;
hVE = hVP;
Ferr = sqrt((F.AvgLeg.(n) .*F.AvgLeg.([n '_SEtot'])).^2 + F.AvgLeg.([n '_REwave']).^2);
for j=1:nL
    h = ploterr(F.AvgLeg.(n)(j),F.AvgLeg.AltG(j),Ferr(j),[],'+');
    set(h,'color',c(j,:))
    hVT(j) = text(F.AvgLeg.(n)(j),F.AvgLeg.AltG(j)*1.01,['L' num2str(j)],'color',c(j,:),'Units','data');
    hVP(j) = h(1);
    hVE(j) = h(2); %error bar
end
xlabel('Flux')
ylabel('Alt (m)')
text(0.01,0.1,'o: divLeg','units','normalized')

divLegs = F.Div.divLegs;
hVD = hVP;
if ~isempty(divLegs)
    
    % add divergence subset
    u = F.Div.([n '_usedpts']);
    flux_div = BinAvg(F.Avg1s.Leg(u),F.Avg1s.(n)(u),1:nL);
    for i=1:nL
        if isnan(flux_div(i)), continue; end
        hVD(i) = plot(flux_div(i),F.AvgLeg.AltG(i),'o','lineWidth',5,'color',c(i,:));
    end
    
    % add fits
    fit = F.Div.([n '_fit']);
    z = [0 max(ylim)];
    for i=1:size(fit,1)
        plot(fit(i,1).*z + fit(i,2),z,'-','linewidth',3,'color',c(F.Div.divLegs{i,1}(1),:))
    end
end

%% ERRORS
set(gcf,'CurrentAxes',handles.Errors);
hold on
box on
plot(F.AvgLeg.Leg,F.AvgLeg.([n '_REfs01'])./abs(F.AvgLeg.(n)),'+','color',[1 0.6 0])
plot(F.AvgLeg.Leg,F.AvgLeg.([n '_REwave'])./abs(F.AvgLeg.(n)),'x','color',[0.5 0.5 0.5])
plot(F.AvgLeg.Leg,F.AvgLeg.([n '_REnoise'])./abs(F.AvgLeg.(n)),'b^')
plot(F.AvgLeg.Leg,F.AvgLeg.REturb,'cv')
plot(F.AvgLeg.Leg-0.1,F.AvgLeg.([n '_SErt']),'mo')
plot(F.AvgLeg.Leg-0.1,F.AvgLeg.SEturb,'rs')
plot(F.AvgLeg.Leg+0.1,sqrt((F.AvgLeg.([n '_REwave'])./abs(F.AvgLeg.(n))).^2 + F.AvgLeg.([n '_SEtot']).^2),'kh') %total error
errorbar(F.AvgLeg.Leg+0.1,F.AvgLeg.([n '_div'])-1,F.AvgLeg.([n '_div_err']),'gp')
plot(F.AvgLeg.Leg,F.L1.(n).data.SEacc.*ones(nL,1),'k-')
xlabel('Leg')
ylabel('Error/Flux')
h = legend('REfs01','REwav','REnoi','REtur','SErt','SEtur','Etot','DIV','Acc');
set(h,'fontsize',10)
set(gca,'ygrid','on')

%% LAGS
set(gcf,'CurrentAxes',handles.Lags);
hLG = nan(nL,1); %lags handles
hLO = nan(nL,1); %optimum lags handles
hold on
box on
optLag = fluxVarGrabber(F,{n,'options','xLag'},2);
for j=1:nL
    Ln = ['L' num2str(j)];
    hLG(j) = plot(F.(Ln).(n).lags,F.(Ln).(n).cov_wxdt,'-','color',c(j,:));
    hLO(j) = plot(optLag(j),F.AvgLeg.(n)(j),'p','color',c(j,:));
end
plot([0 0],get(gca,'ylim'),'k--')
plot(get(gca,'xlim'),[0 0],'k--')
xlabel('X Lag (# pts)')
ylabel('<w''x''>')
axis tight
xlim([-200 200])
grid on
set(gca,'ytick',[])
text(0.02,0.95,'line: EC','color','k','Units','normalized')
text(0.02,0.85,'star: WV','color','k','Units','normalized')
        
%% COSPECTRA
set(gcf,'CurrentAxes',handles.Cospectra);
hCO = nan(nL,1); %handles
hold on
box on
for j=1:nL
    Ln = ['L' num2str(j)];
    hCO(j) = semilogx(F.(Ln).(n).wave.freq,F.(Ln).(n).wave.co_nocoi,'-','color',c(j,:));
end
xlabel('Freq (Hz)')
ylabel('Wave Co(w,x)*f')
set(gca,'XScale','log','Ytick',[])
axis tight
set(gca,'xtick',[1e-4 1e-3 1e-2 0.1 1 10])
set(gca,'xticklabel',get(gca,'xtick'))
plot(get(gca,'xlim'),[0 0],'k--')

%% SPECTRA
set(gcf,'CurrentAxes',handles.Spectra);
hSX = nan(nL,1); %handles
hSW = nan(nL,1);
hold on
box on
for j=1:nL
    Ln = ['L' num2str(j)];
    hSW(j) = loglog(F.(Ln).(n).wave.freq,F.(Ln).(n).wave.w_psd./F.(Ln).(n).stats.w_r(3),'--','color',c(j,:));
    hSX(j) = loglog(F.(Ln).(n).wave.freq,F.(Ln).(n).wave.x_psd./F.(Ln).(n).stats.x_dtl(3),'-','color',c(j,:));
end
xlabel('Freq (Hz)')
ylabel('Wave PSD*f/\sigma^2')
set(gca,'XScale','log','YScale','log','YTick',[])
axis tight
set(gca,'xtick',[1e-4 1e-3 1e-2 0.1 1 10])
set(gca,'xticklabel',get(gca,'xtick'))
text(0.3,0.15,'- x'', -- w''','color','k','Units','normalized')

% fslope lines
f23 = fslopeline(xlim,max(ylim),-2/3,0);
loglog(xlim,f23*1,':','color',[0.6 0.6 0.6])
loglog(xlim,f23*10,':','color',[0.6 0.6 0.6])
loglog(xlim,f23*0.1,':','color',[0.6 0.6 0.6])
loglog(xlim,f23*0.3,':','color',[0.6 0.6 0.6])
loglog(xlim,f23*3,':','color',[0.6 0.6 0.6])
text(0.8,0.95,'-2/3','color',[0.6 0.6 0.6],'Units','normalized')

%% OTHER INITIALIZATION TASKS

% initialize buttons
for j = 1:30
    Ln = ['L' num2str(j)];
    if j<=nL
        set(handles.(Ln),'Value',get(handles.(Ln),'Max'))
    else
        delete(handles.(Ln))
    end
end
handles.plotme = true(nL,1);

%stuff into handles for later use
handles.F  = F; 
handles.n  = n;
handles.nL = nL;
handles.c  = c;
handles.hTX = hTX; %x time series data
handles.hTW = hTW; %w time series data
handles.hTF = hTF; %flux time series data
handles.hTT = hTT; %time series labels
handles.hTE = hTE; %time series errors
handles.hTC = hTC; %time series coi
handles.hTD = hTD; %time series divergence points
handles.hMA = hMA; %map
handles.hMD = hMD; %map divergence points
handles.hVP = hVP; %vertical profile data
handles.hVE = hVE; %vertical profile errors
handles.hVT = hVT; %vertical profile labels
handles.hVD = hVD; %vertical profile divergence subset
handles.hLG = hLG; %lags
handles.hLO = hLO; %lags optima
handles.hCO = hCO; %cospectra
handles.hSX = hSX; %spectra X
handles.hSW = hSW; %spectra W
guidata(hObject,handles);


%% INTERACTION FUNCTIONS

% --- Executes on button press in UPDATEPLOTS
function UPDATEPLOTS_Callback(hObject, eventdata, handles)
% change visibility of plots
hALL = [handles.hTX, handles.hTW, handles.hTF, handles.hTE, handles.hTT, handles.hTC, handles.hMA, ...
    handles.hVP, handles.hVE, handles.hVT, handles.hVD, ...
    handles.hLG, handles.hLO, handles.hCO, handles.hSX, handles.hSW, handles.hTD, handles.hMD];
h2plot = hALL(handles.plotme,:);
h2hide = hALL(~handles.plotme,:);
set(h2plot(:),'Visible','on')
set(h2hide(:),'Visible','off')
set(handles.FluxTimeSeries,'xlimmode','auto')

% --- Executes on button press in ALLNONE.
function ALLNONE_Callback(hObject, eventdata, handles)
%toggle between all and none button states
if all(handles.plotme)
    handles.plotme(:) = 0;
else
    handles.plotme(:) = 1;
end

% update buttons
for i=1:handles.nL
    h = handles.(['L' num2str(i)]);
    set(h,'Value',handles.plotme(i))
end
guidata(hObject,handles)

% --- Update plotme flag on button presses
function L1_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(1)=1; else handles.plotme(1)=0; end; guidata(hObject,handles)
function L2_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(2)=1; else handles.plotme(2)=0; end; guidata(hObject,handles)
function L3_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(3)=1; else handles.plotme(3)=0; end; guidata(hObject,handles)
function L4_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(4)=1; else handles.plotme(4)=0; end; guidata(hObject,handles)
function L5_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(5)=1; else handles.plotme(5)=0; end; guidata(hObject,handles)
function L6_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(6)=1; else handles.plotme(6)=0; end; guidata(hObject,handles)
function L7_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(7)=1; else handles.plotme(7)=0; end; guidata(hObject,handles)
function L8_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(8)=1; else handles.plotme(8)=0; end; guidata(hObject,handles)
function L9_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(9)=1; else handles.plotme(9)=0; end; guidata(hObject,handles)
function L10_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(10)=1; else handles.plotme(10)=0; end; guidata(hObject,handles)
function L11_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(11)=1; else handles.plotme(11)=0; end; guidata(hObject,handles)
function L12_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(12)=1; else handles.plotme(12)=0; end; guidata(hObject,handles)
function L13_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(13)=1; else handles.plotme(13)=0; end; guidata(hObject,handles)
function L14_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(14)=1; else handles.plotme(14)=0; end; guidata(hObject,handles)
function L15_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(15)=1; else handles.plotme(15)=0; end; guidata(hObject,handles)
function L16_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(16)=1; else handles.plotme(16)=0; end; guidata(hObject,handles)
function L17_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(17)=1; else handles.plotme(17)=0; end; guidata(hObject,handles)
function L18_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(18)=1; else handles.plotme(18)=0; end; guidata(hObject,handles)
function L19_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(19)=1; else handles.plotme(19)=0; end; guidata(hObject,handles)
function L20_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(20)=1; else handles.plotme(20)=0; end; guidata(hObject,handles)
function L21_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(21)=1; else handles.plotme(21)=0; end; guidata(hObject,handles)
function L22_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(22)=1; else handles.plotme(22)=0; end; guidata(hObject,handles)
function L23_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(23)=1; else handles.plotme(23)=0; end; guidata(hObject,handles)
function L24_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(24)=1; else handles.plotme(24)=0; end; guidata(hObject,handles)
function L25_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(25)=1; else handles.plotme(25)=0; end; guidata(hObject,handles)
function L26_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(26)=1; else handles.plotme(26)=0; end; guidata(hObject,handles)
function L27_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(27)=1; else handles.plotme(27)=0; end; guidata(hObject,handles)
function L28_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(28)=1; else handles.plotme(28)=0; end; guidata(hObject,handles)
function L29_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(29)=1; else handles.plotme(29)=0; end; guidata(hObject,handles)
function L30_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==get(hObject,'Max'), handles.plotme(30)=1; else handles.plotme(30)=0; end; guidata(hObject,handles)

% --- Outputs from this function are returned to the command line.
function varargout = FluxInspector_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles;

% delete fns (need this to avoid error...don't ask why)
function L1_DeleteFcn(hObject, eventdata, handles)

% That's all folks
