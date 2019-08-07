function varargout = GUI_PLOT_DLSOCT(varargin)
% GUI_PLOT_DLSOCT MATLAB code for GUI_PLOT_DLSOCT.fig
%      GUI_PLOT_DLSOCT, by itself, creates a new GUI_PLOT_DLSOCT or raises the existing
%      singleton*.
%
%      H = GUI_PLOT_DLSOCT returns the handle to a new GUI_PLOT_DLSOCT or the handle to
%      the existing singleton*.
%
%      GUI_PLOT_DLSOCT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_PLOT_DLSOCT.M with the given input arguments.
%
%      GUI_PLOT_DLSOCT('Property','Value',...) creates a new GUI_PLOT_DLSOCT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_PLOT_DLSOCT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_PLOT_DLSOCT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_PLOT_DLSOCT

% Last Modified by GUIDE v2.5 26-Jul-2019 11:37:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_PLOT_DLSOCT_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_PLOT_DLSOCT_OutputFcn, ...
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


% --- Executes just before GUI_PLOT_DLSOCT is made visible.
function GUI_PLOT_DLSOCT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_PLOT_DLSOCT (see VARARGIN)
%% add MATLAB functions' path
handles.CodePath=pwd;
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
handles.defpath='H:\BU';
handles.startZ=1;
handles.stackZ=100;
handles.ShowSide=0;
handles.cRange=2;
% Choose default command line output for GUI_PLOT_DLSOCT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_PLOT_DLSOCT wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_PLOT_DLSOCT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_LoadData.
function btn_LoadData_Callback(hObject, eventdata, handles)
% hObject    handle to btn_LoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
%% select file path %%%%%%%%%
defaultpath=handles.defpath;
[filename,datapath]=uigetfile(defaultpath);
handles.defpath=datapath;
handles.filename=filename;
guidata(hObject, handles);
%% load data %%%%%%%%%%
%%%%% input number of sub GG to be loaded %%%%%%%%
lding=msgbox(['Loading data...  ',datestr(now,'DD:HH:MM')]);
% load([datapath,'V.mat']);
% handles.V=imgaussfilt3(V,0.5);
% load([datapath,'Vz.mat']);
% handles.Vz=imgaussfilt3(Vz,0.5);
% load([datapath,'Vt.mat']);
% handles.Vt=imgaussfilt3(Vt,0.5);
% load([datapath,'D.mat']);
% handles.D=imgaussfilt3(D,0.5);
% load([datapath,'R.mat']);
% handles.R=imgaussfilt3(R,0.5);
% load([datapath,'Mf.mat']);
% handles.Mf=imgaussfilt3(Mf,0.5);
% load([datapath,'Ms.mat']);
% handles.Ms=imgaussfilt3(Ms,0.5);
load([datapath,filename]);
handles.Vz=imgaussfilt3(Vz,0.5);
handles.Vt=imgaussfilt3(Vt,0.5);
handles.D=imgaussfilt3(D,0.5);
handles.R=imgaussfilt3(R,0.5);
handles.Mf=imgaussfilt3(Mf,0.5);
handles.Ms=imgaussfilt3(Ms,0.5);
handles.V=imgaussfilt3(V,0.5);
% temp 
% Vsign=sign(Vz);
% Vsign(Vsign~=-1)=1;
% V=sqrt(Vz.^2+Vt.^2).*Vsign; % Total velocity
V=sqrt(Vz.^2+Vt.^2);
handles.V=imgaussfilt3(V,0.5);
lded=msgbox(['Data loaded. ',datestr(now,'DD:HH:MM')]);
pause(1);
delete(lding); delete(lded);
[nz,nx,ny]=size(V);
guidata(hObject, handles);
%%%%%%%%%%%%%
prompt={'dX (um): ', 'dY(um): ', 'dZ size (um)'};
name='Enter Imaging info';
defaultvalue={'1.5','1.5','2.9'};
dXYZinput=inputdlg(prompt,name, 1, defaultvalue);
handles.Xcoor=[1:nx]*str2num(dXYZinput{1});
handles.Ycoor=[1:ny]*str2num(dXYZinput{2});
handles.Zcoor=[1:nz]*str2num(dXYZinput{3});

%% plot Vz enface MIP
[Vcmap, Vzcmap, Dcmap, Mfcmap, Rcmap, g1OCTAcmap]=Colormaps_DLSOCT;
vSign=sign(squeeze(mean(handles.Vz(:,:,:),1)));
vSign(vSign~=-1)=1;
axes(handles.axes1)
imagesc(squeeze(max(abs(handles.Vz(:,:,:)),[],1)).*vSign); 
colormap(Vzcmap); caxis([-2 2]); colorbar
axis equal tight
title('Vz [mm/s]')

handles.DataShow=handles.Vz;
handles.cMap=Vzcmap;
handles.Title='Vz [mm/s]';
handles.Caxis=[-handles.cRange,handles.cRange];
handles.DataSlt='Vz'; 
guidata(hObject, handles);

% --- Executes on selection change in POP_DataSlect.
function POP_DataSlect_Callback(hObject, eventdata, handles)
% hObject    handle to POP_DataSlect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns POP_DataSlect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from POP_DataSlect
[Vcmap, Vzcmap, Dcmap, Mfcmap, Rcmap, g1OCTAcmap]=Colormaps_DLSOCT;
contents = cellstr(get(hObject,'String'));
handles.DataSlt=contents{get(hObject,'Value')};
if (strcmp(handles.DataSlt,'V'))
    handles.DataShow=handles.V;
    handles.cMap=Vzcmap;
    handles.Title='V [mm/s]';
    handles.Caxis=[-handles.cRange,handles.cRange];
elseif (strcmp(handles.DataSlt,'Vt'))
    handles.DataShow=handles.Vt;
    handles.cMap=Vzcmap;
    handles.Title='Vt [mm/s]';
    handles.Caxis=[-handles.cRange,handles.cRange];
elseif (strcmp(handles.DataSlt,'D'))
    handles.DataShow=handles.D;
    handles.cMap=Dcmap;
    handles.Title='D [um^2/s]';
    handles.Caxis=[0,100];
elseif (strcmp(handles.DataSlt,'R'))
    handles.DataShow=handles.R;
    handles.cMap=Rcmap;
    handles.Title='R';
    handles.Caxis=[0, 1];
elseif (strcmp(handles.DataSlt,'Mf'))
    handles.DataShow=handles.Mf;
    handles.cMap=Mfcmap;
    handles.Title='Mf';
    handles.Caxis=[0, 1];
elseif (strcmp(handles.DataSlt,'Ms'))
    handles.DataShow=handles.Ms;
    handles.cMap=Rcmap;
    handles.Title='Ms';
    handles.Caxis=[0, 1];
else
    handles.DataShow=handles.Vz;
    handles.cMap=Vzcmap;
    handles.Title='Vz [mm/s]';
    handles.Caxis=[-handles.cRange,handles.cRange];
end
guidata(hObject, handles); 
if (strcmp(handles.DataSlt,'Ms'))
    axes(handles.axes1)
    imagesc(squeeze(min(abs(handles.DataShow(:,:,:)),[],1)));
    colormap(handles.cMap); caxis(handles.Caxis); colorbar
    axis equal tight
    title(handles.Title)
elseif (strcmp(handles.DataSlt,'R')) || (strcmp(handles.DataSlt,'Mf')) || (strcmp(handles.DataSlt,'D'))
    axes(handles.axes1)
    imagesc(squeeze(max(abs(handles.DataShow(:,:,:)),[],1)));
    colormap(handles.cMap); caxis(handles.Caxis); colorbar
    axis equal tight
    title(handles.Title)
elseif (strcmp(handles.DataSlt,'V')) || (strcmp(handles.DataSlt,'Vt')) || (strcmp(handles.DataSlt,'Vz')) || (strcmp(handles.DataSlt,'D'))
    vSign=sign(squeeze(mean(handles.Vz(:,:,:),1)));
    vSign(vSign~=-1)=1;
    R=handles.R;
    R(R>0)=1;
    R(R<1)=0;
    axes(handles.axes1)
    imagesc(squeeze(max(abs(handles.DataShow(:,:,:)).*R,[],1)).*vSign);
    colormap(handles.cMap); caxis(handles.Caxis); colorbar
    axis equal tight
    title(handles.Title)
end
guidata(hObject, handles); 

% --- Executes on button press in btn_plot.
function btn_plot_Callback(hObject, eventdata, handles)
% hObject    handle to btn_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
[nz,nx,ny]=size(handles.V);

prompt={'SideView(N:0;Y:1)','MIP zStart', 'MIP zStack','cRange'};
name='Enter Imaging info';
defaultvalue={num2str(handles.ShowSide),num2str(handles.startZ),num2str(min(nz,nz-handles.startZ)),num2str(handles.cRange)};
dXYZinput=inputdlg(prompt,name, 1, defaultvalue);
handles.ShowSide=str2num(dXYZinput{1});
zStart=str2num(dXYZinput{2});
zStack=str2num(dXYZinput{3});
handles.cRange=str2num(dXYZinput{4});

handles.startZ=zStart;
handles.endZ=min((zStart+zStack)-1,nz);
guidata(hObject, handles);  

%% PLOT
[Vcmap, Vzcmap, Dcmap, Mfcmap, Rcmap, g1OCTAcmap]=Colormaps_DLSOCT;
if handles.ShowSide==1 % plot SideView figures
    axes(handles.axes1);
    handles.slt(1)=str2num(get(handles.zStart,'string'));
    [handles.slt(3), handles.slt(2)]=ginput(1); % [y x]
    handles.slt=round(handles.slt);
    handles.Show=['3DSP-','[',num2str(handles.slt), ']']; % 3D single plane show
    
    handles.fig=figure;
    set(handles.fig,'Position',[400 100 1400 800])
    if (strcmp(handles.DataSlt,'Ms'))
        subplot(2,2,1) % XY single plan enface view
        Showxy=(squeeze(max(abs(handles.DataShow(handles.slt(1),:,:)),[],1)));
        Showxy(handles.slt(2)+[0 1],:)=max(Showxy(:));
        Showxy(:,handles.slt(3)+[0 1])=max(Showxy(:));
        imagesc(handles.Xcoor,handles.Ycoor,Showxy);
        colormap(handles.cMap); caxis(handles.Caxis); colorbar
        axis equal tight
        title([handles.Title, ', XY, iz=',num2str(handles.slt(1))])
        xlabel('X [um]'); ylabel('Y [um]')
        
        subplot(2,2,2) % xz side view
        Showxz=squeeze(handles.DataShow(:,handles.slt(2),:));
        Showxz(handles.slt(1),:)=max(Showxz(:));
        Showxz(:,handles.slt(3))=max(Showxz(:));
        imagesc(handles.Ycoor,handles.Zcoor,Showxz);
        colormap(handles.cMap); caxis(handles.Caxis); colorbar
        axis equal tight
        title([handles.Title, ', XZ, iY=',num2str(handles.slt(2))])
        xlabel('X [um]'); ylabel('Z [um]')
        
        subplot(2,2,3) % XY enface view MIP
        showMIP=(squeeze(min(abs(handles.DataShow(handles.startZ:handles.endZ,:,:)),[],1)));
        showMIP(handles.slt(2)+[0 1],:)=max(showMIP(:));
        showMIP(:,handles.slt(3)+[0 1])=max(showMIP(:));
        imagesc(handles.Xcoor,handles.Ycoor,showMIP);
        colormap(handles.cMap); caxis(handles.Caxis); colorbar
        axis equal tight
        title([handles.Title, ', V MIP [',num2str([handles.startZ handles.endZ]), '] [mm/s]'])
        xlabel('X [um]'); ylabel('Y [um]')
        
        subplot(2,2,4) %yz side view
        Showyz=squeeze(handles.DataShow(:,:,handles.slt(3)));
        Showyz(handles.slt(1),:)=max(Showyz(:));
        Showyz(:,handles.slt(2))=max(Showyz(:));
        imagesc(handles.Xcoor,handles.Zcoor,Showyz);
        colormap(handles.cMap); caxis(handles.Caxis); colorbar
        axis equal tight
        title([handles.Title, ', YZ, iX=',num2str(handles.slt(3))])
        xlabel('Y [um]'); ylabel('Z [um]')
    elseif (strcmp(handles.DataSlt,'R')) || (strcmp(handles.DataSlt,'Mf')) || (strcmp(handles.DataSlt,'D'))
        subplot(2,2,1) % XY single plan enface view
        Showxy=(squeeze(max(abs(handles.DataShow(handles.slt(1),:,:)),[],1)));
        Showxy(handles.slt(2)+[0 1],:)=max(Showxy(:));
        Showxy(:,handles.slt(3)+[0 1])=max(Showxy(:));
        imagesc(handles.Xcoor,handles.Ycoor,Showxy);
        colormap(handles.cMap); caxis(handles.Caxis); colorbar
        axis equal tight
        title([handles.Title, ', XY, iz=',num2str(handles.slt(1))])
        xlabel('X [um]'); ylabel('Y [um]')
        
        subplot(2,2,2) % xz side view
        Showxz=squeeze(handles.DataShow(:,handles.slt(2),:));
        Showxz(handles.slt(1),:)=max(Showxz(:));
        Showxz(:,handles.slt(3))=max(Showxz(:));
        imagesc(handles.Ycoor,handles.Zcoor,Showxz);
        colormap(handles.cMap); caxis(handles.Caxis); colorbar
        axis equal tight
        title([handles.Title, ', XZ, iY=',num2str(handles.slt(2))])
        xlabel('X [um]'); ylabel('Z [um]')
        
        subplot(2,2,3) % XY enface view MIP
        showMIP=(squeeze(max(abs(handles.DataShow(handles.startZ:handles.endZ,:,:)),[],1)));
        showMIP(handles.slt(2)+[0 1],:)=max(showMIP(:));
        showMIP(:,handles.slt(3)+[0 1])=max(showMIP(:));
        imagesc(handles.Xcoor,handles.Ycoor,showMIP);
        colormap(handles.cMap); caxis(handles.Caxis); colorbar
        axis equal tight
        title([handles.Title, ', V MIP [',num2str([handles.startZ handles.endZ]), '] [mm/s]'])
        xlabel('X [um]'); ylabel('Y [um]')
        
        subplot(2,2,4) %yz side view
        Showyz=squeeze(handles.DataShow(:,:,handles.slt(3)));
        Showyz(handles.slt(1),:)=max(Showyz(:));
        Showyz(:,handles.slt(2))=max(Showyz(:));
        imagesc(handles.Xcoor,handles.Zcoor,Showyz);
        colormap(handles.cMap); caxis(handles.Caxis); colorbar
        axis equal tight
        title([handles.Title, ', YZ, iX=',num2str(handles.slt(3))])
        xlabel('Y [um]'); ylabel('Z [um]')
    elseif (strcmp(handles.DataSlt,'V')) || (strcmp(handles.DataSlt,'Vt')) || (strcmp(handles.DataSlt,'Vz')) || (strcmp(handles.DataSlt,'D'))
        subplot(2,2,1) % XY single plan enface view
        Showxy=(squeeze(max(abs(handles.DataShow(handles.slt(1),:,:)),[],1)).*sign(squeeze(mean(handles.Vz(handles.slt(1),:,:),1))));
        Showxy(handles.slt(2)+[0 1],:)=max(Showxy(:));
        Showxy(:,handles.slt(3)+[0 1])=max(Showxy(:));
        imagesc(handles.Xcoor,handles.Ycoor,Showxy);
        colormap(handles.cMap); caxis([-handles.cRange handles.cRange]); colorbar
        axis equal tight
        title([handles.Title, ', XY, iz=',num2str(handles.slt(1))])
        xlabel('X [um]'); ylabel('Y [um]')
        
        subplot(2,2,2) % xz side view
        Showxz=squeeze(handles.DataShow(:,handles.slt(2),:)).*sign(squeeze(handles.Vz(:,handles.slt(2),:)));
        Showxz(handles.slt(1),:)=max(Showxz(:));
        Showxz(:,handles.slt(3))=max(Showxz(:));
        imagesc(handles.Ycoor,handles.Zcoor,Showxz);
        colormap(handles.cMap); caxis([-handles.cRange handles.cRange]); colorbar
        axis equal tight
        title([handles.Title, ', XZ, iY=',num2str(handles.slt(2))])
        xlabel('X [um]'); ylabel('Z [um]')
        
        subplot(2,2,3) % XY enface view MIP
        Vsign=sign(squeeze(mean(handles.Vz(handles.startZ:handles.endZ,:,:),1)));
        Vsign(Vsign~=-1)=1;
        R=handles.R;
        R(R>0.1)=1;
        R(R<1)=0;
        showMIP=(squeeze(max(abs(handles.DataShow(handles.startZ:handles.endZ,:,:).*R(handles.startZ:handles.endZ,:,:)),[],1)).*Vsign);
        showMIP(handles.slt(2)+[0 1],:)=max(showMIP(:));
        showMIP(:,handles.slt(3)+[0 1])=max(showMIP(:));
        imagesc(handles.Xcoor,handles.Ycoor,showMIP);
        colormap(handles.cMap); caxis([-handles.cRange handles.cRange]); colorbar
        axis equal tight
        title([handles.Title, ', V MIP [',num2str([handles.startZ handles.endZ]), '] [mm/s]'])
        xlabel('X [um]'); ylabel('Y [um]')
        
        subplot(2,2,4) %yz side view
        Showyz=squeeze(handles.DataShow(:,:,handles.slt(3))).*sign(squeeze(handles.Vz(:,:,handles.slt(3))));
        Showyz(handles.slt(1),:)=max(Showyz(:));
        Showyz(:,handles.slt(2))=max(Showyz(:));
        imagesc(handles.Xcoor,handles.Zcoor,Showyz);
        colormap(handles.cMap); caxis([-handles.cRange handles.cRange]); colorbar
        axis equal tight
        title([handles.Title, ', YZ, iX=',num2str(handles.slt(3))])
        xlabel('Y [um]'); ylabel('Z [um]')
    end
else %% plot enface MIP only
    handles.Show=['MIP-[',num2str([handles.startZ handles.endZ]), ']']; % en face MIP show
    handles.fig=figure;
    if (strcmp(handles.DataSlt,'Ms'))
        imagesc(handles.Xcoor,handles.Ycoor,squeeze(min(abs(handles.DataShow(handles.startZ:handles.endZ,:,:)),[],1)));
        colormap(handles.cMap); caxis(handles.Caxis); colorbar
        axis equal tight
        xlabel('Y [um]'); ylabel('Z [um]')
        title([handles.Title, ' [',num2str([handles.startZ handles.endZ]), ']'])
    elseif (strcmp(handles.DataSlt,'R')) || (strcmp(handles.DataSlt,'Mf')) || (strcmp(handles.DataSlt,'D'))
        imagesc(handles.Xcoor,handles.Ycoor,squeeze(max(abs(handles.DataShow(handles.startZ:handles.endZ,:,:)),[],1)));
        colormap(handles.cMap); caxis(handles.Caxis); colorbar
        axis equal tight
        xlabel('Y [um]'); ylabel('Z [um]')
        title([handles.Title, ' [',num2str([handles.startZ handles.endZ]), ']'])
    elseif (strcmp(handles.DataSlt,'V')) || (strcmp(handles.DataSlt,'Vt')) || (strcmp(handles.DataSlt,'Vz')) || (strcmp(handles.DataSlt,'D'))
        vSign=sign(squeeze(mean(handles.Vz(handles.startZ:handles.endZ,:,:),1)));
        vSign(vSign~=-1)=1;
        R=handles.R;
        R(R>0)=1;
        R(R<1)=0;
        imagesc(handles.Xcoor,handles.Ycoor,squeeze(max(abs(handles.DataShow(handles.startZ:handles.endZ,:,:).*R(handles.startZ:handles.endZ,:,:)),[],1)).*vSign);
        colormap(handles.cMap); caxis([-handles.cRange handles.cRange]); colorbar
        axis equal tight
        xlabel('Y [um]'); ylabel('Z [um]')
        title([handles.Title, ' [',num2str([handles.startZ handles.endZ]), ']'])
    end
end
guidata(hObject, handles);

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
clc
[Vcmap, Vzcmap, Dcmap, Mfcmap, Rcmap, g1OCTAcmap]=Colormaps_DLSOCT;
zStart=str2num(get(handles.zStart,'string'));
zStack=str2num(get(handles.zStack,'string'));

[nz,nx,ny] = size(handles.V);
set(hObject,'SliderStep',[1/(nz-1), 3/(nz-1)])
set(hObject,'Max',nz)
zStart=nz-min(round(get(hObject,'Value')),nz-1);
set(handles.zStart,'string',zStart);

if (strcmp(handles.DataSlt,'Ms'))
    axes(handles.axes1)
    imagesc(squeeze(min(abs(handles.DataShow(zStart:min(zStart+zStack-1,nz),:,:)),[],1)));
    colormap(handles.cMap); caxis(handles.Caxis); colorbar
    axis equal tight
    xlabel('X [pix]'); ylabel('Y [pix]')
    title([handles.Title, ', XY, z=[',num2str(zStart),'-',num2str(min(zStart+zStack-1,nz)),']'])
elseif (strcmp(handles.DataSlt,'R')) || (strcmp(handles.DataSlt,'Mf')) || (strcmp(handles.DataSlt,'D'))
    axes(handles.axes1)
    imagesc(squeeze(max(abs(handles.DataShow(zStart:min(zStart+zStack-1,nz),:,:)),[],1)));
    colormap(handles.cMap); caxis(handles.Caxis); colorbar
    axis equal tight
    xlabel('X [pix]'); ylabel('Y [pix]')
    title([handles.Title, ', XY, z=[',num2str(zStart),'-',num2str(min(zStart+zStack-1,nz)),']'])
elseif (strcmp(handles.DataSlt,'V')) || (strcmp(handles.DataSlt,'Vt')) || (strcmp(handles.DataSlt,'Vz')) || (strcmp(handles.DataSlt,'D'))
    vSign=sign(squeeze(mean(handles.Vz(zStart:min(zStart+zStack-1,nz),:,:),1)));
    vSign(vSign~=-1)=1;
    axes(handles.axes1)
    imagesc(squeeze(max(abs(handles.DataShow(zStart:min(zStart+zStack-1,nz),:,:)),[],1)).*vSign);
    colormap(handles.cMap); caxis(handles.Caxis); colorbar
    axis equal tight
    xlabel('X [pix]'); ylabel('Y [pix]')
    title([handles.Title, ', XY, z=[',num2str(zStart),'-',num2str(min(zStart+zStack-1,nz)),']'])
end
guidata(hObject, handles);

function zStart_Callback(hObject, eventdata, handles)
% hObject    handle to zStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zStart as text
%        str2double(get(hObject,'String')) returns contents of zStart as a double
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
[Vcmap, Vzcmap, Dcmap, Mfcmap, Rcmap, g1OCTAcmap]=Colormaps_DLSOCT;
zStart=str2num(get(handles.zStart,'string'));
zStack=str2num(get(handles.zStack,'string'));
[nz,nx,ny]=size(handles.V);

if (strcmp(handles.DataSlt,'Ms'))
    axes(handles.axes1)
    imagesc(squeeze(min(abs(handles.DataShow(zStart:min(zStart+zStack-1,nz),:,:)),[],1)));
    colormap(handles.cMap); caxis(handles.Caxis); colorbar
    axis equal tight
    xlabel('X [pix]'); ylabel('Y [pix]')
    title([handles.Title, ', XY, z=[',num2str(zStart),'-',num2str(min(zStart+zStack-1,nz)),']'])
elseif (strcmp(handles.DataSlt,'R')) || (strcmp(handles.DataSlt,'Mf')) || (strcmp(handles.DataSlt,'D'))
    axes(handles.axes1)
    imagesc(squeeze(max(abs(handles.DataShow(zStart:min(zStart+zStack-1,nz),:,:)),[],1)));
    colormap(handles.cMap); caxis(handles.Caxis); colorbar
    axis equal tight
    xlabel('X [pix]'); ylabel('Y [pix]')
    title([handles.Title, ', XY, z=[',num2str(zStart),'-',num2str(min(zStart+zStack-1,nz)),']'])
elseif (strcmp(handles.DataSlt,'V')) || (strcmp(handles.DataSlt,'Vt')) || (strcmp(handles.DataSlt,'Vz')) || (strcmp(handles.DataSlt,'D'))
    vSign=sign(squeeze(mean(handles.Vz(zStart:min(zStart+zStack-1,nz),:,:),1)));
    vSign(vSign~=-1)=1;
    axes(handles.axes1)
    imagesc(squeeze(max(abs(handles.DataShow(zStart:min(zStart+zStack-1,nz),:,:)),[],1)).*vSign);
    colormap(handles.cMap); caxis(handles.Caxis); colorbar
    axis equal tight
    xlabel('X [pix]'); ylabel('Y [pix]')
    title([handles.Title, ', XY, z=[',num2str(zStart),'-',num2str(min(zStart+zStack-1,nz)),']'])
end

function zStack_Callback(hObject, eventdata, handles)
% hObject    handle to zStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zStack as text
%        str2double(get(hObject,'String')) returns contents of zStack as a double
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
[Vcmap, Vzcmap, Dcmap, Mfcmap, Rcmap, g1OCTAcmap]=Colormaps_DLSOCT;
zStart=str2num(get(handles.zStart,'string'));
zStack=str2num(get(handles.zStack,'string'));
[nz,nx,ny]=size(handles.V);

if (strcmp(handles.DataSlt,'Ms'))
    axes(handles.axes1)
    imagesc(squeeze(min(abs(handles.DataShow(zStart:min(zStart+zStack-1,nz),:,:)),[],1)));
    colormap(handles.cMap); caxis(handles.Caxis); colorbar
    axis equal tight
    xlabel('X [pix]'); ylabel('Y [pix]')
    title([handles.Title, ', XY, z=[',num2str(zStart),'-',num2str(min(zStart+zStack-1,nz)),']'])
elseif (strcmp(handles.DataSlt,'R')) || (strcmp(handles.DataSlt,'Mf')) || (strcmp(handles.DataSlt,'D'))
    axes(handles.axes1)
    imagesc(squeeze(max(abs(handles.DataShow(zStart:min(zStart+zStack-1,nz),:,:)),[],1)));
    colormap(handles.cMap); caxis(handles.Caxis); colorbar
    axis equal tight
    xlabel('X [pix]'); ylabel('Y [pix]')
    title([handles.Title, ', XY, z=[',num2str(zStart),'-',num2str(min(zStart+zStack-1,nz)),']'])
elseif (strcmp(handles.DataSlt,'V')) || (strcmp(handles.DataSlt,'Vt')) || (strcmp(handles.DataSlt,'Vz')) || (strcmp(handles.DataSlt,'D'))
    vSign=sign(squeeze(mean(handles.Vz(zStart:min(zStart+zStack-1,nz),:,:),1)));
    vSign(vSign~=-1)=1;
    axes(handles.axes1)
    imagesc(squeeze(max(abs(handles.DataShow(zStart:min(zStart+zStack-1,nz),:,:)),[],1)).*vSign);
    colormap(handles.cMap); caxis(handles.Caxis); colorbar
    axis equal tight
    xlabel('X [pix]'); ylabel('Y [pix]')
    title([handles.Title, ', XY, z=[',num2str(zStart),'-',num2str(min(zStart+zStack-1,nz)),']'])
end


function btn_Save_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
%% select file path %%%%%%%%%
defaultpath=handles.defpath;
disp('Saving data......')
saveas(handles.fig,[defaultpath,handles.DataSlt,'-', handles.Show, '.fig']);
saveas(handles.fig,[defaultpath,handles.DataSlt,'-', handles.Show, '.jpg']);
disp('Data saved!')



% --- Executes on button press in btn_reset.
function btn_reset_Callback(hObject, eventdata, handles)
% hObject    handle to btn_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

addpath(handles.CodePath);
addpath([handles.CodePath, '\SubFunctions'])
handles.defpath='H:\BU';

handles.startZ=1;
handles.stackZ=100;
% Choose default command line output for GUI_PLOT_DLSOCT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function zStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function zStack_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function POP_DataSlect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to POP_DataSlect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
