function varargout = MuPI_A(varargin)
% MUPI_A M-file for MuPI_A.fig
%      MUPI_A, by itself, creates a new MUPI_A or raises the existing
%      singleton*.
%
%      H = MUPI_A returns the handle to a new MUPI_A or the handle to
%      the existing singleton*.
%
%      MUPI_A('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MUPI_A.M with the given input arguments.
%
%      MUPI_A('Property','Value',...) creates a new MUPI_A or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MuPI_A_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MuPI_A_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MuPI_A

% Last Modified by GUIDE v2.5 20-Jun-2025 08:11:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MuPI_A_OpeningFcn, ...
                   'gui_OutputFcn',  @MuPI_A_OutputFcn, ...
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


% --- Executes just before MuPI_A is made visible.
function MuPI_A_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MuPI_A (see VARARGIN)

% Choose default command line output for MuPI_A
set(hObject,'visible','off');
ekran=get(0,'screensize');
if ekran(3)>=800 
    prozor=get(hObject,'position');
    prozor(1)=(ekran(3)-prozor(3))/2;
    prozor(2)=(ekran(4)-prozor(4))/2;
else
    prozor=get(hObject,'position');
    prozor=0.8*ekran;
    prozor(1)=(ekran(3)-prozor(3))/2;
    prozor(2)=(ekran(4)-prozor(4))/2;    
end
figd = uifigure;
figd.WindowStyle = 'modal';
dijalog = uiprogressdlg(figd,'Title','Please Wait',...
    'Message','Opening the next application');
pause(.5);
handles.output = hObject;
load('banka_def_tmp.mat');
fs=ul_sig.fs;
ulaz=ul_sig.x;
wc_uk=banka.fc/fs*2*pi;
izlaz_tmp=banka_all_pass(wc_uk,banka.wp,banka.As,ulaz,0,fs,banka.tip_banke-1,0);
snoise_tmp=banka_all_pass(wc_uk,banka.wp,banka.As,randn(size(ulaz)),0,fs,banka.tip_banke-1,0);
[duz_sig,br_sig]=size(izlaz_tmp);
t=(0:duz_sig-1)'/fs;
for brojac=1:br_sig
    dijalog.Value =floor(100*brojac/br_sig)/100; 
    dijalog.Message =['Processing of channel ' num2str(brojac)];
    pause(0.1);
    izlaz(:,brojac)=fazna_korekcija(wc_uk,banka.wp,banka.As,izlaz_tmp(:,brojac),banka.tip_banke-1);
    snoise(:,brojac)=fazna_korekcija(wc_uk,banka.wp,banka.As,snoise_tmp(:,brojac),banka.tip_banke-1);
    anv(:,brojac)=abs(hilbert(izlaz(:,brojac)));
    fitfun=' a1*((x).^(a2-1)).*(a3.^(-(x-a4)/a5));';
    [f1_fitovanje,gof,pomocni]=fit(t,anv(:,brojac),fitfun,'Startpoint',[0.5 2 exp(1) 0 2],'Lower',[1e-6 1.001 0.1 -0.5 1e-6],'Upper',[10 Inf Inf 0.5 Inf]);
    anv_fit(:,brojac)=f1_fitovanje.a1*((t).^(f1_fitovanje.a2-1)).*(f1_fitovanje.a3.^(-(t-f1_fitovanje.a4)/f1_fitovanje.a5));    
    anv_podaci(brojac,1)=f1_fitovanje.a1;
    anv_podaci(brojac,2)=f1_fitovanje.a2;
    anv_podaci(brojac,3)=f1_fitovanje.a3;
    anv_podaci(brojac,4)=f1_fitovanje.a4;
    anv_podaci(brojac,5)=f1_fitovanje.a5;
end
close(dijalog);
close(figd);
set(hObject,'position',prozor);
anv_podaci_p=anv_podaci;
dboja=1/br_sig;
boja=[(0:br_sig-1)'*dboja+dboja/2 (br_sig-1:-1:0)'*dboja+dboja/2 zeros(br_sig,1) ];
str_pom{1}=[num2str(0) ' (' num2str(0) '-' num2str(banka.fc(1)) ' [Hz])'];
for brojac=2:br_sig-1
    str_pom{brojac}=[num2str(brojac-1) ' (' num2str(banka.fc(brojac-1)) '-' num2str(banka.fc(brojac)) ' ' num2str(sqrt(banka.fc(brojac-1)*banka.fc(brojac))) ' [Hz])'];
end
if br_sig>2
    str_pom{brojac+1}=[num2str(brojac) ' (' num2str(banka.fc(brojac)) '-' num2str(fs/2) ' [Hz])'];
else
    str_pom{2}=[num2str(1) ' (' num2str(banka.fc(1)) '-' num2str(fs/2) ' [Hz])'];
end
set(handles.komponente,'string',str_pom);

axes(handles.ose_vreme);
plot(t,izlaz(:,1),'color',boja(1,:)); hold on
plot(t,anv(:,1),'k');
plot(t,-anv(:,1),'k');
plot(t,anv_fit(:,1),'k-');
plot(t,-anv_fit(:,1),'k-');
xlim([0 t(end)]),xlabel('t [s]');
hold off
IZLAZ=abs(fft(izlaz));
skal_spec=max(max(IZLAZ));
IZLAZ=IZLAZ/max(max(IZLAZ));
f=(0:length(IZLAZ)-1)/length(IZLAZ)*fs;
ukljucen_1=[ones(br_sig,1) ones(br_sig,1) zeros(br_sig,1) zeros(br_sig,1) ones(br_sig,1)]; % ukljucen, 
ukljucen_2=[zeros(br_sig,1) ones(br_sig,1) zeros(br_sig,1) zeros(br_sig,1) ones(br_sig,1)];
axes(handles.ose_sve);
hold on
for brojac=1:br_sig
    plot(f(1:round(length(f)/2))/1000,(IZLAZ(1:round(length(IZLAZ)/2),brojac)),'color',boja(brojac,:));
    Spom=fft(izlaz(:,brojac))/skal_spec;
    [~,indmax]=max(abs(Spom));
    fi(brojac)=atan2(imag(Spom(indmax)),real(Spom(indmax)));
    f0(brojac)=f(indmax);
    za_string_1{brojac}='Included into sig. 1 OS';
    za_string_2{brojac}='Include into sig. 2';
end
f0_p=f0;
xlim([0 fs/2000]),xlabel('f [kHz]');
% legend(str_pom);
axes(handles.ose_spektar);
plot(f(1:round(length(f)/2))/1000,(IZLAZ(1:round(length(IZLAZ)/2),1)),'color',boja(1,:)),xlim([0 fs/2000]),xlabel('f [kHz]');
izlaz_ceo=struct('fs',fs,'t',t,'f',f,'vr_obl',izlaz,'sp_obl',IZLAZ,'skal_spec',skal_spec,'snoise',snoise,'f0',f0,'f0_p',f0_p,'fi',fi,'anv',anv,'anv_fit',anv_fit,'anv_podaci',anv_podaci,'anv_podaci_p',anv_podaci_p,'ukljucen_1',ukljucen_1,'ukljucen_2',ukljucen_2,'boja',boja);
skal=zeros(1,br_sig);
set(handles.skal_snoise,'value',skal);
set(handles.skal_snoise,'string','0');
OScal=ones(1,br_sig);
set(handles.OScal,'value',OScal);
set(handles.OScal,'string','1');
M2M=[0;M2M;0];
set(handles.M2M,'value',M2M);
set(handles.M2M,'string',num2str(M2M(1)));
set(handles.lin_log_ose_sve,'userdata',izlaz_ceo);
set(handles.A1,'string',num2str(anv_podaci(1,1)));
set(handles.A2,'string',num2str(anv_podaci(1,2)));
set(handles.A3,'string',num2str(anv_podaci(1,3)));
set(handles.A4,'string',num2str(anv_podaci(1,4)));
set(handles.A5,'string',num2str(anv_podaci(1,5)));
set(handles.f0,'string',num2str(f0(1)));
set(handles.A1_labela,'string',['A1=' num2str(anv_podaci_p(1,1))]);
set(handles.A2_labela,'string',['A2=' num2str(anv_podaci_p(1,2))]);
set(handles.A3_labela,'string',['A3=' num2str(anv_podaci_p(1,3))]);
set(handles.A4_labela,'string',['A4=' num2str(anv_podaci_p(1,4))]);
set(handles.A5_labela,'string',['A5=' num2str(anv_podaci_p(1,5))]);
set(handles.f0_labela,'string',['f0=' num2str(f0_p(1))]);
set(handles.skal_snoise_labela,'string','Scal N');
set(handles.OScal_labela,'string','OScal');
set(handles.ukljuci_u_zbir_1,'userdata',za_string_1);
set(handles.ukljuci_u_zbir_1,'value',1);
set(handles.ukljuci_u_zbir_2,'userdata',za_string_2);
set(hObject,'visible','on');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MuPI_A wait for user response (see UIRESUME)
% uiwait(handles.MuPI_A);


% --- Outputs from this function are returned to the command line.
function varargout = MuPI_A_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in lin_log_ose_sve.
function lin_log_ose_sve_Callback(hObject, eventdata, handles)
% hObject    handle to lin_log_ose_sve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lin_log_ose_sve
lin_log=get(hObject,'value');
if lin_log==0
    set(hObject,'string','Lin/Log -> Lin');
else
    set(hObject,'string','Lin/Log -> Log');
end

izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
f=izlaz_ceo.f;
IZLAZ=izlaz_ceo.sp_obl;
[duz_sig br_sig]=size(IZLAZ);
boja=izlaz_ceo.boja;
axes(handles.ose_sve);
x_za_lim=get(handles.ose_sve,'xlim');
hold off
if lin_log==0
	plot(f(1:round(length(f)/2))/1000,(IZLAZ(1:round(length(IZLAZ)/2),1)),'color',boja(1,:));
    hold on
    for brojac=2:br_sig
        plot(f(1:round(length(f)/2))/1000,(IZLAZ(1:round(length(IZLAZ)/2),brojac)),'color',boja(brojac,:));
    end
else
	plot(f(1:round(length(f)/2))/1000,20*log10(IZLAZ(1:round(length(IZLAZ)/2),1)),'color',boja(1,:));
    hold on
    for brojac=2:br_sig
        plot(f(1:round(length(f)/2))/1000,20*log10(IZLAZ(1:round(length(IZLAZ)/2),brojac)),'color',boja(brojac,:));
    end
end
xlim(x_za_lim),xlabel('f [kHz]');


% --- Executes on selection change in komponente.
function komponente_Callback(hObject, eventdata, handles)
% hObject    handle to komponente (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns komponente contents as cell array
%        contents{get(hObject,'Value')} returns selected item from komponente
izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
ukljucen_1=izlaz_ceo.ukljucen_1;
ukljucen_2=izlaz_ceo.ukljucen_2;
br_comp=get(hObject,'value');
set(handles.ukljuci_u_zbir_1,'value',ukljucen_1(br_comp,1));
set(handles.ukljuci_u_zbir_2,'value',ukljucen_2(br_comp,1));
set(handles.crtaj_izbor,'value',1);

anv_podaci=izlaz_ceo.anv_podaci;
anv_podaci_p=izlaz_ceo.anv_podaci_p;
f0=izlaz_ceo.f0;
f0_p=izlaz_ceo.f0_p;
set(handles.A1,'string',num2str(anv_podaci(br_comp,1)));
set(handles.A2,'string',num2str(anv_podaci(br_comp,2)));
set(handles.A3,'string',num2str(anv_podaci(br_comp,3)));
set(handles.A4,'string',num2str(anv_podaci(br_comp,4)));
set(handles.A5,'string',num2str(anv_podaci(br_comp,5)));
set(handles.A5,'string',num2str(anv_podaci(br_comp,5)));
set(handles.f0,'string',num2str(f0(br_comp)));
set(handles.A1_labela,'string',['A1=' num2str(anv_podaci_p(br_comp,1))]);
set(handles.A2_labela,'string',['A2=' num2str(anv_podaci_p(br_comp,2))]);
set(handles.A3_labela,'string',['A3=' num2str(anv_podaci_p(br_comp,3))]);
set(handles.A4_labela,'string',['A4=' num2str(anv_podaci_p(br_comp,4))]);
set(handles.A5_labela,'string',['A5=' num2str(anv_podaci_p(br_comp,5))]);
set(handles.f0_labela,'string',['f0=' num2str(f0_p(br_comp))]);
skal_snoise=get(handles.skal_snoise,'value');
set(handles.skal_snoise,'string',num2str(skal_snoise(br_comp)));
M2M=get(handles.M2M,'value');
set(handles.M2M,'string',num2str(M2M(br_comp)));
OScal=get(handles.OScal,'value');
set(handles.OScal,'string',num2str(OScal(br_comp)));
za_string=get(handles.ukljuci_u_zbir_1,'userdata');
set(handles.ukljuci_u_zbir_1,'string',za_string{br_comp});
za_string=get(handles.ukljuci_u_zbir_2,'userdata');
set(handles.ukljuci_u_zbir_2,'string',za_string{br_comp});
crtaj_jedan(handles);

% --- Executes during object creation, after setting all properties.
function komponente_CreateFcn(hObject, eventdata, handles)
% hObject    handle to komponente (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ukljuci_u_zbir_1.
function ukljuci_u_zbir_1_Callback(hObject, eventdata, handles)
% hObject    handle to ukljuci_u_zbir_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ukljuci_u_zbir_1
izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
ukljucen_1=izlaz_ceo.ukljucen_1;
f0=izlaz_ceo.f0;
br_comp=get(handles.komponente,'value');
skal_snoise=get(handles.skal_snoise,'value');
OScal=get(handles.OScal,'value');
ukljucen_1(br_comp,1)=get(hObject,'value');
za_string=get(handles.ukljuci_u_zbir_1,'userdata');
if ukljucen_1(br_comp,1)==1
    ukljucen_1(br_comp,2)=get(handles.izbor_za_zbir,'value');
    ukljucen_1(br_comp,3)=f0(br_comp);
    ukljucen_1(br_comp,4)=skal_snoise(br_comp);
    ukljucen_1(br_comp,5)=OScal(br_comp);
    switch ukljucen_1(br_comp,2)
        case 1
            za_string{br_comp}='Included into sig. 1 OS';
        case 2
            za_string{br_comp}=['Included into sig. 1 OScal, Scal=' num2str(OScal(br_comp))];            
        case 3
            za_string{br_comp}=['Included into sig. 1 AT, f0=' num2str(f0(br_comp))];
        case 4
            za_string{br_comp}=['Included into sig. 1 ATN, f0=' num2str(f0(br_comp)) ', SN=' num2str(skal_snoise(br_comp))];
        case 5
            za_string{br_comp}=['Included into sig. 1 N, SN=' num2str(skal_snoise(br_comp))];
    end
else
    za_string{br_comp}='Include into sig. 1';
end
izlaz_ceo.ukljucen_1=ukljucen_1;
set(hObject,'string',za_string{br_comp});
set(hObject,'userdata',za_string);
set(handles.lin_log_ose_sve,'userdata',izlaz_ceo);


% --- Executes on button press in ukljuci_u_zbir_2.
function ukljuci_u_zbir_2_Callback(hObject, eventdata, handles)
% hObject    handle to ukljuci_u_zbir_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ukljuci_u_zbir_2
izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
ukljucen_2=izlaz_ceo.ukljucen_2;
f0=izlaz_ceo.f0;
br_comp=get(handles.komponente,'value');
skal_snoise=get(handles.skal_snoise,'value');
OScal=get(handles.OScal,'value');
ukljucen_2(br_comp,1)=get(hObject,'value');
za_string=get(handles.ukljuci_u_zbir_2,'userdata');
if ukljucen_2(br_comp,1)==1
    ukljucen_2(br_comp,2)=get(handles.izbor_za_zbir,'value');
    ukljucen_2(br_comp,3)=f0(br_comp);
    ukljucen_2(br_comp,4)=skal_snoise(br_comp);
    ukljucen_2(br_comp,5)=OScal(br_comp);
    switch ukljucen_2(br_comp,2)
        case 1
            za_string{br_comp}=['Included into sig. 2 OS'];
        case 2
            za_string{br_comp}=['Included into sig. 2 OScal, Scal=' num2str(OScal(br_comp))];              
        case 3
            za_string{br_comp}=['Included into sig. 2 AT, f0=' num2str(f0(br_comp))];
        case 4
            za_string{br_comp}=['Included into sig. 2 ATN, f0=' num2str(f0(br_comp)) ', SN=' num2str(skal_snoise(br_comp))];
        case 5
            za_string{br_comp}=['Included into sig. 2 N, SN=' num2str(skal_snoise(br_comp))];
    end
else
    za_string{br_comp}=['Include into sig. 2'];
end
izlaz_ceo.ukljucen_2=ukljucen_2;
set(hObject,'string',za_string{br_comp});
set(hObject,'userdata',za_string);
set(handles.lin_log_ose_sve,'userdata',izlaz_ceo);

% --- Executes on button press in normalizovan.
function normalizovan_Callback(hObject, eventdata, handles)
% hObject    handle to normalizovan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normalizovan
ne_norm=get(hObject,'value');
if ne_norm==0
    set(hObject,'string','No-Scal');
else
    set(hObject,'string','Scal');
end

% --- Executes on button press in pusti_jedan.
function pusti_jedan_Callback(hObject, eventdata, handles)
% hObject    handle to pusti_jedan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
crtaj_jedan(handles);
br_comp=get(handles.komponente,'value');
izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
fs=izlaz_ceo.fs;
f=izlaz_ceo.f;
f0=izlaz_ceo.f0;
fi=izlaz_ceo.fi;
izlaz=izlaz_ceo.vr_obl;
anv=izlaz_ceo.anv;
anv_fit=izlaz_ceo.anv_fit;
t=(0:length(izlaz(:,br_comp))-1)'/fs;
x0=anv_fit(:,br_comp).*cos(2*pi*f0(br_comp)*t+fi(br_comp));
snoise=izlaz_ceo.snoise;
skal=get(handles.skal_snoise,'value');
OScal=get(handles.OScal,'value');
ne_norm=get(handles.normalizovan,'value');
tip_signala=get(handles.izbor_za_zbir,'value');
if ne_norm==0
    if tip_signala==1
        p1=audioplayer(izlaz(:,br_comp),fs);
        playblocking(p1);
    elseif tip_signala==2
        p1=audioplayer(OScal(br_comp)*izlaz(:,br_comp),fs);
        playblocking(p1);    
    elseif tip_signala==3
        p1=audioplayer(x0,fs);
        playblocking(p1);
    elseif tip_signala==4
        p1=audioplayer(x0+skal(br_comp)*snoise(:,br_comp),fs);
        playblocking(p1);
    elseif tip_signala==5
        p1=audioplayer(skal(br_comp)*snoise(:,br_comp),fs);
        playblocking(p1);        
    end
else
    if tip_signala==1
        p1=audioplayer(izlaz(:,br_comp)/max(abs(izlaz(:,br_comp))),fs);
        playblocking(p1);
    elseif tip_signala==2
        p1=audioplayer(OScal(br_comp)*izlaz(:,br_comp)/max(abs(OScal(br_comp)*izlaz(:,br_comp))),fs);
        playblocking(p1);        
    elseif tip_signala==3
        p1=audioplayer(x0/max(abs(x0)),fs);
        playblocking(p1);
    elseif tip_signala==4
        p1=audioplayer((x0+skal(br_comp)*snoise(:,br_comp))/max(abs(x0+skal(br_comp)*snoise(:,br_comp))),fs);
        playblocking(p1);
    elseif tip_signala==5
        p1=audioplayer(skal(br_comp)*snoise(:,br_comp)/max(abs(skal(br_comp)*snoise(:,br_comp))),fs);
        playblocking(p1);        
    end
end

% --- Executes on button press in pusti_zbir_1.
function pusti_zbir_1_Callback(hObject, eventdata, handles)
% hObject    handle to pusti_zbir_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
fs=izlaz_ceo.fs;
izlaz=izlaz_ceo.vr_obl;
f0=izlaz_ceo.f0;
fi=izlaz_ceo.fi;
anv_fit=izlaz_ceo.anv_fit;

snoise=izlaz_ceo.snoise;

ukljucen_1=izlaz_ceo.ukljucen_1;
izlaz_uk=zeros(length(izlaz),1);
for brojac=1:length(ukljucen_1)
    if ukljucen_1(brojac,1)==1
        switch ukljucen_1(brojac,2)
            case 1
                izlaz_uk=izlaz_uk+izlaz(:,brojac);
            case 2
                izlaz_uk=izlaz_uk+ukljucen_1(brojac,5)*izlaz(:,brojac);                
            case 3
                t=(0:length(izlaz(:,brojac))-1)'/fs;
                x0=anv_fit(:,brojac).*cos(2*pi*ukljucen_1(brojac,3)*t+fi(brojac));
                izlaz_uk=izlaz_uk+x0;
            case 4
                t=(0:length(izlaz(:,brojac))-1)'/fs;
                x0=anv_fit(:,brojac).*cos(2*pi*ukljucen_1(brojac,3)*t+fi(brojac));
                izlaz_uk=izlaz_uk+x0+ukljucen_1(brojac,4)*snoise(:,brojac);
            case 5
                izlaz_uk=izlaz_uk+ukljucen_1(brojac,4)*snoise(:,brojac);
        end
    end
end
ne_norm=get(handles.normalizovan,'value');
if ne_norm==0
    p1=audioplayer(izlaz_uk,fs);
    playblocking(p1);   
else
    p1=audioplayer(izlaz_uk/max(abs(izlaz_uk)),fs);
    playblocking(p1);
end

% --- Executes on button press in pusti_zbir_2.
function pusti_zbir_2_Callback(hObject, eventdata, handles)
% hObject    handle to pusti_zbir_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
fs=izlaz_ceo.fs;
izlaz=izlaz_ceo.vr_obl;
f0=izlaz_ceo.f0;
fi=izlaz_ceo.fi;
anv_fit=izlaz_ceo.anv_fit;

snoise=izlaz_ceo.snoise;

ukljucen_2=izlaz_ceo.ukljucen_2;
izlaz_uk=zeros(length(izlaz),1);
for brojac=1:length(ukljucen_2)
    if ukljucen_2(brojac,1)==1
        switch ukljucen_2(brojac,2)
            case 1
                izlaz_uk=izlaz_uk+izlaz(:,brojac);
            case 2
                izlaz_uk=izlaz_uk+ukljucen_2(brojac,5)*izlaz(:,brojac);                
            case 3
                t=(0:length(izlaz(:,brojac))-1)'/fs;
                x0=anv_fit(:,brojac).*cos(2*pi*ukljucen_2(brojac,3)*t+fi(brojac));
                izlaz_uk=izlaz_uk+x0;
            case 4
                t=(0:length(izlaz(:,brojac))-1)'/fs;
                x0=anv_fit(:,brojac).*cos(2*pi*ukljucen_2(brojac,3)*t+fi(brojac));
                izlaz_uk=izlaz_uk+x0+ukljucen_2(brojac,4)*snoise(:,brojac);
            case 5
                izlaz_uk=izlaz_uk+ukljucen_2(brojac,4)*snoise(:,brojac);
        end
    end
end
ne_norm=get(handles.normalizovan,'value');
if ne_norm==0
    p2=audioplayer(izlaz_uk,fs);
    playblocking(p2);   
else
    p2=audioplayer(izlaz_uk/max(abs(izlaz_uk)),fs);
    playblocking(p2);
end


% --- Executes on selection change in crtaj_izbor.
function crtaj_izbor_Callback(hObject, eventdata, handles)
% hObject    handle to crtaj_izbor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns crtaj_izbor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from crtaj_izbor

switch get(hObject,'value')
    case 1
        crtaj_jedan(handles);
    case 2
        crtaj_zbir_1(handles);
    case 3
        crtaj_zbir_2(handles);
    case 4
        crtaj_snimljen(handles);
end

% --- Executes during object creation, after setting all properties.
function crtaj_izbor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to crtaj_izbor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in lin_log_ose_spektar.
function lin_log_ose_spektar_Callback(hObject, eventdata, handles)
% hObject    handle to lin_log_ose_spektar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lin_log_ose_spektar
lin_log=get(hObject,'value');
if lin_log==0
    set(hObject,'string','Lin/Log -> Lin');
else
    set(hObject,'string','Lin/Log -> Log');
end
switch get(handles.crtaj_izbor,'value')
    case 1
        crtaj_jedan(handles);
    case 2
        crtaj_zbir_1(handles);
    case 3
        crtaj_zbir_2(handles);
end

% --- Executes on button press in sacuvaj_sve.
function sacuvaj_sve_Callback(hObject, eventdata, handles)
% hObject    handle to sacuvaj_sve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('banka_def_tmp.mat');
fs=ul_sig.fs;
ime_file=ul_sig.ime;
ind_f=strfind(ime_file,'\');
if ~isempty(ind_f)
    ime_path=ime_file(1:max(ind_f));
    ime_file=ime_file(max(ind_f)+1:end-4);
    ime_file_tmp=[ime_path ime_file '_harmonici/' ime_file];
    if ~exist([ime_path ime_file '_harmonici/'],'dir')
        mkdir(ime_path,[ime_file '_harmonici/']);
    end
else
    ime_file_tmp=ime_file(1:end-4);
end
izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
izlaz=izlaz_ceo.vr_obl;
for brojac=1:size(izlaz,2)
    ime_file_s=[ime_file_tmp '_' datestr(now,'dd_mm_yy_HH_MM') '_' num2str(brojac,'%02d') '.wav'];
    audiowrite(ime_file_s,izlaz(:,brojac),fs,'BitsPerSample',24);
end

% --- Executes on button press in sacuvaj_zbir_1.
function sacuvaj_zbir_1_Callback(hObject, eventdata, handles)
% hObject    handle to sacuvaj_zbir_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('banka_def_tmp.mat');
fs=ul_sig.fs;
ime_file=ul_sig.ime;
izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
fs=izlaz_ceo.fs;
izlaz=izlaz_ceo.vr_obl;
f0=izlaz_ceo.f0;
fi=izlaz_ceo.fi;
anv_fit=izlaz_ceo.anv_fit;

snoise=izlaz_ceo.snoise;

ukljucen_1=izlaz_ceo.ukljucen_1;
izlaz_uk=zeros(length(izlaz),1);
for brojac=1:length(ukljucen_1)
    if ukljucen_1(brojac,1)==1
        switch ukljucen_1(brojac,2)
            case 1
                izlaz_uk=izlaz_uk+izlaz(:,brojac);
            case 2
                izlaz_uk=izlaz_uk+ukljucen_1(brojac,5)*izlaz(:,brojac);                
            case 3
                t=(0:length(izlaz(:,brojac))-1)'/fs;
                x0=anv_fit(:,brojac).*cos(2*pi*ukljucen_1(brojac,3)*t+fi(brojac));
                izlaz_uk=izlaz_uk+x0;
            case 4
                t=(0:length(izlaz(:,brojac))-1)'/fs;
                x0=anv_fit(:,brojac).*cos(2*pi*ukljucen_1(brojac,3)*t+fi(brojac));
                izlaz_uk=izlaz_uk+x0+ukljucen_1(brojac,4)*snoise(:,brojac);
            case 5
                izlaz_uk=izlaz_uk+ukljucen_1(brojac,4)*snoise(:,brojac);
        end
    end
end
ime_file_tmp=ime_file(1:end-4);
ime_file=[ime_file_tmp '_' datestr(now,'dd_mm_yy_HH_MM') '_sum_1.wav'];
audiowrite(ime_file,izlaz_uk,fs,'BitsPerSample',24);


% --- Executes on button press in sacuvaj_zbir_2.
function sacuvaj_zbir_2_Callback(hObject, eventdata, handles)
% hObject    handle to sacuvaj_zbir_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('banka_def_tmp.mat');
fs=ul_sig.fs;
ime_file=ul_sig.ime;
izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
fs=izlaz_ceo.fs;
izlaz=izlaz_ceo.vr_obl;
fi=izlaz_ceo.fi;
anv_fit=izlaz_ceo.anv_fit;

snoise=izlaz_ceo.snoise;

ukljucen_2=izlaz_ceo.ukljucen_2;
izlaz_uk=zeros(length(izlaz),1);
for brojac=1:length(ukljucen_2)
    if ukljucen_2(brojac,1)==1
        switch ukljucen_2(brojac,2)
            case 1
                izlaz_uk=izlaz_uk+izlaz(:,brojac);
            case 2
                izlaz_uk=izlaz_uk+ukljucen_2(brojac,5)*izlaz(:,brojac);                
            case 3
                t=(0:length(izlaz(:,brojac))-1)'/fs;
                x0=anv_fit(:,brojac).*cos(2*pi*ukljucen_2(brojac,3)*t+fi(brojac));
                izlaz_uk=izlaz_uk+x0;
            case 4
                t=(0:length(izlaz(:,brojac))-1)'/fs;
                x0=anv_fit(:,brojac).*cos(2*pi*ukljucen_2(brojac,3)*t+fi(brojac));
                izlaz_uk=izlaz_uk+x0+ukljucen_2(brojac,4)*snoise(:,brojac);
            case 5
                izlaz_uk=izlaz_uk+ukljucen_2(brojac,4)*snoise(:,brojac);
        end
    end
end
ime_file_tmp=ime_file(1:end-4);
ime_file=[ime_file_tmp '_' datestr(now,'dd_mm_yy_HH_MM') '_sum_2.wav'];
audiowrite(ime_file,izlaz_uk,fs,'BitsPerSample',24);

% --- Executes on button press in sacuvaj_zbir_uk.
function sacuvaj_zbir_uk_Callback(hObject, eventdata, handles)
% hObject    handle to sacuvaj_zbir_uk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load('banka_def_tmp.mat');
fs=ul_sig.fs;
ime_file=ul_sig.ime;
br_comp=get(handles.komponente,'value');
izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
izlaz=izlaz_ceo.vr_obl;
izlaz_uk=zeros(size(izlaz,1),1);
for brojac=1:size(izlaz,2)
    izlaz_uk=izlaz_uk+izlaz(:,brojac);
end
ime_file_tmp=ime_file(1:end-4);
ime_file=[ime_file_tmp '_' datestr(now,'dd_mm_yy_HH_MM') '_zbir_uk.wav'];
audiowrite(ime_file,izlaz_uk,fs,'BitsPerSample',24);


function crtaj_jedan(handles)
br_comp=get(handles.komponente,'value');
izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
t=izlaz_ceo.t;
f=izlaz_ceo.f;
izlaz=izlaz_ceo.vr_obl;
IZLAZ=izlaz_ceo.sp_obl;
skal_spec=izlaz_ceo.skal_spec;
anv=izlaz_ceo.anv;
anv_fit=izlaz_ceo.anv_fit;
boja=izlaz_ceo.boja;
snoise=izlaz_ceo.snoise;
f0=izlaz_ceo.f0;
fi=izlaz_ceo.fi;
axes(handles.ose_vreme);
x_za_lim=get(handles.ose_vreme,'xlim');
x0=anv_fit(:,br_comp).*cos(2*pi*f0(br_comp)*t+fi(br_comp));
skal=get(handles.skal_snoise,'value');
OScal=get(handles.OScal,'value');
x0n=x0+skal(br_comp)*snoise(:,br_comp);
Spomx=abs(fft(x0))/skal_spec;
Spomxn=abs(fft(x0n))/skal_spec;
switch get(handles.izbor_za_zbir,'value')
    case 1
        plot(t,izlaz(:,br_comp),'color',boja(br_comp,:)); hold on
        plot(t,anv(:,br_comp),'k');
        plot(t,-anv(:,br_comp),'k');
        plot(t,anv_fit(:,br_comp),'color',boja(br_comp,:));
        plot(t,-anv_fit(:,br_comp),'color',boja(br_comp,:));        
        plot(t,anv_fit(:,br_comp),'k--');
        plot(t,-anv_fit(:,br_comp),'k--');
        xlim(x_za_lim),xlabel('t [s]');
        hold off
    case 2
        plot(t,OScal(br_comp)*izlaz(:,br_comp),'color',boja(br_comp,:)); hold on
        plot(t,anv(:,br_comp),'k');
        plot(t,-anv(:,br_comp),'k');
        plot(t,anv_fit(:,br_comp),'color',boja(br_comp,:));
        plot(t,-anv_fit(:,br_comp),'color',boja(br_comp,:));        
        plot(t,anv_fit(:,br_comp),'k--');
        plot(t,-anv_fit(:,br_comp),'k--');
        xlim(x_za_lim),xlabel('t [s]');
        hold off        
    case 3
        plot(t,x0,'color',boja(br_comp,:)); hold on
        plot(t,anv(:,br_comp),'k');
        plot(t,-anv(:,br_comp),'k');
        plot(t,anv_fit(:,br_comp),'color',boja(br_comp,:));
        plot(t,-anv_fit(:,br_comp),'color',boja(br_comp,:));        
        plot(t,anv_fit(:,br_comp),'k--');
        plot(t,-anv_fit(:,br_comp),'k--');
        xlim(x_za_lim),xlabel('t [s]');
        hold off
    case 4       
        plot(t,x0n,'color',boja(br_comp,:)); hold on
        plot(t,anv(:,br_comp),'k');
        plot(t,-anv(:,br_comp),'k');
        plot(t,anv_fit(:,br_comp),'color',boja(br_comp,:));
        plot(t,-anv_fit(:,br_comp),'color',boja(br_comp,:));        
        plot(t,anv_fit(:,br_comp),'k--');
        plot(t,-anv_fit(:,br_comp),'k--');
        xlim(x_za_lim),xlabel('t [s]');
        hold off
    case 5
        plot(t,skal(br_comp)*snoise(:,br_comp),'color',boja(br_comp,:)); hold on
        plot(t,anv(:,br_comp),'k');
        plot(t,-anv(:,br_comp),'k');
        plot(t,anv_fit(:,br_comp),'color',boja(br_comp,:));
        plot(t,-anv_fit(:,br_comp),'color',boja(br_comp,:));        
        plot(t,anv_fit(:,br_comp),'k--');
        plot(t,-anv_fit(:,br_comp),'k--');
        xlim(x_za_lim),xlabel('t [s]');
        hold off
end
axes(handles.ose_spektar);
lin_log=get(handles.lin_log_ose_spektar,'value');
if lin_log==0
    plot(f(1:round(length(f)/2))/1000,(IZLAZ(1:round(length(IZLAZ)/2),br_comp)),'Color',boja(br_comp,:)); hold on    
    plot(f(1:round(length(f)/2))/1000,(Spomx(1:round(length(Spomx)/2))),'Color',0.5*boja(br_comp,:));
    plot(f(1:round(length(f)/2))/1000,(Spomx(1:round(length(Spomx)/2))),'k--');
    plot(f(1:round(length(f)/2))/1000,(abs(Spomxn(1:round(length(Spomx)/2)))),'color',boja(br_comp,:));
    plot(f(1:round(length(f)/2))/1000,(abs(Spomxn(1:round(length(Spomxn)/2)))),'k--');
    hold off
else
    plot(f(1:round(length(f)/2))/1000,20*log10(IZLAZ(1:round(length(IZLAZ)/2),br_comp)),'Color',boja(br_comp,:)); hold on
    plot(f(1:round(length(f)/2))/1000,20*log10(Spomx(1:round(length(Spomx)/2))),'Color',0.5*boja(br_comp,:));
    plot(f(1:round(length(f)/2))/1000,20*log10(Spomx(1:round(length(Spomx)/2))),'k--');
    plot(f(1:round(length(f)/2))/1000,20*log10(abs(Spomxn(1:round(length(Spomxn)/2)))),'Color',boja(br_comp,:));
    plot(f(1:round(length(f)/2))/1000,20*log10(abs(Spomxn(1:round(length(Spomxn)/2)))),'k--');
    hold off
end
xlabel('f [kHz]');

function crtaj_zbir_1(handles)
load('banka_def_tmp.mat');
fs=ul_sig.fs;
ime_file=ul_sig.ime;
izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
fs=izlaz_ceo.fs;
izlaz=izlaz_ceo.vr_obl;
f0=izlaz_ceo.f0;
fi=izlaz_ceo.fi;
anv_fit=izlaz_ceo.anv_fit;

snoise=izlaz_ceo.snoise;

ukljucen_1=izlaz_ceo.ukljucen_1;
izlaz_uk=zeros(length(izlaz),1);
t=(0:length(izlaz(:,1))-1)'/fs;
for brojac=1:length(ukljucen_1)
    if ukljucen_1(brojac,1)==1
        switch ukljucen_1(brojac,2)
            case 1
                izlaz_uk=izlaz_uk+izlaz(:,brojac);
            case 2
                izlaz_uk=izlaz_uk+ukljucen_1(brojac,5)*izlaz(:,brojac);                
            case 3
                x0=anv_fit(:,brojac).*cos(2*pi*ukljucen_1(brojac,3)*t+fi(brojac));
                izlaz_uk=izlaz_uk+x0;
            case 4                
                x0=anv_fit(:,brojac).*cos(2*pi*ukljucen_1(brojac,3)*t+fi(brojac));
                izlaz_uk=izlaz_uk+x0+ukljucen_1(brojac,4)*snoise(:,brojac);
            case 5
                izlaz_uk=izlaz_uk+ukljucen_1(brojac,4)*snoise(:,brojac);
        end
    end
end
IZLAZ=abs(fft(izlaz_uk));
IZLAZ=IZLAZ/max(IZLAZ);
f=(0:length(IZLAZ)-1)/length(IZLAZ)*fs;
boja=[1 0 0];
axes(handles.ose_vreme);
x_za_lim=get(handles.ose_vreme,'xlim');
plot(t,izlaz_uk,'color',boja);
xlim(x_za_lim),xlabel('t [s]');
axes(handles.ose_spektar);
x_za_lim=get(handles.ose_spektar,'xlim');
lin_log=get(handles.lin_log_ose_spektar,'value');
if lin_log==0
    plot(f(1:round(length(f)/2))/1000,(IZLAZ(1:round(length(IZLAZ)/2))),'color',boja);
else
    plot(f(1:round(length(f)/2))/1000,20*log10(IZLAZ(1:round(length(IZLAZ)/2))),'color',boja);
end
xlim(x_za_lim),xlabel('f [kHz]');

function crtaj_zbir_2(handles)
load('banka_def_tmp.mat');
fs=ul_sig.fs;
ime_file=ul_sig.ime;
izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
fs=izlaz_ceo.fs;
izlaz=izlaz_ceo.vr_obl;
f0=izlaz_ceo.f0;
fi=izlaz_ceo.fi;
anv_fit=izlaz_ceo.anv_fit;

snoise=izlaz_ceo.snoise;

ukljucen_2=izlaz_ceo.ukljucen_2;
izlaz_uk=zeros(length(izlaz),1);
t=(0:length(izlaz(:,1))-1)'/fs;
for brojac=1:length(ukljucen_2)
    if ukljucen_2(brojac,1)==1
        switch ukljucen_2(brojac,2)
            case 1
                izlaz_uk=izlaz_uk+izlaz(:,brojac);
            case 2
                izlaz_uk=izlaz_uk+ukljucen_2(brojac,5)*izlaz(:,brojac);                
            case 3
                x0=anv_fit(:,brojac).*cos(2*pi*ukljucen_2(brojac,3)*t+fi(brojac));
                izlaz_uk=izlaz_uk+x0;
            case 4                
                x0=anv_fit(:,brojac).*cos(2*pi*ukljucen_2(brojac,3)*t+fi(brojac));
                izlaz_uk=izlaz_uk+x0+ukljucen_2(brojac,4)*snoise(:,brojac);
            case 5
                izlaz_uk=izlaz_uk+ukljucen_2(brojac,4)*snoise(:,brojac);
        end
    end
end
IZLAZ=abs(fft(izlaz_uk));
IZLAZ=IZLAZ/max(IZLAZ);
f=(0:length(IZLAZ)-1)/length(IZLAZ)*fs;
boja=[1 0 0];
axes(handles.ose_vreme);
x_za_lim=get(handles.ose_vreme,'xlim');
plot(t,izlaz_uk,'color',boja);
xlim(x_za_lim),xlabel('t [s]');
axes(handles.ose_spektar);
x_za_lim=get(handles.ose_spektar,'xlim');
lin_log=get(handles.lin_log_ose_spektar,'value');
if lin_log==0
    plot(f(1:round(length(f)/2))/1000,(IZLAZ(1:round(length(IZLAZ)/2))),'color',boja);
else
    plot(f(1:round(length(f)/2))/1000,20*log10(IZLAZ(1:round(length(IZLAZ)/2))),'color',boja);
end
xlim(x_za_lim),xlabel('f [kHz]');

function crtaj_snimljen(handles)
load('banka_def_tmp.mat');
fs=ul_sig.fs;
ulaz=ul_sig.x;
t=(0:length(ulaz)-1)'/fs;
axes(handles.ose_vreme);
x_za_lim=get(handles.ose_vreme,'xlim');
hold off        
plot(t,ulaz);
xlim(x_za_lim),xlabel('t [s]');
        

% --- Executes on button press in pusti_snimljen.
function pusti_snimljen_Callback(hObject, eventdata, handles)
% hObject    handle to pusti_snimljen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pusti_snimljen
load('banka_def_tmp.mat');
fs=ul_sig.fs;
ulaz=ul_sig.x;
ne_norm=get(handles.normalizovan,'value');
if ne_norm==0
    p1=audioplayer(ulaz,fs);
    playblocking(p1);
else
    p1=audioplayer(ulaz/max(abs(ulaz)),fs);
    playblocking(p1);
end

function A1_Callback(hObject, eventdata, handles)
% hObject    handle to A1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A1 as text
%        str2double(get(hObject,'String')) returns contents of A1 as a double
izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
anv_podaci=izlaz_ceo.anv_podaci;
anv_fit=izlaz_ceo.anv_fit;
br_comp=get(handles.komponente,'value');
str_pom=get(hObject,'string');
anv_podaci(br_comp,1)=str2double(str_pom);
izlaz_ceo.anv_podaci=anv_podaci;
t=izlaz_ceo.t;
anv_fit(:,br_comp)=anv_podaci(br_comp,1)*((t).^(anv_podaci(br_comp,2)-1)).*(anv_podaci(br_comp,3).^(-(t-anv_podaci(br_comp,4))/anv_podaci(br_comp,5))); 
izlaz_ceo.anv_fit=anv_fit;
set(handles.lin_log_ose_sve,'userdata',izlaz_ceo);
crtaj_jedan(handles);

% --- Executes during object creation, after setting all properties.
function A1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A2_Callback(hObject, eventdata, handles)
% hObject    handle to A2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A2 as text
%        str2double(get(hObject,'String')) returns contents of A2 as a double
izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
anv_podaci=izlaz_ceo.anv_podaci;
anv_fit=izlaz_ceo.anv_fit;
br_comp=get(handles.komponente,'value');
str_pom=get(hObject,'string');
anv_podaci(br_comp,2)=str2double(str_pom);
izlaz_ceo.anv_podaci=anv_podaci;
t=izlaz_ceo.t;
anv_fit(:,br_comp)=anv_podaci(br_comp,1)*((t).^(anv_podaci(br_comp,2)-1)).*(anv_podaci(br_comp,3).^(-(t-anv_podaci(br_comp,4))/anv_podaci(br_comp,5))); 
izlaz_ceo.anv_fit=anv_fit;
set(handles.lin_log_ose_sve,'userdata',izlaz_ceo);
crtaj_jedan(handles);

% --- Executes during object creation, after setting all properties.
function A2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A3_Callback(hObject, eventdata, handles)
% hObject    handle to A3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A3 as text
%        str2double(get(hObject,'String')) returns contents of A3 as a double
izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
anv_podaci=izlaz_ceo.anv_podaci;
anv_fit=izlaz_ceo.anv_fit;
br_comp=get(handles.komponente,'value');
str_pom=get(hObject,'string');
anv_podaci(br_comp,3)=str2double(str_pom);
izlaz_ceo.anv_podaci=anv_podaci;
t=izlaz_ceo.t;
anv_fit(:,br_comp)=anv_podaci(br_comp,1)*((t).^(anv_podaci(br_comp,2)-1)).*(anv_podaci(br_comp,3).^(-(t-anv_podaci(br_comp,4))/anv_podaci(br_comp,5))); 
izlaz_ceo.anv_fit=anv_fit;
set(handles.lin_log_ose_sve,'userdata',izlaz_ceo);
crtaj_jedan(handles);

% --- Executes during object creation, after setting all properties.
function A3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A4_Callback(hObject, eventdata, handles)
% hObject    handle to A4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A4 as text
%        str2double(get(hObject,'String')) returns contents of A4 as a double
izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
anv_podaci=izlaz_ceo.anv_podaci;
anv_fit=izlaz_ceo.anv_fit;
br_comp=get(handles.komponente,'value');
str_pom=get(hObject,'string');
anv_podaci(br_comp,4)=str2double(str_pom);
izlaz_ceo.anv_podaci=anv_podaci;
t=izlaz_ceo.t;
anv_fit(:,br_comp)=anv_podaci(br_comp,1)*((t).^(anv_podaci(br_comp,2)-1)).*(anv_podaci(br_comp,3).^(-(t-anv_podaci(br_comp,4))/anv_podaci(br_comp,5))); 
izlaz_ceo.anv_fit=anv_fit;
set(handles.lin_log_ose_sve,'userdata',izlaz_ceo);
crtaj_jedan(handles);

% --- Executes during object creation, after setting all properties.
function A4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A5_Callback(hObject, eventdata, handles)
% hObject    handle to A5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A5 as text
%        str2double(get(hObject,'String')) returns contents of A5 as a double
izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
anv_podaci=izlaz_ceo.anv_podaci;
anv_fit=izlaz_ceo.anv_fit;
br_comp=get(handles.komponente,'value');
str_pom=get(hObject,'string');
anv_podaci(br_comp,5)=str2double(str_pom);
izlaz_ceo.anv_podaci=anv_podaci;
t=izlaz_ceo.t;
anv_fit(:,br_comp)=anv_podaci(br_comp,1)*((t).^(anv_podaci(br_comp,2)-1)).*(anv_podaci(br_comp,3).^(-(t-anv_podaci(br_comp,4))/anv_podaci(br_comp,5))); 
izlaz_ceo.anv_fit=anv_fit;
set(handles.lin_log_ose_sve,'userdata',izlaz_ceo);
crtaj_jedan(handles);

% --- Executes during object creation, after setting all properties.
function A5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function A1_labela_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A1_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function A2_labela_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A2_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function A3_labela_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A3_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function A4_labela_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A4_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function A5_labela_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A5_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function f0_Callback(hObject, eventdata, handles)
% hObject    handle to f0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f0 as text
%        str2double(get(hObject,'String')) returns contents of f0 as a double
izlaz_ceo=get(handles.lin_log_ose_sve,'userdata');
f0=izlaz_ceo.f0;
br_comp=get(handles.komponente,'value');
str_pom=get(hObject,'string');
f0(br_comp)=str2double(str_pom);
izlaz_ceo.f0=f0;
set(handles.lin_log_ose_sve,'userdata',izlaz_ceo);
crtaj_jedan(handles);

% --- Executes during object creation, after setting all properties.
function f0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function skal_snoise_Callback(hObject, eventdata, handles)
% hObject    handle to skal_snoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of skal_snoise as text
%        str2double(get(hObject,'String')) returns contents of skal_snoise as a double
str_pom=get(hObject,'string');
br_comp=get(handles.komponente,'value');
skal=get(hObject,'value');
skal(br_comp)=str2double(str_pom);
set(hObject,'value',skal);
crtaj_jedan(handles);

% --- Executes during object creation, after setting all properties.
function skal_snoise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to skal_snoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in izbor_za_zbir.
function izbor_za_zbir_Callback(hObject, eventdata, handles)
% hObject    handle to izbor_za_zbir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns izbor_za_zbir contents as cell array
%        contents{get(hObject,'Value')} returns selected item from izbor_za_zbir
crtaj_jedan(handles);

% --- Executes during object creation, after setting all properties.
function izbor_za_zbir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to izbor_za_zbir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A1_labela_Callback(hObject, eventdata, handles)
% hObject    handle to A1_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A1_labela as text
%        str2double(get(hObject,'String')) returns contents of A1_labela as a double



function A2_labela_Callback(hObject, eventdata, handles)
% hObject    handle to A2_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A2_labela as text
%        str2double(get(hObject,'String')) returns contents of A2_labela as a double



function A3_labela_Callback(hObject, eventdata, handles)
% hObject    handle to A3_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A3_labela as text
%        str2double(get(hObject,'String')) returns contents of A3_labela as a double



function A4_labela_Callback(hObject, eventdata, handles)
% hObject    handle to A4_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A4_labela as text
%        str2double(get(hObject,'String')) returns contents of A4_labela as a double



function A5_labela_Callback(hObject, eventdata, handles)
% hObject    handle to A5_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A5_labela as text
%        str2double(get(hObject,'String')) returns contents of A5_labela as a double



function M2M_Callback(hObject, eventdata, handles)
% hObject    handle to M2M (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of M2M as text
%        str2double(get(hObject,'String')) returns contents of M2M as a double


% --- Executes during object creation, after setting all properties.
function M2M_CreateFcn(hObject, eventdata, handles)
% hObject    handle to M2M (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function OScal_Callback(hObject, eventdata, handles)
% hObject    handle to OScal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OScal as text
%        str2double(get(hObject,'String')) returns contents of OScal as a double
str_pom=get(hObject,'string');
br_comp=get(handles.komponente,'value');
OScal=get(hObject,'value');
OScal(br_comp)=str2double(str_pom);
set(hObject,'value',OScal);
crtaj_jedan(handles);

% --- Executes during object creation, after setting all properties.
function OScal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OScal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function OScal_labela_Callback(hObject, eventdata, handles)
% hObject    handle to OScal_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OScal_labela as text
%        str2double(get(hObject,'String')) returns contents of OScal_labela as a double


% --- Executes during object creation, after setting all properties.
function OScal_labela_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OScal_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in eksp_spektar.
function eksp_spektar_Callback(hObject, eventdata, handles)
% hObject    handle to eksp_slika (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fignew = figure('Visible','on'); % Invisible figure
newAxes = copyobj(handles.ose_spektar,fignew); % Copy the appropriate axes
set(newAxes,'units','normalized');
set(newAxes,'Position',get(groot,'DefaultAxesPosition')); % The original position is copied too, so adjust it.

% --- Executes on button press in eksp_vreme.
function eksp_vreme_Callback(hObject, eventdata, handles)
% hObject    handle to eksp_slika (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fignew = figure('Visible','on'); % Invisible figure
newAxes = copyobj(handles.ose_vreme,fignew); % Copy the appropriate axes
set(newAxes,'units','normalized');
set(newAxes,'Position',get(groot,'DefaultAxesPosition')); % The original position is copied too, so adjust it.
