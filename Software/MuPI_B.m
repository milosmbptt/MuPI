function varargout = MuPI_B(varargin)
% MUPI_B M-file for MuPI_B.fig
%      MUPI_B, by itself, creates a new MUPI_B or raises the existing
%      singleton*.
%
%      H = MUPI_B returns the handle to a new MUPI_B or the handle to
%      the existing singleton*.
%
%      MUPI_B('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MUPI_B.M with the given input arguments.
%
%      MUPI_B('Property','Value',...) creates a new MUPI_B or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MuPI_B_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MuPI_B_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MuPI_B

% Last Modified by GUIDE v2.5 17-Jun-2025 07:55:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MuPI_B_OpeningFcn, ...
                   'gui_OutputFcn',  @MuPI_B_OutputFcn, ...
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


% --- Executes just before MuPI_B is made visible.
function MuPI_B_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MuPI_B (see VARARGIN)

% Choose default command line output for MuPI_B
set(hObject,'visible','off');
ekran=get(0,'screensize');
if ekran(3)>=800 
    prozor_0=get(hObject,'position');
    prozor_0(1)=(ekran(3)-prozor_0(3))/2;
    prozor_0(2)=(ekran(4)-prozor_0(4))/2;
else
    prozor_0=get(hObject,'position');
    prozor_0=0.8*ekran;
    prozor_0(1)=(ekran(3)-prozor_0(3))/2;
    prozor_0(2)=(ekran(4)-prozor_0(4))/2;    
end

handles.output = hObject;

load asa_banka;
ul_sig=struct('ime','ime','fs',44100,'x',1:100);
sig_banka=[];
set(handles.Projektuj_banku,'userdata',sig_banka);
fname=signal_fajl.fname;
pname=signal_fajl.pname;
ime_filea=strcat(pname,fname);
ul_sig.ime=ime_filea;
set(handles.signal_fajl,'string',ime_filea);
[x,fs] = audioread(ime_filea);
x=x-mean(x);
if rem(length(x),2)~=0
 x(end+1,:)=0;
end
ul_sig.fs=fs;
ul_sig.x=x;
ul_sig.ime=ime_filea;




f0=100;
set(handles.podesi_f0,'value',f0);
set(handles.podesi_f0,'min',f0/sqrt(2));
set(handles.podesi_f0,'max',f0*sqrt(2));
%set(handles.podesi_f0,'sliderstep',[f0/10000 f0/1000]);
set(handles.podesi_f0_str,'value',f0);
set(handles.podesi_f0_str,'string',num2str(f0));
set(handles.oznaci_osnovnu_fr,'userdata',f0);
set(hObject,'position',prozor_0,'visible','on');
set(handles.ucitaj_banku,'enable','on');
set(handles.Nova_banka,'enable','on');
set(handles.Projektuj_banku,'enable','off');
set(handles.As,'enable','on');
set(handles.wp,'enable','on');
set(handles.Sacuvaj_banku,'enable','off');
set(handles.Procesiraj_signal,'enable','off');
set(handles.CBank,'value',0);
if length(p_p)==1
    set(handles.s1,'enable','off');
else
    set(handles.s1,'enable','on');
end
set(handles.oAxD,'enable','off');
set(handles.next,'enable','off');
set(handles.signal_fajl,'userdata',ul_sig);
axes(handles.crtez);
crtaj_signal(handles);
xlim([0 fs/2]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MuPI_B wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MuPI_B_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ucitaj_banku.
function ucitaj_banku_Callback(hObject, eventdata, handles)
% hObject    handle to ucitaj_banku (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,pname]=uigetfile('*.mat', 'Banka');
if ~(fname==0)
	ime_filea=strcat(pname,fname);
    load(ime_filea);
    ul_sig=get(handles.signal_fajl,'userdata');
    fs=ul_sig.fs;
    if banka.fs~=fs
        iz_poruka=Poruka_1('title','Error','string','The sampling frequency of the bank is different from the signal sampling frequency!!!');
    else
        fc=banka.fc;
        set(handles.As,'value',banka.As);
        set(handles.As,'string',num2str(banka.As));
        set(handles.wp,'value',banka.wp);
        set(handles.wp,'string',num2str(banka.wp));
        set(handles.tip_banke,'value',banka.tip_banke);
        sig_banka=[];
        set(handles.Projektuj_banku,'userdata',sig_banka);
        ul_sig=get(handles.signal_fajl,'userdata');
        fc=sort(fc);
        wc_uk=fc/fs*2*pi;
        impuls=zeros(10000,1);
        impuls(1)=1;
        tip_banke=get(handles.tip_banke,'value');
        sig_banka=banka_all_pass(wc_uk,get(handles.wp,'value'),get(handles.As,'value'),impuls,0,fs,tip_banke-1,1);
        SIG_banka=abs(fft(sig_banka));
        set(handles.Projektuj_banku,'userdata',SIG_banka);
        axes(handles.crtez);
        x_za_lim=get(handles.crtez,'xlim');
        hold off
        crtaj_signal(handles);
        hold on
%        crtaj_gr_tac(handles);
        set(handles.CBank,'value',1);        
        crtaj_banka(handles);
        xlim(x_za_lim);
        set(handles.Sacuvaj_banku,'enable','on');
        set(handles.Procesiraj_signal,'enable','on');
        br_harmonika=length(fc)-1;
        set(handles.broj_harmonika,'value',br_harmonika);
        set(handles.broj_harmonika,'string',num2str(br_harmonika));
        f0=fc(1)*sqrt(2);
        B=(-5+sqrt(25-16*(1-fc(2)^4/(4*f0^4))))/8;
        set(handles.podesi_f0,'min',f0/sqrt(2));
        set(handles.podesi_f0,'max',f0*sqrt(2));
        set(handles.podesi_f0,'value',f0);
        set(handles.podesi_f0_str,'string',num2str(f0));
        set(handles.podesi_f0_str,'value',B);
        set(handles.podesi_B,'value',B);
        set(handles.podesi_B_str,'string',num2str(B));
        set(handles.podesi_B_str,'value',B);
        fc_osnovna=get(handles.oznaci_osnovnu_fr,'userdata');
        if ~isempty(fc_osnovna)
            hold on;
            crtaj_osn_fr(handles);
        end        
        if strcmp(get(handles.Projektuj_banku,'enable'),'off')
            set(handles.Projektuj_banku,'enable','on');
            set(handles.As,'enable','on');
            set(handles.wp,'enable','on');
        end
    end
end


% --- Executes on button press in Nova_banka.
function Nova_banka_Callback(hObject, eventdata, handles)
% hObject    handle to Nova_banka (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.tip_banke,'value',2);
sig_banka=[];
set(handles.Projektuj_banku,'userdata',sig_banka);
axes(handles.crtez);
x_za_lim=get(handles.crtez,'xlim');
hold off
crtaj_signal(handles);
fc_osnovna=get(handles.oznaci_osnovnu_fr,'userdata');
if ~isempty(fc_osnovna)
    hold on;
    crtaj_osn_fr(handles);
end
xlim(x_za_lim);
set(handles.Projektuj_banku,'enable','on');
set(handles.As,'enable','on');
set(handles.wp,'enable','on');
set(handles.Sacuvaj_banku,'enable','off');
set(handles.Procesiraj_signal,'enable','off');
set(handles.oAxD,'enable','off');
set(handles.Procesiraj_signal,'userdata',[]);

function As_Callback(hObject, eventdata, handles)
% hObject    handle to As (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of As as text
%        str2double(get(hObject,'String')) returns contents of As as a double
As_str=get(hObject,'string');
As=str2double(As_str);
if (isnan(As) || As<40 || As>90)
    As_str=num2str(get(hObject,'value'));
    set(hObject,'string',As_str);
else
    set(hObject,'value',As);
end


% --- Executes during object creation, after setting all properties.
function As_CreateFcn(hObject, eventdata, handles)
% hObject    handle to As (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function As_labela_Callback(hObject, eventdata, handles)
% hObject    handle to As_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of As_labela as text
%        str2double(get(hObject,'String')) returns contents of As_labela as a double


% --- Executes during object creation, after setting all properties.
function As_labela_CreateFcn(hObject, eventdata, handles)
% hObject    handle to As_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wp_Callback(hObject, eventdata, handles)
% hObject    handle to wp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wp as text
%        str2double(get(hObject,'String')) returns contents of wp as a double
wp_str=get(hObject,'string');
wp=str2double(wp_str);
if (isnan(wp) || wp<0.4 || wp>0.495)
    wp_str=num2str(get(hObject,'value'));
    set(hObject,'string',wp_str);
else
    set(hObject,'value',wp);
end

% --- Executes during object creation, after setting all properties.
function wp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wp_labela_Callback(hObject, eventdata, handles)
% hObject    handle to wp_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wp_labela as text
%        str2double(get(hObject,'String')) returns contents of wp_labela as a double


% --- Executes during object creation, after setting all properties.
function wp_labela_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wp_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Projektuj_banku.
function Projektuj_banku_Callback(hObject, eventdata, handles)
% hObject    handle to Projektuj_banku (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tip_banke=get(handles.tip_banke,'value');
%fc=get(handles.Dodaj_wc,'userdata');
f0=get(handles.podesi_f0,'value');
B=get(handles.podesi_B,'value');
N=get(handles.broj_harmonika,'value');
fcpom=(1:N+1).*f0.*sqrt(1+(1:N+1).^2*B);
fc=sqrt(fcpom(1:N).*fcpom(2:N+1));
fc=[f0/sqrt(2) fc];
ul_sig=get(handles.signal_fajl,'userdata');
fs=ul_sig.fs;
fc=sort(fc);
wc_uk=fc/fs*2*pi;
impuls=zeros(100000,1);
impuls(1)=1;
sig_banka=banka_all_pass(wc_uk,get(handles.wp,'value'),get(handles.As,'value'),impuls,1,fs,tip_banke-1,0);
SIG_banka=abs(fft(sig_banka));
set(hObject,'userdata',SIG_banka);
axes(handles.crtez);
x_za_lim=get(handles.crtez,'xlim');
hold off
crtaj_signal(handles);
hold on
set(handles.CBank,'value',1);
crtaj_banka(handles);
fc_osnovna=get(handles.oznaci_osnovnu_fr,'userdata');
if ~isempty(fc_osnovna)
    hold on;
    crtaj_osn_fr(handles);
end
xlim(x_za_lim);
set(handles.Sacuvaj_banku,'enable','on');
set(handles.Procesiraj_signal,'enable','on');

% --- Executes on button press in Sacuvaj_banku.
function Sacuvaj_banku_Callback(hObject, eventdata, handles)
% hObject    handle to Sacuvaj_banku (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ul_sig=get(handles.signal_fajl,'userdata');
za_pocetno_ime=ul_sig.ime(1:end-4);
fs=ul_sig.fs;
f0=get(handles.podesi_f0,'value');
B=get(handles.podesi_B,'value');
N=get(handles.broj_harmonika,'value');
fcpom=(1:N+1).*f0.*sqrt(1+(1:N+1).^2*B);
fc=sqrt(fcpom(1:N).*fcpom(2:N+1));
fc=[f0/sqrt(2) fc];
fc=sort(fc);
As=get(handles.As,'value');
wp=get(handles.wp,'value');
tip_banke=get(handles.tip_banke,'value');
banka=struct('fs',fs,'fc',fc,'As',As,'wp',wp,'tip_banke',tip_banke);
pocetno_ime=strcat(za_pocetno_ime,'_',datestr(now,'dd_mm_yy'),'.mat');
[fname,pname]=uiputfile('*.mat', 'Banka',pocetno_ime);
if ~(fname==0)
    ime_filea=strcat(pname,fname);
    save(ime_filea, 'banka');
end


% --- Executes on button press in Procesiraj_signal.
function Procesiraj_signal_Callback(hObject, eventdata, handles)
% hObject    handle to Procesiraj_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load asa_banka
ul_sig=get(handles.signal_fajl,'userdata');
ulaz=ul_sig.x;
fs=ul_sig.fs;
t=(0:length(ulaz)-1)'/fs;
L=round(Tprozor*fs);
L=L-rem(L,8);
FP=0.125;
R=0.125*L;
f0=get(handles.podesi_f0,'value');
B=get(handles.podesi_B,'value');
N=get(handles.broj_harmonika,'value');
fcpom=(1:N+1).*f0.*sqrt(1+(1:N+1).^2*B);
fc=sqrt(fcpom(1:N).*fcpom(2:N+1));
fc=[f0/sqrt(2) fc];
fc=sort(fc);
As=get(handles.As,'value');
wp=get(handles.wp,'value');
tip_banke=get(handles.tip_banke,'value');
banka=struct('fs',fs,'fc',fc,'As',As,'wp',wp,'tip_banke',tip_banke);
wc_uk=fc/fs*2*pi;
Np(1)=sum(t<=p_p(1));
Nk(1)=sum(t<=k_p(1));
if length(p_p)>1
    Np(2)=sum(t<=p_p(2));
    Nk(2)=sum(t<=k_p(2));
end
ul_sig.x=ulaz(Np(1):Nk(1));
izlaz_tmp=banka_all_pass(wc_uk,banka.wp,banka.As,ulaz,0,fs,banka.tip_banke-1,0);
[duz_sig,br_sig]=size(izlaz_tmp);
figd = uifigure;
figd.WindowStyle = 'modal';
dijalog = uiprogressdlg(figd,'Title','Please Wait',...
    'Message','Opening the next application');
pause(.5);
for br_comp=2:br_sig-1
    dijalog.Value =floor(100*br_comp/br_sig)/100; 
    dijalog.Message =['Processing of channel ' num2str(br_comp)];
    pause(0.1);
    izlaz(:,br_comp-1)=fazna_korekcija(wc_uk,banka.wp,banka.As,izlaz_tmp(:,br_comp),banka.tip_banke-1);
    Xo=fft(izlaz(t>=p_p(1) & t<=k_p(1),br_comp-1));
    Xo=Xo/length(Xo);
    f=(0:length(Xo)-1)/length(Xo)*fs;
    Xo_deo=Xo(f>fc(br_comp-1) & f<fc(br_comp)); 
    f_deo=f(f>fc(br_comp-1) & f<fc(br_comp));
    za_ukljuci_full(br_comp-1,1)=max(abs(Xo_deo))/mean(abs(Xo_deo));
    ukljuci_full(br_comp-1,1)=za_ukljuci_full(br_comp-1,1)>5;    
    [max_Xo_sort,ind_Xo_sort]=sort(abs(Xo_deo),'descend');
    Xo_max_full(br_comp-1,1)=max_Xo_sort(1);
    fr_harm_full(br_comp-1,1)=f_deo(ind_Xo_sort(1));
    if length(p_p)>1
        izlaz_win=izlaz(t>=p_p(2) & t<=k_p(2),br_comp-1);
        if prozor==2
            PF=window(@hann,length(izlaz_win),'periodic');
        else
            PF=window(@blackmanharris,length(izlaz_win),'periodic');
        end
        Xo=fft(izlaz_win.*PF,max(length(izlaz_win),Nfft_pom));   
        Xo=Xo/sum(PF);
        f=(0:length(Xo)-1)/length(Xo)*fs;
        Xo_deo=Xo(f>fc(br_comp-1) & f<fc(br_comp)); 
        f_deo=f(f>fc(br_comp-1) & f<fc(br_comp));  
        za_ukljuci_full(br_comp-1,2)=max(abs(Xo_deo))/mean(abs(Xo_deo));
        ukljuci_full(br_comp-1,2)=za_ukljuci_full(br_comp-1,2)>5; 
        [max_Xo_sort,ind_Xo_sort]=sort(abs(Xo_deo),'descend');
        Xo_max_full(br_comp-1,2)=max_Xo_sort(1);
        fr_harm_full(br_comp-1,2)=f_deo(ind_Xo_sort(1));
    end

    for br1=1:length(p_p)
        if prozor==1
            PF=window(@hann,L,'periodic');
        else
            PF=window(@blackmanharris,L,'periodic');
        end
        brpr(br1)=min((floor((Nk(br1)-Np(br1))/L)-1)/FP+1,50);
        for br2=1:brpr(br1)
            Nfft=max(Nfft_pom,L);
            xp=izlaz(Np(br1)+(br2-1)*R:Np(br1)+(br2-1)*R+L-1,br_comp-1);
            Xp=abs(fft(xp.*PF,Nfft));
            Xp=Xp/sum(PF);
            f=(0:length(Xp)-1)/length(Xp)*fs;
            f_deo=f(f>fc(br_comp-1) & f<fc(br_comp));
            Xp_deo=Xp(f>fc(br_comp-1) & f<fc(br_comp)); 
            [max_Xp_sort,ind_Xp_sort]=sort(abs(Xp_deo),'descend');
            Xp_max(br_comp-1,br1,br2)=max_Xp_sort(1);
            fr_harm(br_comp-1,br1,br2)=f_deo(ind_Xp_sort(1));
        end
    end
end
close(dijalog);
close(figd);
B1_uk=procena_B(1e-4,fr_harm_full(1:br_sig-2,1),br_sig-2);
for br=1:brpr(1)
    niz_za_B=squeeze(fr_harm(:,1,br));
    B1(br)=procena_B(1e-4,niz_za_B,br_sig-2);
end

if length(p_p)>1
    B2_uk=procena_B(1e-4,fr_harm_full(1:br_sig-2,2),br_sig-2);
    for br=1:brpr(2)
        niz_za_B=squeeze(fr_harm(:,2,br));
        B2(br)=procena_B(1e-4,niz_za_B,br_sig-2);
    end
end


figure
za_crt=squeeze(fr_harm(:,1,1:brpr(1)));
plot((Np(1)+(0:brpr(1)-1)*R+L/2)/fs,za_crt.*ukljuci_full(1:br_sig-2,1),'Color',[0 0.45 0.74]); hold on
plot((Np(1)+(0:brpr(1)-1)*R+L/2)/fs,za_crt.*(1-ukljuci_full(1:br_sig-2,1)),'--k');
if length(p_p)>1
    za_crt=squeeze(fr_harm(:,2,1:brpr(2)));
    plot((Np(2)+(0:brpr(2)-1)*R+L/2)/fs,za_crt.*ukljuci_full(1:br_sig-2,2),'Color',[0.8500 0.3250 0.0980]);
    plot((Np(2)+(0:brpr(2)-1)*R+L/2)/fs,za_crt.*(1-ukljuci_full(1:br_sig-2,2)),'--k');
end
xlabel('t [s]');
ylabel('f [Hz]');
%legend('Seg 1','Seg2');
title('Partial frequencies');

figure
plot((Np(1)+(0:brpr(1)-1)*R+L/2)/fs,B1,'--','Color',[0 0.45 0.74]); hold on
p1=plot((Np(1)+(0:brpr(1)-1)*R+L/2)/fs,ones(size((Np(1)+(0:brpr(1)-1)*R+L/2)/fs))*B1_uk,'Color',[0 0.45 0.74]);
if length(p_p)>1
    plot((Np(2)+(0:brpr(2)-1)*R+L/2)/fs,B2,'--','Color',[0.8500 0.3250 0.0980]);
    p2=plot((Np(2)+(0:brpr(2)-1)*R+L/2)/fs,ones(size((Np(2)+(0:brpr(2)-1)*R+L/2)/fs))*B2_uk,'Color',[0.8500 0.3250 0.0980]);
    legend([p1,p2],'Seg 1','Seg2');
end
xlabel('t [s]');
ylabel('B');
title('Estimated B');

maksimumi=struct('f_full',fr_harm_full,'X_full',Xo_max_full,'f_p',fr_harm,'X_p',Xp_max,'br_pr',brpr);
set(handles.Procesiraj_signal,'userdata',maksimumi);
set(handles.oAxD,'enable','on');
set(handles.next,'enable','on');
M2M=za_ukljuci_full(:,1);
save('banka_def_tmp.mat','ul_sig','banka','B','f0','M2M');

% --- Executes on button press in lin_log.
function lin_log_Callback(hObject, eventdata, handles)
% hObject    handle to lin_log (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lin_log=get(hObject,'value');
if lin_log==0
    set(hObject,'string','Lin/Log -> Lin');
else
    set(hObject,'string','Lin/Log -> Log');
end
axes(handles.crtez);
x_za_lim=get(handles.crtez,'xlim');
hold off
crtaj_signal(handles);
fc_osnovna=get(handles.oznaci_osnovnu_fr,'userdata');
if ~isempty(fc_osnovna)
    hold on;
    crtaj_osn_fr(handles);
end
sig_banka=get(handles.Projektuj_banku,'userdata');
if ~isempty(sig_banka)
    hold on;
    crtaj_banka(handles);
end
fr_oAxD=get(handles.Procesiraj_signal,'userdata');
if ~isempty(fr_oAxD)
    hold on;
    crtaj_oAxD(handles);
end
xlim(x_za_lim);

function crtaj_signal(handles)
load asa_banka
ul_sig=get(handles.signal_fajl,'userdata');
x0=ul_sig.x;
fs=ul_sig.fs;
t=(0:length(x0)-1)/fs;
x=x0(t>=p_p(1) & t<=k_p(1));
if get(handles.s1,'value')==1
    x_1=x0(t>=p_p(2) & t<=k_p(2));
else
    x_1=0;
end
X=abs(fft(x));
X=X/length(X);
XScal=max(X);
X=X/XScal;
switch prozor
    case 1
        X_1=abs(fft(x_1));
    case 2
        win=window(@hann,length(x),'periodic');
        X_p=abs(fft(x.*win,max(length(x),Nfft_pom)));
        X_p=X_p/sum(win);
        X_p=X_p/XScal;
        win=window(@hann,length(x_1),'periodic');
        X_1=abs(fft(x_1.*win,max(length(x_1),Nfft_pom)));
        X_1=X_1/sum(win);
    case 3
        win=window(@blackmanharris,length(x),'periodic');
        X_p=abs(fft(x.*win,max(length(x),Nfft_pom))); 
        X_p=X_p/sum(win);
        X_p=X_p/XScal;
        win=window(@blackmanharris,length(x_1),'periodic');
        X_1=abs(fft(x_1.*win,max(length(x_1),Nfft_pom)));
        X_1=X_1/sum(win);
end
X_1=X_1/XScal;
f=(0:length(X)-1)/length(X)*fs;
f_p=(0:length(X_p)-1)/length(X_p)*fs;
f_1=(0:length(X_1)-1)/length(X_1)*fs;
axes(handles.crtez);
lin_log=get(handles.lin_log,'value');
if lin_log==0
    plot(f(1:round(length(f)/2)),X(1:round(length(f)/2)));
    if get(handles.win,'value')==1 && prozor>1
        hold on
        uk_w=plot(f_p(1:round(length(f_p)/2)),X_p(1:round(length(f_p)/2)),'m'); hold off
        set(handles.win,'userdata',uk_w);
    end
    if get(handles.s1,'value')==1
        hold on
        uk_1s=plot(f_1(1:round(length(f_1)/2)),X_1(1:round(length(f_1)/2)),'Color',[0.8500 0.3250 0.0980]); hold off
        set(handles.s1,'userdata',uk_1s);
    end
else
    plot(f(1:round(length(f)/2)),20*log10(X(1:round(length(f)/2))));
    if get(handles.win,'value')==1 && prozor>1
        hold on
        uk_w=plot(f_p(1:round(length(f_p)/2)),20*log10(X_p(1:round(length(f_p)/2))),'m'); hold off
        set(handles.win,'userdata',uk_w);
    end
    if get(handles.s1,'value')==1 
        hold on;
        uk_1s=plot(f_1(1:round(length(f_1)/2)),20*log10(X_1(1:round(length(f_1)/2))),'Color',[0.8500 0.3250 0.0980]); hold off
        set(handles.s1,'userdata',uk_1s);
    end
end
xlabel('fs [Hz]');

function crtaj_osn_fr(handles)
axes(handles.crtez);
x_za_lim=get(handles.crtez,'xlim');
y_za_lim=get(handles.crtez,'ylim');
crt_harm=get(handles.crtez,'userdata');
delete(crt_harm);
if exist('crt_harm','var')==1
    clear crt_harm
end
f0=get(handles.podesi_f0,'value');
B=get(handles.podesi_B,'value');
N=get(handles.broj_harmonika,'value');
fc=(1:N).*f0.*sqrt(1+(1:N).^2*B);
if N>2   
    br=1;
    while fc(br)+fc(br+1)<=fc(N)
        fc_d1(br)=fc(br)+fc(br+1);
        br=br+1;
    end
else
    fc_d1=fc(1);
end
if N>2   
    br=1;
    while 2*fc(br)<=fc(N)
        fc_x2(br)=2*fc(br);
        br=br+1;
    end
else
    fc_x2=2*fc(1);
end
fcc(1,:)=fc;
fcc(2,:)=fc;
fcc_d1(1,:)=fc_d1;
fcc_d1(2,:)=fc_d1;
fcc_x2(1,:)=fc_x2;
fcc_x2(2,:)=fc_x2;
lin_log=get(handles.lin_log,'value');
% if lin_log==0
	crt=y_za_lim(1)*ones(size(fcc));
	crt(2,:)=y_za_lim(2);
    crt_d1=y_za_lim(1)*ones(size(fcc_d1));
	crt_d1(2,:)=y_za_lim(2);
    crt_x2=y_za_lim(1)*ones(size(fcc_x2));
	crt_x2(2,:)=y_za_lim(2);
	crt_harm(1:length(fcc),1)=plot(fcc,crt,'r--');
    if get(handles.d1,'value')==1
        crt_harm(1:length(fcc_d1),2)=plot(fcc_d1,crt_d1,'g--');
    end
    if get(handles.x2,'value')==1
        crt_harm(1:length(fcc_x2),4)=plot(fcc_x2,crt_x2,'c--');
    end
set(handles.crtez,'userdata',crt_harm);
xlim(x_za_lim);


function crtaj_banka(handles)
if get(handles.CBank,'value')==1
ul_sig=get(handles.signal_fajl,'userdata');
SIG_banka=get(handles.Projektuj_banku,'userdata');
fs=ul_sig.fs;
f=(0:length(SIG_banka)-1)/length(SIG_banka)*fs;
axes(handles.crtez);
lin_log=get(handles.lin_log,'value');
if lin_log==0
    plot(f(1:round(length(f)/2)),SIG_banka(1:round(length(f)/2),:),'k--');
else
    plot(f(1:round(length(f)/2)),20*log10(SIG_banka(1:round(length(f)/2),:)),'k--');
end
end
xlabel('fs [Hz]');


% --- Executes on button press in oznaci_osnovnu_fr.
function oznaci_osnovnu_fr_Callback(hObject, eventdata, handles)
% hObject    handle to oznaci_osnovnu_fr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.crtez);
x_za_lim=get(handles.crtez,'xlim');
hold off
crtaj_signal(handles);
sig_banka=get(handles.Projektuj_banku,'userdata');
if ~isempty(sig_banka)
    hold on;
    crtaj_banka(handles);
end
xlim(x_za_lim);
[fc_osnovna,nista] = ginput(1);
set(hObject,'userdata',fc_osnovna);
rbr_h=get(handles.ozn_harmonik,'value');
B=get(handles.podesi_B,'value');
f0=fc_osnovna/rbr_h/sqrt(1+rbr_h^2*B);
set(handles.podesi_f0,'value',f0);
set(handles.podesi_f0,'min',f0/sqrt(2));
set(handles.podesi_f0,'max',f0*sqrt(2));
set(handles.podesi_f0_str,'value',f0);
set(handles.podesi_f0_str,'string',num2str(f0));
axes(handles.crtez);
hold on;
crtaj_osn_fr(handles);
f_ispis=(1:5).*f0.*sqrt(1+(1:5).^2*B);
ispis=['f1=' num2str(f_ispis(1)) ', f2=' num2str(f_ispis(2)) ', f3=' num2str(f_ispis(3)) ', f4=' num2str(f_ispis(4)) ', f5=' num2str(f_ispis(5)) ' [Hz]'];  
set(handles.fr_list,'string',ispis);
f_C3=440*2^(-9/12);
fc=fc_osnovna;
n=round(12*log(fc/f_C3)/log(2));
if n>=0
    okt=floor(n/12);
    ton=n-okt*12;
    okt=okt+4;
else
    okt=floor(abs(n)/12);
    ton=12+n+okt*12;
    okt=3-okt;
    if ton==12
        ton=0;
        okt=okt+1;
    end
end
note_sve{1}='C';
note_sve{2}='C#';
note_sve{3}='D';
note_sve{4}='D#';
note_sve{5}='E';
note_sve{6}='F';
note_sve{7}='F#';
note_sve{8}='G';
note_sve{9}='G#';
note_sve{10}='A';
note_sve{11}='A#';
note_sve{12}='H';
note_ispis=[' ' note_sve{ton+1} num2str(okt)];
ispis=['Frequency ' note_ispis];
set(hObject,'string',ispis);

% --- Executes on selection change in ozn_harmonik.
function ozn_harmonik_Callback(hObject, eventdata, handles)
% hObject    handle to ozn_harmonik (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ozn_harmonik contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ozn_harmonik
axes(handles.crtez);
x_za_lim=get(handles.crtez,'xlim');
hold off
crtaj_signal(handles);
sig_banka=get(handles.Projektuj_banku,'userdata');
if ~isempty(sig_banka)
    hold on;
    crtaj_banka(handles);
end
xlim(x_za_lim);
fc_osnovna = get(handles.oznaci_osnovnu_fr,'userdata');
rbr_h=get(hObject,'value');
B=get(handles.podesi_B,'value');
f0=fc_osnovna/rbr_h/sqrt(1+rbr_h^2*B);
set(handles.podesi_f0,'value',f0);
set(handles.podesi_f0,'min',f0/sqrt(2));
set(handles.podesi_f0,'max',f0*sqrt(2));
%set(handles.podesi_f0,'sliderstep',[f0/10000 f0/1000]);
set(handles.podesi_f0,'enable','on');
set(handles.podesi_f0,'visible','on');
set(handles.podesi_f0_str,'value',fc_osnovna);
set(handles.podesi_f0_str,'string',num2str(fc_osnovna));
set(handles.podesi_f0_str,'enable','on');
set(handles.podesi_f0_str,'visible','on');
axes(handles.crtez);
hold on;
crtaj_osn_fr(handles);
f_ispis=(1:5).*f0.*sqrt(1+(1:5).^2*B);
ispis=['f1=' num2str(f_ispis(1)) ', f2=' num2str(f_ispis(2)) ', f3=' num2str(f_ispis(3)) ', f4=' num2str(f_ispis(4)) ', f5=' num2str(f_ispis(5)) ' [Hz]'];  
set(handles.fr_list,'string',ispis);
f_C3=440*2^(-9/12);
n=round(12*log(fc_osnovna/f_C3)/log(2));
if n>=0
    okt=floor(n/12);
    ton=n-okt*12;
    okt=okt+4;
else
    okt=floor(abs(n)/12);
    ton=12+n+okt*12;
    okt=3-okt;
    if ton==12
        ton=0;
        okt=okt+1;
    end
end
note_sve{1}='C';
note_sve{2}='C#';
note_sve{3}='D';
note_sve{4}='D#';
note_sve{5}='E';
note_sve{6}='F';
note_sve{7}='F#';
note_sve{8}='G';
note_sve{9}='G#';
note_sve{10}='A';
note_sve{11}='A#';
note_sve{12}='H';
note_ispis=[' ' note_sve{ton+1} num2str(okt)];
ispis=get(hObject,'string');
ispis=['Frekvencija ' note_ispis];
set(handles.oznaci_osnovnu_fr,'string',ispis);

% --- Executes during object creation, after setting all properties.
function ozn_harmonik_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ozn_harmonik (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fr_list_Callback(hObject, eventdata, handles)
% hObject    handle to fr_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fr_list as text
%        str2double(get(hObject,'String')) returns contents of fr_list as a double


% --- Executes during object creation, after setting all properties.
function fr_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fr_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tip_banke.
function tip_banke_Callback(hObject, eventdata, handles)
% hObject    handle to tip_banke (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tip_banke contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tip_banke


% --- Executes during object creation, after setting all properties.
function tip_banke_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tip_banke (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function podesi_f0_Callback(hObject, eventdata, handles)
% hObject    handle to podesi_f0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fc_osnovna=get(handles.oznaci_osnovnu_fr,'userdata');
if ~isempty(fc_osnovna)
    hold on;
    crtaj_osn_fr(handles);
end
f0=get(hObject,'value');
set(handles.podesi_f0_str,'string',num2str(f0));
set(handles.podesi_f0_str,'value',f0);
B=get(handles.podesi_B,'value');
f_ispis=(1:5).*f0.*sqrt(1+(1:5).^2*B);
ispis=['f1=' num2str(f_ispis(1)) ', f2=' num2str(f_ispis(2)) ', f3=' num2str(f_ispis(3)) ', f4=' num2str(f_ispis(4)) ', f5=' num2str(f_ispis(5)) ' [Hz]'];  
set(handles.fr_list,'string',ispis);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function podesi_f0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to podesi_f0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function broj_harmonika_Callback(hObject, eventdata, handles)
% hObject    handle to broj_harmonika (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ul_sig=get(handles.signal_fajl,'userdata');
fs=ul_sig.fs;
f0=get(handles.podesi_f0,'value');
br_harmonika_str=get(hObject,'string');
br_harmonika=floor(str2double(br_harmonika_str));
if br_harmonika>fs/2/f0
    br_harmonika=floor(fs/2/f0);
end
set(hObject,'string',num2str(br_harmonika));
set(hObject,'value',br_harmonika);
fc_osnovna=get(handles.oznaci_osnovnu_fr,'userdata');
if ~isempty(fc_osnovna)
    hold on;
    crtaj_osn_fr(handles);
end
% Hints: get(hObject,'String') returns contents of broj_harmonika as text
%        str2double(get(hObject,'String')) returns contents of broj_harmonika as a double


% --- Executes during object creation, after setting all properties.
function broj_harmonika_CreateFcn(hObject, eventdata, handles)
% hObject    handle to broj_harmonika (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function podesi_B_Callback(hObject, eventdata, handles)
% hObject    handle to podesi_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fc_osnovna=get(handles.oznaci_osnovnu_fr,'userdata');
if ~isempty(fc_osnovna)
    hold on;
    crtaj_osn_fr(handles);
end
B=get(hObject,'value');
set(handles.podesi_B_str,'string',num2str(B));
set(handles.podesi_B_str,'value',B);
f0=get(handles.podesi_f0,'value');
f_ispis=(1:5).*f0.*sqrt(1+(1:5).^2*B);
ispis=['f1=' num2str(f_ispis(1)) ', f2=' num2str(f_ispis(2)) ', f3=' num2str(f_ispis(3)) ', f4=' num2str(f_ispis(4)) ', f5=' num2str(f_ispis(5)) ' [Hz]'];  
set(handles.fr_list,'string',ispis);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function podesi_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to podesi_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in X_rezoom.
function X_rezoom_Callback(hObject, eventdata, handles)
% hObject    handle to X_rezoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.crtez);
hold off
crtaj_signal(handles);
fc_osnovna=get(handles.oznaci_osnovnu_fr,'userdata');
if ~isempty(fc_osnovna)
    hold on;
    crtaj_osn_fr(handles);
end
sig_banka=get(handles.Projektuj_banku,'userdata');
if ~isempty(sig_banka)
    hold on;
    crtaj_banka(handles);
end
fr_oAxD=get(handles.Procesiraj_signal,'userdata');
if ~isempty(fr_oAxD)
    hold on;
    crtaj_oAxD(handles);
end
ul_sig=get(handles.signal_fajl,'userdata');
xlim([0 ul_sig.fs/2]);



function podesi_f0_str_Callback(hObject, eventdata, handles)
% hObject    handle to podesi_f0_str (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of podesi_f0_str as text
%        str2double(get(hObject,'String')) returns contents of podesi_f0_str as a double
f0=str2double(get(hObject,'string'));
ul_sig=get(handles.signal_fajl,'userdata');
if (f0<=0 || f0>=0.8*ul_sig.fs/2 || isnan(f0))
    f0=get(handles.podesi_f0,'value');
    set(hObject,'string',num2str(f0));
else
    set(hObject,'value',f0);
    set(handles.podesi_f0,'min',f0/sqrt(2));
    set(handles.podesi_f0,'max',f0*sqrt(2));
    set(handles.podesi_f0,'value',f0);
    set(handles.podesi_f0,'string',num2str(f0));
end
fc_osnovna=get(handles.oznaci_osnovnu_fr,'userdata');
if ~isempty(fc_osnovna)
    hold on;
    crtaj_osn_fr(handles);
end

B=get(handles.podesi_B,'value');
f_ispis=(1:5).*f0.*sqrt(1+(1:5).^2*B);
ispis=['f1=' num2str(f_ispis(1)) ', f2=' num2str(f_ispis(2)) ', f3=' num2str(f_ispis(3)) ', f4=' num2str(f_ispis(4)) ', f5=' num2str(f_ispis(5)) ' [Hz]'];  
set(handles.fr_list,'string',ispis);



function podesi_B_str_Callback(hObject, eventdata, handles)
% hObject    handle to podesi_B_str (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of podesi_B_str as text
%        str2double(get(hObject,'String')) returns contents of podesi_B_str as a double
B=str2double(get(hObject,'string'));
if (B<0 || B>get(handles.podesi_B,'max') || isnan(B))
    B=get(handles.podesi_B,'value');
    set(hObject,'string','num2str(B)');
else
set(hObject,'value',B);
set(handles.podesi_B,'value',B);
set(handles.podesi_B,'string',num2str(B));
end
fc_osnovna=get(handles.oznaci_osnovnu_fr,'userdata');
if ~isempty(fc_osnovna)
    hold on;
    crtaj_osn_fr(handles);
end

B=get(hObject,'value');
set(handles.podesi_B_str,'string',num2str(B));
f0=get(handles.podesi_f0,'value');
f_ispis=(1:5).*f0.*sqrt(1+(1:5).^2*B);
ispis=['f1=' num2str(f_ispis(1)) ', f2=' num2str(f_ispis(2)) ', f3=' num2str(f_ispis(3)) ', f4=' num2str(f_ispis(4)) ', f5=' num2str(f_ispis(5)) ' [Hz]'];  
set(handles.fr_list,'string',ispis);


% --- Executes on button press in cuvaj.
function cuvaj_Callback(hObject, eventdata, handles)
% hObject    handle to cuvaj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ul_sig=get(handles.signal_fajl,'userdata');
za_pocetno_ime=ul_sig.ime(1:end-4);
f0=get(handles.podesi_f0,'value');
B=get(handles.podesi_B,'value');
N=get(handles.broj_harmonika,'value');
za_cuvanje(1,1)=f0;
za_cuvanje(2,1)=B;
za_cuvanje(3,1)=N;
pocetno_ime=strcat(za_pocetno_ime,'_',datestr(now,'dd_mm_yy'),'_podaci.csv');
[fname,pname]=uiputfile('*.csv', 'Podaci',pocetno_ime);
if ~(fname==0)
    ime_filea=strcat(pname,fname);
    dlmwrite(ime_filea,za_cuvanje,'delimiter',';','precision',20);
end


% --- Executes on button press in checkbox1.
function d1_Callback(hObject, eventdata, handles)
% hObject    handle to d1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of d1
crtaj_osn_fr(handles);

% --- Executes on button press in win.
function win_Callback(hObject, eventdata, handles)
% hObject    handle to win (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of win
if get(hObject,'value')==1
    axes(handles.crtez);
    x_za_lim=get(handles.crtez,'xlim');
    hold off
    crtaj_signal(handles);
    hold on
    sig_banka=get(handles.Projektuj_banku,'userdata');
    if ~isempty(sig_banka)
        hold on;
        crtaj_banka(handles);
    end
    fr_oAxD=get(handles.Procesiraj_signal,'userdata');
    if ~isempty(fr_oAxD)
        hold on;
        crtaj_oAxD(handles);
    end
    fc_osnovna=get(handles.oznaci_osnovnu_fr,'userdata');
    if ~isempty(fc_osnovna)
        hold on;
        crtaj_osn_fr(handles);
    end
    xlim(x_za_lim);
else
    uk_w=get(hObject,'userdata');
    delete(uk_w);
end

% --- Executes on button press in x2.
function x2_Callback(hObject, eventdata, handles)
% hObject    handle to x2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of x2
crtaj_osn_fr(handles);


% --- Executes on button press in s1.
function s1_Callback(hObject, eventdata, handles)
% hObject    handle to s1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of s1
if get(hObject,'value')==1
    axes(handles.crtez);
    x_za_lim=get(handles.crtez,'xlim');
    hold off
    crtaj_signal(handles);
    hold on
    sig_banka=get(handles.Projektuj_banku,'userdata');
    if ~isempty(sig_banka)
        hold on;
        crtaj_banka(handles);
    end
    fr_oAxD=get(handles.Procesiraj_signal,'userdata');
    if ~isempty(fr_oAxD)
        hold on;
        crtaj_oAxD(handles);
    end
    fc_osnovna=get(handles.oznaci_osnovnu_fr,'userdata');
    if ~isempty(fc_osnovna)
        hold on;
        crtaj_osn_fr(handles);
    end
    xlim(x_za_lim);
else
    uk_1s=get(hObject,'userdata');
    delete(uk_1s);
end


% --- Executes on button press in CBank.
function CBank_Callback(hObject, eventdata, handles)
% hObject    handle to CBank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.crtez);
x_za_lim=get(handles.crtez,'xlim');
hold off
crtaj_signal(handles);
fc_osnovna=get(handles.oznaci_osnovnu_fr,'userdata');
if ~isempty(fc_osnovna)
    hold on;
    crtaj_osn_fr(handles);
end
sig_banka=get(handles.Projektuj_banku,'userdata');
if ~isempty(sig_banka)
    hold on;
    crtaj_banka(handles);
end
fr_oAxD=get(handles.Procesiraj_signal,'userdata');
if ~isempty(fr_oAxD)
    hold on;
    crtaj_oAxD(handles);
end
xlim(x_za_lim);
% Hint: get(hObject,'Value') returns toggle state of CBank


% --- Executes on button press in oAxD.
function oAxD_Callback(hObject, eventdata, handles)
% hObject    handle to oAxD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.crtez);
x_za_lim=get(handles.crtez,'xlim');
hold off
crtaj_signal(handles);
fc_osnovna=get(handles.oznaci_osnovnu_fr,'userdata');
if ~isempty(fc_osnovna)
    hold on;
    crtaj_osn_fr(handles);
end
sig_banka=get(handles.Projektuj_banku,'userdata');
if ~isempty(sig_banka)
    hold on;
    crtaj_banka(handles);
end
fr_oAxD=get(handles.Procesiraj_signal,'userdata');
if ~isempty(fr_oAxD)
    hold on;
    crtaj_oAxD(handles);
end
xlim(x_za_lim);
% Hint: get(hObject,'Value') returns toggle state of oAxD


function crtaj_oAxD(handles)
maksimumi=get(handles.Procesiraj_signal,'userdata');
%maksimumi=struct('f_full',fr_harm_full,'X_full',Xo_max_full,'f_p',fr_harm,'X_p',Xp_max);
fm=maksimumi.f_full;
Xm=maksimumi.X_full;
fm_p=maksimumi.f_p;
Xm_p=maksimumi.X_p;
br_pr=maksimumi.br_pr;
lin_log=get(handles.lin_log,'value');
if get(handles.oAxD,'value')==1
    if lin_log==0        
        plot(fm(:,1),Xm(:,1)/max(Xm(:,1)),'*','Color',[0 0.45 0.74]);
        fpm=squeeze(fm_p(:,1,1:br_pr(1)));
        Xpm=squeeze(Xm_p(:,1,1:br_pr(1)));
        vmax=6;
        vmin=2;
        scal_velicina=(vmin-vmax)/(br_pr(1)-1)*(1:br_pr(1))+vmax-(vmin-vmax)/(br_pr(1)-1);
        vmax=0.8;
        vmin=0.2;
        scal_boja=(vmin-vmax)/(br_pr(1)-1)*(1:br_pr(1))+vmax-(vmin-vmax)/(br_pr(1)-1);
        for br=1:br_pr(1)
            plot(fpm(:,br),Xpm(:,br)/max(Xm(:,1)),'o','Markersize',scal_velicina(br),'Color',scal_boja(br)*[0 0.45 0.74]);
        end
        if get(handles.s1,'value')==1
            plot(fm(:,2),Xm(:,2)/max(Xm(:,1)),'*','Color',[0.8500 0.3250 0.0980]);
            fpm=squeeze(fm_p(:,2,1:br_pr(2)));
            Xpm=squeeze(Xm_p(:,2,1:br_pr(2)));
            vmax=6;
            vmin=2;
            scal_velicina=(vmin-vmax)/(br_pr(1)-1)*(1:br_pr(1))+vmax-(vmin-vmax)/(br_pr(1)-1);
            vmax=0.8;
            vmin=0.2;
            scal_boja=(vmin-vmax)/(br_pr(1)-1)*(1:br_pr(1))+vmax-(vmin-vmax)/(br_pr(1)-1);
            for br=1:br_pr(2)
                plot(fpm(:,br),Xpm(:,br)/max(Xm(:,1)),'d','Markersize',scal_velicina(br),'Color',scal_boja(br)*[0.8500 0.3250 0.0980]);
            end
        end
    else
        plot(fm(:,1),20*log10(Xm(:,1)/max(Xm(:,1))),'*','Color',[0 0.4470 0.7410]);
        fpm=squeeze(fm_p(:,1,1:br_pr(1)));
        Xpm=squeeze(Xm_p(:,1,1:br_pr(1)));
        vmax=6;
        vmin=2;
        scal_velicina=(vmin-vmax)/(br_pr(1)-1)*(1:br_pr(1))+vmax-(vmin-vmax)/(br_pr(1)-1);
        vmax=0.8;
        vmin=0.2;
        scal_boja=(vmin-vmax)/(br_pr(1)-1)*(1:br_pr(1))+vmax-(vmin-vmax)/(br_pr(1)-1);
        for br=1:br_pr(1)
            plot(fpm(:,br),20*log10(Xpm(:,br)/max(Xm(:,1))),'o','Markersize',scal_velicina(br),'Color',scal_boja(br)*[0 0.45 0.74]);
        end
        [r,c]=size(fpm);
        for br=1:r
            plot(fpm(br,:),20*log10(Xpm(br,:)/max(Xm(:,1))),'--','Color',[0 0.45 0.74]);
        end
        if get(handles.s1,'value')==1
            plot(fm(:,2),20*log10(Xm(:,2)/max(Xm(:,1))),'*','Color',[0.8500 0.3250 0.0980]);
            fpm=squeeze(fm_p(:,2,1:br_pr(2)));
            Xpm=squeeze(Xm_p(:,2,1:br_pr(2)));
            vmax=6;
            vmin=2;
            scal_velicina=(vmin-vmax)/(br_pr(1)-1)*(1:br_pr(1))+vmax-(vmin-vmax)/(br_pr(1)-1);
            vmax=0.8;
            vmin=0.2;
            scal_boja=(vmin-vmax)/(br_pr(1)-1)*(1:br_pr(1))+vmax-(vmin-vmax)/(br_pr(1)-1);
            for br=1:br_pr(2)
                plot(fpm(:,br),20*log10(Xpm(:,br)/max(Xm(:,1))),'d','Markersize',scal_velicina(br),'Color',scal_boja(br)*[0.8500 0.3250 0.0980]);
            end
        end
    end
end
xlabel('fs [Hz]');
%'Color',[0.8500 0.3250 0.0980]


% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MuPI_A();

% --- Executes on button press in eksp_slika.
function eksp_slika_Callback(hObject, eventdata, handles)
% hObject    handle to eksp_slika (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fignew = figure('Visible','on'); % Invisible figure
newAxes = copyobj(handles.crtez,fignew); % Copy the appropriate axes
set(newAxes,'Position',get(groot,'DefaultAxesPosition')); % The original position is copied too, so adjust it.
