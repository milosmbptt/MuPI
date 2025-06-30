function varargout = MuPI_S(varargin)
% MuPI_S M-file for MuPI_S.fig
%      MuPI_S, by itself, creates a new MuPI_S or raises the existing
%      singleton*.
%
%      H = MuPI_S returns the handle to a new MuPI_S or the handle to
%      the existing singleton*.
%
%      MuPI_S('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MuPI_S.M with the given input arguments.
%
%      MuPI_S('Property','Value',...) creates a new MuPI_S or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MuPI_S_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MuPI_S_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MuPI_S

% Last Modified by GUIDE v2.5 20-Jun-2025 08:02:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MuPI_S_OpeningFcn, ...
                   'gui_OutputFcn',  @MuPI_S_OutputFcn, ...
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


% --- Executes just before MuPI_S is made visible.
function MuPI_S_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MuPI_S (see VARARGIN)

% Choose default command line output for MuPI_S
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

handles.output = hObject;
set(handles.dodaj_segment,'userdata',1);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MuPI_S wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MuPI_S_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ucitaj_signal.
function ucitaj_signal_Callback(hObject, eventdata, handles)
% hObject    handle to ucitaj_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%[fname,pname]=uigetfile('*.wav', 'Ulazni signal');
pfname=get(handles.signal_fajl,'userdata');
if isempty(pfname)
    [fname,pname]=uigetfile({'*.wav', '*.mp3'}, 'Input signal');
else
    pname=pfname.pname;
    [fname,pname]=uigetfile({'*.wav', '*.mp3'}, 'Input signal', pname);
end
if ~(fname==0)
	ime_filea=strcat(pname,fname);
    ul_sig.ime=ime_filea;
    set(hObject,'userdata',ul_sig);
    [x,fs] = audioread(ime_filea);
    x=x-mean(x);
    if rem(length(x),2)~=0
     x(end+1,:)=0;
    end
    pfname.fname=fname;
    pfname.pname=pname;
    set(handles.signal_fajl,'userdata',pfname);
    set(handles.signal_fajl,'string',ime_filea);
    ul_sig=get(hObject,'userdata');
    ul_sig.fs=fs;
    ainfo=audioinfo(ime_filea);        
    ul_sig.nbits=ainfo.BitsPerSample;
    ul_sig.x=x;
    set(hObject,'userdata',ul_sig);
    axes(handles.crtez);
    hold off
    x_za_lim=[0 (length(x)-1)/fs];
    set(handles.crtez,'xlim',x_za_lim);
    crtaj_signal(handles);
    set(handles.dodaj_segment,'userdata',1);
    izbor_segmenta_str='1';
    set(handles.izbor_segmenta,'string',izbor_segmenta_str);
    set(handles.izbor_segmenta,'value',1);
    t=(0:length(x)-1)'/fs;
    p=t(1);
    k=t(end);
    set(handles.Pocetak,'userdata',p);
    set(handles.Kraj,'userdata',k);
    hold on;
    crtaj_p_k(handles); 
    x_za_lim=[0 fs/2];
    set(handles.spektar,'xlim',x_za_lim);
    crtaj_spektar(handles);
end


% --- Executes on button press in DFT_Spec.
function DFT_Spec_Callback(hObject, eventdata, handles)
% hObject    handle to DFT_Spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DFT_Spec=get(hObject,'value');
if DFT_Spec==0
    set(hObject,'string','DFT/Spec -> DFT');
    y_za_lim=[0 1.1];
    set(handles.spektar,'ylim',y_za_lim);
else
    set(hObject,'string','DFT/Spec -> Spec');
    y_za_lim=[-120 10];
    set(handles.spektar,'ylim',y_za_lim);
end
crtaj_spektar(handles);

function crtaj_spektar(handles)
ul_sig=get(handles.ucitaj_signal,'userdata');
x=ul_sig.x;
fs=ul_sig.fs;
t=(0:length(x)-1)'/fs;
p=get(handles.Pocetak,'userdata');
duzina=length(p);
k=get(handles.Kraj,'userdata');
DFT_Spec=get(handles.DFT_Spec,'value');
switch get(handles.Nfft,'value')
    case 1
        Nfft_pom=0;
    case 2
        Nfft_pom=fs;
    case 3
        Nfft_pom=10*fs;
end
boje=[0.00,0.45,0.74; 0.85,0.33,0.10; 0.93,0.69,0.13; 0.49,0.18,0.56];
axes(handles.spektar);
x_za_lim=get(handles.spektar,'xlim');
cla(handles.spektar);
if DFT_Spec==0
    br_leg=0;
    for br=1:duzina
        xo=x(t>=p(br) & t<=k(br));
        if rem(length(xo),2)~=0
            xo(end+1,:)=0;
        end
        Nfft=max(Nfft_pom,length(xo));
        Xrect=abs(fft(xo,Nfft));
        Xrect=Xrect/max(Xrect);
        prozor=window(@hann,length(xo),'periodic');
        Xhann=abs(fft(xo.*prozor,Nfft));
        Xhann=Xhann/max(Xhann);
        prozor=window(@blackmanharris,length(xo),'periodic');
        Xbh=abs(fft(xo.*prozor,Nfft));
        Xbh=Xbh/max(Xbh);
        f=(0:length(Xrect)-1)/length(Xrect)*fs;
        plot(f(1:round(length(f)/2)),20*log10(Xrect(1:round(length(f)/2))),'color',boje(br,:)); hold on
        br_leg=br_leg+1;
        za_legendu{br_leg}=['Segment ' num2str(br), ' Rect win'];
        plot(f(1:round(length(f)/2)),20*log10(Xhann(1:round(length(f)/2))),'color',boje(br,:),'LineStyle','--'); hold on
        br_leg=br_leg+1;
        za_legendu{br_leg}=['Segment ' num2str(br), ' Hann win'];
        plot(f(1:round(length(f)/2)),20*log10(Xbh(1:round(length(f)/2))),'color',boje(br,:),'LineStyle','-.'); hold on
        br_leg=br_leg+1;
        za_legendu{br_leg}=['Segment ' num2str(br), ' BlckHar win'];
    end
    if br_leg>0
        legend(za_legendu);
    end  
else
    izbor_segmenta=get(handles.izbor_segmenta,'value');
    Np=sum(t<=p(izbor_segmenta));
    Nk=sum(t<=k(izbor_segmenta));
    Tprozor=get(handles.TWIN,'value')/1000;
    L=round(Tprozor*fs);
    L=L-rem(L,8);
    FP=0.125;
    R=0.125*L;
    if get(handles.win_tip,'value')==1
        PF=window(@hann,L,'periodic');
    else
        PF=window(@blackmanharris,L,'periodic');
    end
    brpr=min((floor((Nk-Np)/L)-1)/FP+1,50);
    scal_boje=1-(1-0.2)/(brpr-1)*(0:brpr-1);
    Nfft=max(Nfft_pom,L);
    for br1=1:brpr
        xp=x(Np+(br1-1)*R:Np+(br1-1)*R+L-1);
        Xp=abs(fft(xp.*PF,Nfft));
        if br1==1
            Scal_Xp=max(Xp);
            f=(0:length(Xp)-1)/length(Xp)*fs;
        end
        Xp=Xp/Scal_Xp;
        f=(0:length(Xp)-1)/length(Xp)*fs;
        plot(f(1:round(length(f)/2)),20*log10(Xp(1:round(length(f)/2))),'color',scal_boje(br1)*boje(izbor_segmenta,:)); hold on
    end
end
hold off
xlabel('fs [Hz]');
y_za_lim=[-120 10];
set(handles.spektar,'ylim',y_za_lim);
set(handles.spektar,'xlim',x_za_lim);

function crtaj_p_k(handles)
% ul_sig=get(handles.ucitaj_signal,'userdata');
% x=ul_sig.x;
% fs=ul_sig.fs;
% t=(0:length(x)-1)'/fs;
axes(handles.crtez);
x_za_lim=get(handles.crtez,'xlim');
y_za_lim=get(handles.crtez,'ylim');
p=get(handles.Pocetak,'userdata');
k=get(handles.Kraj,'userdata');
boja=[0.00,0.45,0.74; 0.85,0.33,0.10; 0.93,0.69,0.13; 0.49,0.18,0.56];
br_segmenata=length(p);
za_legendu{1}='Signal';
for br=1:br_segmenata
    pk=[p(br) k(br)];
    ppkk(1,:)=pk;
    ppkk(2,:)=pk;
    crt=y_za_lim(1)*ones(size(ppkk));
    crt(2,:)=y_za_lim(2);
    plot(ppkk,crt,'color',boja(br,:),'LineStyle','--'); hold on
    za_legendu{2*(br-1)+2}=['Seg. ' num2str(br) ' ' num2str(p(br))];
    za_legendu{2*(br-1)+3}=['Seg. ' num2str(br) ' ' num2str(k(br))];
    % Np=sum(t<=p(br_segmenata));
    % Nk=sum(t<=k(br_segmenata));
    % t_tmp=t(Np:Nk);
    % PF=window(@hann,Nk-Np+1,'periodic');
    % plot(t_tmp,PF,'color',boja(br,:),'LineStyle',':');
end
legend(za_legendu);
hold off
xlabel('t [t]');
ylabel('x(t)');
set(handles.crtez,'xlim',x_za_lim);
set(handles.crtez,'ylim',y_za_lim);


% --- Executes on button press in Pocetak.
function Pocetak_Callback(hObject, eventdata, handles)
% hObject    handle to Pocetak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.crtez);
ul_sig=get(handles.ucitaj_signal,'userdata');
fs=ul_sig.fs;
x=ul_sig.x;
t=(0:length(x)-1)'/fs;
[tp,nista] = ginput(1);
if tp>0 && tp<t(end)
    izbor_segmenta=get(handles.izbor_segmenta,'value');
    pocetak=get(hObject,'userdata');
    pocetak(izbor_segmenta)=tp;
    set(hObject,'userdata',pocetak);
    kraj=get(handles.Kraj,'userdata');
    if kraj(izbor_segmenta)<=tp
        kraj(izbor_segmenta)=t(end);
        set(handles.Kraj,'userdata',kraj);
    end
    crtaj_signal(handles);
    hold on;
    crtaj_p_k(handles);
end
crtaj_spektar(handles);

% --- Executes on button press in Kraj.
function Kraj_Callback(hObject, eventdata, handles)
% hObject    handle to Kraj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.crtez);
ul_sig=get(handles.ucitaj_signal,'userdata');
fs=ul_sig.fs;
x=ul_sig.x;
t=(0:length(x)-1)'/fs;
[tk,nista] = ginput(1);
if tk>0 && tk<t(end)
    izbor_segmenta=get(handles.izbor_segmenta,'value');
    kraj=get(hObject,'userdata');
    kraj(izbor_segmenta)=tk;
    set(hObject,'userdata',kraj);
    pocetak=get(handles.Pocetak,'userdata');
    % if length(pocetak)<izbor_segmenta
    %     pocetak(izbor_segmenta)=t(1);
    %     set(handles.Pocetak,'userdata',pocetak);
    % end
    if pocetak(izbor_segmenta)>=tk
        pocetak(izbor_segmenta)=t(1);
        set(handles.Pocetak,'userdata',pocetak);
    end
    crtaj_signal(handles);
    hold on;
    crtaj_p_k(handles);
end
crtaj_spektar(handles);

function crtaj_signal(handles)
axes(handles.crtez);
x_za_lim=get(handles.crtez,'xlim');
ul_sig=get(handles.ucitaj_signal,'userdata');
fs=ul_sig.fs;
x=ul_sig.x;
t=(0:length(x)-1)'/fs;
hold off
plot(t,x,'k');
grid on
xlim(x_za_lim);


% --- Executes on button press in X_rezoom.
function X_rezoom_Callback(hObject, eventdata, handles)
% hObject    handle to X_rezoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ul_sig=get(handles.ucitaj_signal,'userdata');
fs=ul_sig.fs;
x=ul_sig.x;
t=(0:length(x)-1)'/fs;
axes(handles.crtez);
x_za_lim=[0 t(end)];
hold off
set(handles.crtez,'xlim',x_za_lim);
crtaj_signal(handles);
hold on;
crtaj_p_k(handles);


% --- Executes on button press in PKX.
function PKX_Callback(hObject, eventdata, handles)
% hObject    handle to PKX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.dodaj_segment,'userdata',1);
izbor_segmenta_str='1';
set(handles.izbor_segmenta,'string',izbor_segmenta_str);
set(handles.izbor_segmenta,'value',1);
set(handles.Pocetak,'userdata',0);
ul_sig=get(handles.ucitaj_signal,'userdata');
fs=ul_sig.fs;
x=ul_sig.x;
t=(0:length(x)-1)'/fs;
k=t(end);
set(handles.Kraj,'userdata',k);
crtaj_signal(handles);
hold on;
crtaj_p_k(handles);
crtaj_spektar(handles);
set(handles.dodaj_segment,'enable','on');

% --- Executes on button press in Smax.
function Smax_Callback(hObject, eventdata, handles)
% hObject    handle to Smax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.spektar);
[fm,nista] = ginput(1);
set(hObject,'userdata',fm);
f_C4=440*2^(-9/12);
fc=fm;
n=round(12*log(fc/f_C4)/log(2));
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
note_ispis=[note_sve{ton+1} num2str(okt)];
set(handles.string_nota,'string',note_ispis);



% --- Executes on button press in dodaj_segment.
function dodaj_segment_Callback(hObject, eventdata, handles)
% hObject    handle to dodaj_segment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ul_sig=get(handles.ucitaj_signal,'userdata');
fs=ul_sig.fs;
x=ul_sig.x;
t=(0:length(x)-1)'/fs;
br_segmenata=get(hObject,'userdata');
br_segmenata=br_segmenata+1;
if br_segmenata>3
    set(hObject,'enable','off');
end
set(hObject,'userdata',br_segmenata);
izbor_segmenta_str=get(handles.izbor_segmenta,'string');
l=length(izbor_segmenta_str);
izbor_segmenta_str(l+1,1)=num2str(br_segmenata);
set(handles.izbor_segmenta,'string',izbor_segmenta_str);
set(handles.izbor_segmenta,'value',l+1);
pocetak=get(handles.Pocetak,'userdata');
if length(pocetak)<br_segmenata
    pocetak(br_segmenata)=t(1);
    set(handles.Pocetak,'userdata',pocetak);
end
kraj=get(handles.Kraj,'userdata');
if length(kraj)<br_segmenata
    kraj(br_segmenata)=t(end);
    set(handles.Kraj,'userdata',kraj);
end
crtaj_signal(handles);
hold on;
crtaj_p_k(handles);


% --- Executes on selection change in izbor_segmenta.
function izbor_segmenta_Callback(hObject, eventdata, handles)
% hObject    handle to izbor_segmenta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns izbor_segmenta contents as cell array
%        contents{get(hObject,'Value')} returns selected item from izbor_segmenta
crtaj_spektar(handles);

% --- Executes during object creation, after setting all properties.
function izbor_segmenta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to izbor_segmenta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in X_spektar.
function X_spektar_Callback(hObject, eventdata, handles)
% hObject    handle to X_spektar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ul_sig=get(handles.ucitaj_signal,'userdata');
fs=ul_sig.fs;
axes(handles.spektar);
x_za_lim=[0 fs/2];
hold off
set(handles.spektar,'xlim',x_za_lim);
crtaj_spektar(handles);



function string_nota_Callback(hObject, eventdata, handles)
% hObject    handle to string_nota (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of string_nota as text
%        str2double(get(hObject,'String')) returns contents of string_nota as a double


% --- Executes during object creation, after setting all properties.
function string_nota_CreateFcn(hObject, eventdata, handles)
% hObject    handle to string_nota (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
da_ne=Poruka_2('title','Warning','string','For furhter analysis the segment 1 and another currently selected segment will be used, if more than one segment is defined.','d1','Proceed.','d2','Return.');
if da_ne(1)=='Y'
    br_sel=0;
    if get(handles.win_tip,'value')==1
        br_sel=br_sel+1;
        prozor=2;
    end
    if get(handles.win_tip,'value')==2
        br_sel=br_sel+1;
        prozor=3;
    end
    br_segmenta=get(handles.izbor_segmenta,'value');
    ul_sig=get(handles.ucitaj_signal,'userdata');
    fs=ul_sig.fs;
    signal_fajl=get(handles.signal_fajl,'userdata');
    p=get(handles.Pocetak,'userdata');
    k=get(handles.Kraj,'userdata');
    if br_segmenta==1 && length(p)>1
        da_ne=Poruka_1('title','Warning','string','For furhter analysis the segment 1 and another currently selected segment will be used, select another segment.','d1','Proceed.','d2','Return.');
    end
    if da_ne(1)=='Y'
        p_p(1)=p(1);
        k_p(1)=k(1);
        if br_segmenta>1
            p_p(2)=p(br_segmenta);
            k_p(2)=k(br_segmenta);
        end
        switch get(handles.Nfft,'value')
            case 1
                Nfft_pom=0;
            case 2
                Nfft_pom=fs;
            case 3
                Nfft_pom=10*fs;
        end
        Tprozor=get(handles.TWIN,'value')/1000;
        save asa_banka signal_fajl p_p k_p prozor Nfft_pom Tprozor
        MuPI_B();
    end
end



function Nfft_labela_Callback(hObject, eventdata, handles)
% hObject    handle to Nfft_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nfft_labela as text
%        str2double(get(hObject,'String')) returns contents of Nfft_labela as a double


% --- Executes on selection change in Nfft.
function Nfft_Callback(hObject, eventdata, handles)
% hObject    handle to Nfft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Nfft contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Nfft
crtaj_spektar(handles);

% --- Executes during object creation, after setting all properties.
function Nfft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nfft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in win_tip.
function win_tip_Callback(hObject, eventdata, handles)
% hObject    handle to win_tip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of win_tip
crtaj_spektar(handles);



function TWIN_labela_Callback(hObject, eventdata, handles)
% hObject    handle to TWIN_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TWIN_labela as text
%        str2double(get(hObject,'String')) returns contents of TWIN_labela as a double


% --- Executes during object creation, after setting all properties.
function TWIN_labela_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TWIN_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TWIN_Callback(hObject, eventdata, handles)
% hObject    handle to TWIN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TWIN as text
%        str2double(get(hObject,'String')) returns contents of TWIN as a double
tmp=get(hObject,'string');
tmpv=str2double(tmp);
set(hObject,'value',tmpv);
crtaj_spektar(handles);

% --- Executes during object creation, after setting all properties.
function TWIN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TWIN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function win_tip_labela_Callback(hObject, eventdata, handles)
% hObject    handle to win_tip_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of win_tip_labela as text
%        str2double(get(hObject,'String')) returns contents of win_tip_labela as a double


% --- Executes during object creation, after setting all properties.
function win_tip_labela_CreateFcn(hObject, eventdata, handles)
% hObject    handle to win_tip_labela (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in eksp_spektar.
function eksp_spektar_Callback(hObject, eventdata, handles)
% hObject    handle to eksp_spektar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fignew = figure('Visible','on'); % Invisible figure
newAxes = copyobj(handles.spektar,fignew); % Copy the appropriate axes
set(newAxes,'Position',get(groot,'DefaultAxesPosition')); % The original position is copied too, so adjust it.


% --- Executes on button press in eksp_vreme.
function eksp_vreme_Callback(hObject, eventdata, handles)
% hObject    handle to eksp_vreme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fignew = figure('Visible','on'); % Invisible figure
newAxes = copyobj(handles.crtez,fignew); % Copy the appropriate axes
set(newAxes,'Position',get(groot,'DefaultAxesPosition')); % The original position is copied too, so adjust it.


% --- Executes on button press in sviraj.
function sviraj_Callback(hObject, eventdata, handles)
% hObject    handle to sviraj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ul_sig=get(handles.ucitaj_signal,'userdata');fs=ul_sig.fs;
x=ul_sig.x;
t=(0:length(x)-1)'/fs;
p=get(handles.Pocetak,'userdata');
k=get(handles.Kraj,'userdata');
izbor_segmenta=get(handles.izbor_segmenta,'value');
xo=x(t>=p(izbor_segmenta) & t<=k(izbor_segmenta));
to=(0:length(x)-1)'/fs;
p=audioplayer(xo,fs);
playblocking(p); 
