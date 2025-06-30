function izlaz=banka_all_pass(wc_uk,fp,as,signal,crt,fsamp,b_ili_e,redosled)

if b_ili_e==1
    N=emqf_ord_HB(pi*fp,as);
else
    N=bmqf_ord_HB(pi*fp,as);
end

if b_ili_e==1
    beta=emqf_HB(N,pi*fp);
else
    beta=bmqf_HB(N);
end

[beta_HB,prazno]=sort(beta,'ascend');
beta0_HB=beta_HB(1:2:length(beta_HB));
beta1_HB=beta_HB(2:2:length(beta_HB));

br_sec_0=length(beta0_HB);
br_sec_1=length(beta1_HB);

alfa=zeros(1,length(wc_uk));
alfa1=zeros(1,length(wc_uk));
beta0=zeros(length(wc_uk),length(beta0_HB));
beta1=zeros(length(wc_uk),length(beta1_HB));

for br=1:length(wc_uk)
    wc=wc_uk(br);
    alfa(br)=-cos(wc);
    alfa1(br)=(1-sqrt(1-alfa(br)^2))/alfa(br);
    beta0(br,:)=(beta0_HB+alfa1(br)^2)./(beta0_HB*alfa1(br)^2+1);
    beta1(br,:)=(beta1_HB+alfa1(br)^2)./(beta1_HB*alfa1(br)^2+1);
end


ulaz=zeros(length(signal),length(wc_uk)+1);
izlaz=zeros(length(signal),length(wc_uk)+1);
ulaz(:,1)=signal;
if redosled==0
    for br=1:length(wc_uk)
        for br_pom=1:br-1
            ulaz_tmp=izlaz(:,br_pom);
            for br_A0=1:br_sec_0
                im=[1 alfa(br)*(1+beta0(br,br_A0)) beta0(br,br_A0)];
                izlaz_tmp=filter(fliplr(im),im,ulaz_tmp);
                ulaz_tmp=izlaz_tmp;
            end
            izlaz(:,br_pom)=ulaz_tmp;
        end
        ulaz_tmp=ulaz(:,br);
        for br_A0=1:br_sec_0
            im=[1 alfa(br)*(1+beta0(br,br_A0)) beta0(br,br_A0)];
            izlaz_tmp=filter(fliplr(im),im,ulaz_tmp);
            ulaz_tmp=izlaz_tmp;
        end
        izlaz_A0=ulaz_tmp;
        im=[1 alfa1(br)];
        ulaz_tmp=filter(fliplr(im),im,ulaz(:,br));
        for br_A1=1:br_sec_1
            im=[1 alfa(br)*(1+beta1(br,br_A1)) beta1(br,br_A1)];
            izlaz_tmp=filter(fliplr(im),im,ulaz_tmp);
            ulaz_tmp=izlaz_tmp;
        end
        izlaz_A1=ulaz_tmp;
        izlaz(:,br)=(izlaz_A0+izlaz_A1)/2;
        ulaz(:,br+1)=(izlaz_A0-izlaz_A1)/2;
    end
    izlaz(:,end)=ulaz(:,end);
else
    for br=1:length(wc_uk)
        for br_pom=1:br-1
            ulaz_tmp=izlaz(:,length(wc_uk)-br_pom+2);
            for br_A0=1:br_sec_0
                im=[1 alfa(length(wc_uk)-br+1)*(1+beta0(length(wc_uk)-br+1,br_A0)) beta0(length(wc_uk)-br+1,br_A0)];
                izlaz_tmp=filter(fliplr(im),im,ulaz_tmp);
                ulaz_tmp=izlaz_tmp;
            end
            izlaz(:,length(wc_uk)-br_pom+2)=ulaz_tmp;
        end
        ulaz_tmp=ulaz(:,br);
        for br_A0=1:br_sec_0
            im=[1 alfa(length(wc_uk)-br+1)*(1+beta0(length(wc_uk)-br+1,br_A0)) beta0(length(wc_uk)-br+1,br_A0)];
            izlaz_tmp=filter(fliplr(im),im,ulaz_tmp);
            ulaz_tmp=izlaz_tmp;
        end
        izlaz_A0=ulaz_tmp;
        im=[1 alfa1(length(wc_uk)-br+1)];
        ulaz_tmp=filter(fliplr(im),im,ulaz(:,br));
        for br_A1=1:br_sec_1
            im=[1 alfa(length(wc_uk)-br+1)*(1+beta1(length(wc_uk)-br+1,br_A1)) beta1(length(wc_uk)-br+1,br_A1)];
            izlaz_tmp=filter(fliplr(im),im,ulaz_tmp);
            ulaz_tmp=izlaz_tmp;
        end
        izlaz_A1=ulaz_tmp;
        izlaz(:,length(wc_uk)-br+2)=(izlaz_A0-izlaz_A1)/2;
        ulaz(:,br+1)=(izlaz_A0+izlaz_A1)/2;
    end
    izlaz(:,1)=ulaz(:,end);    
end



if crt==1
    IZLAZ=fft(izlaz);
    IZLAZ_sum=fft(izlaz);
    IZLAZ_uk=sum(IZLAZ,2);
    IZLAZ_pow=sum(IZLAZ.*conj(IZLAZ),2);
    f=(0:length(IZLAZ)-1)/length(IZLAZ)*fsamp;
    figure,plot(f/1000,20*log10(abs(IZLAZ)),f,20*log10(abs(IZLAZ_uk)),'k',f,20*log10(abs(IZLAZ_pow)),'k:'),xlim([0 fsamp/2/1000]);
    ylim([-120 10]);
    xlabel('f [kHz]');
    ylabel('Gain [dB]');
    figure,semilogx(f,20*log10(abs(IZLAZ)),f,20*log10(abs(IZLAZ_uk)),'k',f,20*log10(abs(IZLAZ_pow)),'k:'),xlim([20 fsamp/2]);
    ylim([-120 10]);
    xlabel('f [kHz]');
    ylabel('Gain [dB]');    
    figure,plot(f,unwrap(angle(IZLAZ)),f,unwrap(angle(IZLAZ_uk)),'k:'),xlim([0 fsamp/2]);
end