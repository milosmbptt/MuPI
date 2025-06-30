function izlaz=fazna_korekcija(wc_uk,fp,as,ulaz,b_ili_e)
izlaz=zeros(size(ulaz));

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



for br_sig=1:size(ulaz,2)
    ulaz_tmp=ulaz(end:-1:1,br_sig);
	for br=1:length(wc_uk)
        for br_A0=1:br_sec_0
            im=[1 alfa(br)*(1+beta0(br,br_A0)) beta0(br,br_A0)];
            izlaz_tmp=filter(fliplr(im),im,ulaz_tmp);
	        ulaz_tmp=izlaz_tmp;
        end
    end
    izlaz(:,br_sig)=izlaz_tmp(end:-1:1);
end