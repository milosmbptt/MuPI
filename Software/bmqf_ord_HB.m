function N=bmqf_ord_HB(wp,As)
Ap=10*log10(1+1/(10^(As/10)-1));
Kp=sqrt(1-10^(-Ap/10))/(10^(-Ap/20));
Ks=sqrt(1-10^(-As/10))/(10^(-As/20));
ksi=1/tan(wp/2)^2;
L=Ks/Kp;
Num=log(L);
Den=log(ksi);
N=round(Num/Den);
if rem(N,2)==0
    N=N+1;
end