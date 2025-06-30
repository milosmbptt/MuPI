function N=emqf_ord_HB(wp,As)
Ap=10*log10(1+1/(10^(As/10)-1));
Kp=sqrt(1-10^(-Ap/10))/(10^(-Ap/20));
Ks=sqrt(1-10^(-As/10))/(10^(-As/20));
ksi=1/tan(wp/2)^2;
L=Ks/Kp;
Num=ellipke(1-1/L^2)/ellipke(1/L^2);
Den=ellipke(1-1/ksi^2)/ellipke(1/ksi^2);
N=round(Num/Den);
if rem(N,2)==0
    N=N+1;
end