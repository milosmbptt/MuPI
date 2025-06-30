function bete=emqf_HB(N,wp)
ksi=1/tan(wp/2)^2;
for br=1:(N-1)/2
    x(br)=ellipj(((2*br-1)/N+1)*ellipke(1/ksi^2),1/ksi^2);
    beta(br)=(ksi+x(br)^2-sqrt((1-x(br)^2)*(ksi^2-x(br)^2)))/(ksi+x(br)^2+sqrt((1-x(br)^2)*(ksi^2-x(br)^2)));
end
bete=sort(beta);
