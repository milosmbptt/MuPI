function bete=bmqf_HB(N)
[z,p,k]=butter(N,0.5);
beta=abs(p).^2;
bete=sort(beta);
bete=bete(2:2:end);