function [koef_B]=procena_B(B_pocetno,fk1,br_harmonika)
k=1;
delta=1;
B=B_pocetno;
while(k<31)
    fk_procenjeno(k,:)=fk1(1)*(1:br_harmonika).*sqrt(1+B(k)*(1:br_harmonika).^2);
    D=(fk1'-fk_procenjeno(k,1:br_harmonika));
    znak(k)=sign(sum(sign(diff(D))));

    if (k>=2) && (sign(znak(k))~= sign(znak(k-1)))
        delta(k+1)=delta(k)/2;
    else
        delta(k+1)=delta(k);
    end

    if znak(k)>=0
        B(k+1)=B(k)*10^(delta(k));
    else
        B(k+1)=B(k)*10^(-delta(k));
    end

    k=k+1;
end
koef_B=B(end);
