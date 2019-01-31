function V_cor=make_vcor(dd,pv,nb,nl,f,Zb)

V_corr=zeros(nb,1);I=zeros(nb,1);I(pv)=dd;
#back sweep
for k = nl:-1:2
i=f(k);I(i) = I(i) + I(k);end
#fore sweep
for k = 2:nl
i=f(k);V_cor(k)=V_cor(i)-Zb(k)*I(k);end