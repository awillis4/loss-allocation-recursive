function Zpv=make_zpv(pv,nb,nl,f,Zstem)
%MAKE_ZPV  Calculates loop impedances for all PV buses.
npv=length(pv);Zpv=zeros(npv);
for ipv=1:npv
V=zeros(nb,1);I=zeros(nb,1);I(pv(ipv))=-1;
for k=nl:-1:2
i=f(k);I(i)=I(i)+I(k);end
for k=2:nl
i=f(k);V(k)=V(i)-Zbranch(k)*I(k);end 
Zpv(:,ipv)=V(pv);end