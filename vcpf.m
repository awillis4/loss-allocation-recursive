function [U,Sslack,Sg,Sf,St,iter]=vcpf(Uslack,nb,nl,f,Zstem,Ystem,Yd,Sd,gen_bus,Sg,Ug,gen_type,e,iter_max)

U=Uslack*ones(nb,1);Uold=U;iter=0;finish=0;
#stems point to nodes
f=[0;f];Zstem=[0;Zstem];nl=nl+1;
# Split generators into 2 types
gen_pq=find(gen_type==1);gen_pv=find(gen_type==2);pv=gen_bus(gen_pv);Ug=Ug(gen_pv);
# For PQ gens add negative of their dy power to bus load
Sd(gen_bus(gen_pq))=Sd(gen_bus(gen_pq))-Sg(gen_pq);
# For PV gens add negative of their dx power to bus load
Sd(gen_bus(gen_pv))=Sd(gen_bus(gen_pv))-dx(Sg(gen_pv));
# Make impedance matrix for all PV buses
Zpv=make_zpv(pv,nb,nl,f,Zstem);
if size(Zpv,1) > 0
  Bpv=(dy(Zpv))^-1;end
# Voltage Correction of Power Flows
npv=length(pv);Qpv=zeros(npv,1);
#Initial branch flows are equal to node first load
while finish==0&&iter<iter_max
iter=iter+1;S=Sd+conj(Yd).*abs(U).^2;St=S;Sf=St;
#Back sweep
for k=NL:-1:2
i=f(k);Sf(k)=St(k)+Zstem(k)*abs(St(k)/U(k))^2;St(i)=St(i)+Sf(k);end
#Fore sweep
for k=2:NL
i=f(k);U(k)=U(i)-Zbranch(k)*conj(Sf(k)/U(i));end
#Checks for convergence
    DU=abs(U-Uold);
if max(D)>e
Uold = U;
if ~isempty(pv)
            % Calculate reactive power correction for PV generators
DE=(Ug./abs(U(pv))-1).*dx(U(pv));DD=Bpv*DE;DC=DD.*dy(U(pv))./dx(U(pv));
            % Make voltage correction
V_cor=make_vcor(DC+1j*DD,pv,nb,nl,f,Zstem);U=U+V_cor;DQ=DD.*abs(U(pv)).^2 ./real(U(pv));
            % Update reactive power for PV generators
            Qpv = Qpv + DQ;Sd(pv) = Sd(pv) - 1j*DQ;end else
												finish = 1;end end

% Save generators reactive powers and branch flows
if ~isempty(gen_pv)
Sg(gen_pv)=dx(Sg(gen_pv))+1j*Qpv;end
Sslack=St(1);Sf=Sf(2:end);St=St(2:end);f=f(2:end);
% Account for branch shunt power flows
Sf=Sf+conj(Ystem).*abs(U(f)).^2/2;St=St-conj(Ystem).*abs(U(2:end)).^2/2;