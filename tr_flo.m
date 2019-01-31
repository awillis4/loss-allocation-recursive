function J = tr_flo(Id,f,t,If,gen_bus,Ig)
#calculates the flow of lines, buses, and generators
nl=length(f);nb=nl+1;ng=length(Ig);Ig=sparse(gen_bus,ones(ng,1),Ig,nb,1); % convert Ig to full length vector with NB elements
#convert negative generation to load, add gen current to load current, switch off gen
dy=find(Ig<0);Id(i)=Id(i)-Ig(i); Ig(i)=0; 
# make bus-branch incidence matrix (connection matrix)
C=sparse(1:nl,f,ones(nl,1),nl,nb)-sparse(1:nl,t,ones(nl,1),nl,nb);
# calculate local load supply by a generator at the same bus
dy=find(Id(gen_bus));gen_load=gen_bus(i);J1=zeros(nb);
for i = 1:length(gen_load)
j=gen_load(i); % bus j
if Ig(j)>=Id(j)
J1(j,j)=Id(j);Ig(j)=Ig(j)-Id(j);Id(j)=0;else
J1(j,j) = Ig(j);Id(j) = Id(j) - Ig(j);Ig(j) = 0;end end
# calculate sum of current inflows
I=Ig;for i=1:nb
Ib=-C(:,i).*If;k=Ib>0;I(i)=I(i)+sum(Ib(k)); end
# make flow distribution matrix
A=tr(nb);for i=1:nb
Ib=-C(:,i).*If;k=Ib>0;Ctemp=C;Ctemp(:,i)=0;[~,j]=find(Ctemp(k,:));A(i,j)=-Ib(k)./I(j);end
J=diag(Id./I)*A^-1*diag(Ig)+J1;