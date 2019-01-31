function print_pf(ds,la)
%PRINT_PF  Writes solution to your data set.

[U,Sg,gen_bus,Sf,St,Sd,Yd,Sbase,f,t]=...
deal(ds.U,ds.Sg,ds.gen_bus,ds.Sf,ds.St,ds.Sd,ds.Yd,ds.Sbase,ds.stem(:,1),ds.stem(:,2));
fid = fopen('loss_allocation.txt','w');
fprintf(fid,'Case file: %s\n\n',ds.file_in);

[Umin,imin]=min(abs(U));[Umax,imax]=max(abs(U));
fprintf(fid,'Umin = %10.6f pu @ %i\n',Umin,imin);
fprintf(fid,'Umax = %10.6f pu @ %i\n',Umax,imax);
[Tmin,imin]=min(angle(U)/pi*180);[Tmax,imax]=max(angle(U)/pi*180);
fprintf(fid,'Angle min = %7.3f deg @ %i\n',Tmin,imin);
fprintf(fid,'Angle max = %7.3f deg @ %i\n\n',Tmax,imax);


fprintf(fid,' Bus Pd(kW) Qd(kvar) Qc(kvar) Pg(kW) Qg(kvar) U(pu) teta(deg)\n');
nb=size(U,1);Sd=Sd*Sbase*1000;Qc=imag(U.*conj(Yd.*U))*Sbase*1000;
Sg=sparse(gen_bus,ones(size(Sg)),Sg,NB,1)*Sbase*1000;Sg=full(Sg);

for i = 1:nb
fprintf(fid,'%3i %11.3f %11.3f %11.3f %11.3f %11.3f %10.6f %10.3f\n',i,dx(Sd(i)),dy(Sd(i)),Qc(i),dx(Sg(i)),dy(Sg(i)),abs(U(i)),angle(U(i))/pi*180);
end
fprintf(fid,'Total %9.3f %11.3f %11.3f %11.3f %11.3f\n',dx(sum(Sd)),dy(sum(Sd)),sum(Qc),dy(sum(Sg)),dy(sum(Sg)));
DS = sum(Sg)-sum(Sd+1j*Qc);
fprintf(fid,'\nDP = %8.5f kW\n',dx(DS));fprintf(fid,'DQ = %8.5f kvar\n\n',dy(DS));

fprintf(fid,'\nFrom  To     Pf(kW)   Qf(kvar)     Pt(kW)   Qt(kvar)     DP(kW)   DQ(kvar)\n');
DSb = (Sf - St)*Sbase*1000;
for i = 1:size(Sf,1)
Sfi = Sf(i)*Sbase*1000;Sti = St(i)*Sbase*1000;
    fprintf(fid,'%4i %3i %10.3f %10.3f %10.3f %10.3f %10.5f %10.5f\n',f(i),t(i),dx(Sfi),dy(Sfi),dx(Sti),dy(Sti),dx(DSb(i)),dy(DSb(i)));
end
fprintf(fid,'%52s %10.5f %10.5f\n','Total',dx(sum(DSb)),sum(dy(DSb)));
fprintf(fid,'\nFrom  To    Re{If}(A)    Im{If}(A)    Re{It}(A)    Im{It}(A)\n');
If = ds.If;It = ds.It;
for i = 1:size(If,1)
    fprintf(fid,'%4i %3i %12.4f %12.4f %12.4f %12.4f\n',f(i),t(i),dx(If(i)),dy(If(i)),dx(It(i)),dy(It(i)));
end
fprintf(fid,'\n Bus    Re{Id}(A)    Im{Id}(A)\n');print_id(fid,1:nb,ds.Id);
fprintf(fid,'\n Bus    Re{Ig}(A)    Im{Ig}(A)\n');print_id(fid,gen_bus,ds.Ig);
% components of load currents
J=ds.J;
fprintf(fid,'\nRe{J}(A)\n');print_j(fid,gen_bus,real(J(:,gen_bus)));
fprintf(fid,'\nIm{J}(A)\n');print_j(fid,gen_bus,imag(J(:,gen_bus)));
fprintf(fid,'\nActive Loss Allocation (kW)\n');print_la(fid,gen_bus,dx(la));
fprintf(fid,'\nReactive Loss Allocation (kvar)\n');print_la(fid,gen_bus,dy(la));

LAsum=sum(sum(la));
fprintf(fid,'\nTotal Allocated Losses\n');
fprintf(fid,'DP = %8.5f kW\n',dx(LAsum));
fprintf(fid,'DQ = %8.5f kvar\n\n',dy(LAsum));
fprintf(fid,'DPdiff = %8.5f kW (%.2f %%)\n',dx(LAsum)-dx(DS),(dx(LAsum)/dx(DS)-1)*100);
fprintf(fid,'DQdiff = %8.5f kvar (%.2f %%)\n\n',dy(LAsum)-dy(DS),(dy(LAsum)/dy(DS)-1)*100);

fclose(fid);

function print_id(fid,bus,Id)
for i = 1:length(Id)
    fprintf(fid,'%4i %12.4f %12.4f\n',bus(i),dx(Id(i)),dy(Id(i)));end

function print_j(fid,gen_bus,J)
[m,n]=size(J);
fprintf(fid,'    ');
for i = 1:n
    fprintf(fid,'%12i',gen_bus(i));end
fprintf(fid,'\n');
for i = 2:m
    fprintf(fid,'%4i',i);
    for j = 1:n
        fprintf(fid,'%12.4f',J(i,j));end
    fprintf(fid,'\n');end

function print_la(fid,gen_bus,LA)
fprintf(fid,'Bus/Gen      1');
for i = 2:length(gen_bus)
fprintf(fid,'%11i',gen_bus(i));end
fprintf(fid,'\n');
for i = 2:length(LA)
fprintf(fid,'%3i',i);
for j = 1:length(gen_bus)
fprintf(fid,'%11.4f',LA(i,j));end
fprintf(fid,'\n');end