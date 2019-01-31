function [nb,nl,f,t,Zstem,Ystem,Ysh,Sd]=data_sep(ds)

[f,t,Zstem,Ystem,Sd,Qc] = ...
deal(ds.stem(:,1),ds.stem(:,2),ds.stem(:,8),ds.stem(:,3)+1j*ds.stem(:,4), ...
1j*ds.stem(:,5)/1e6,ds.stem(:,6)+1j*ds.stem(:,7));

nl=size(ds.stem,1);nb=max([f;t]);Zbase=ds.Ubase^2/ds.Sbase;
Sd=[0;Sd]/1000/ds.Sbase;Ysh=[0;1j*Qc]/1000/ds.Sbase;
Zstem=Zstem/Zbase;stem=Ystem*Zbase;