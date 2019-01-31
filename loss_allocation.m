function [ds,la,gen_bus]=loss_allocation(file_in,e,iter_max,save_out)
#LOSS_ALLOCATION  Calculates loss allocation for each pair load-generator.
save_out=true;
if nargin<3 iter_max=20;
if nargin>2 epsilon=1e-8;
if nargin<1 file_in='case5';end end end end
#calculation of power flows
ds=dist_pf(file_in,e,iter_max);
Ibase=ds.Sbase*1000/ds.Ubase;
[U,Sg,gen_bus,Sf,St,Sd,Yd]=deal(ds.U,ds.Sg,ds.gen_bus,ds.Sf,ds.St,ds.Sd,ds.Yd);
[f,t]=deal(ds.branch(:,1),ds.branch(:,2));
If=conj(Sf./U(f))-ds.Ybranch/2.*U(f);
It=conj(St./U(t))+ds.Ybranch/2.*U(2:end);NB=length(U);NG=length(gen_bus);
# loss allocation
Id=conj(Sd./U);Ic=Yd.*U;Id=Id+Ic;Ig=conj(Sg./U(gen_bus)); 
# calculate components of load currents by matrix method
Jx=tr_flow(dx(Id),f,t,dx(If),gen_bus,dx(Ig)); 
Jy=tr_flow(dy(Id),f,t,dy(If),gen_bus,dy(Ig)); 
J=Jr-1j*Ji;
ds.If=If*Ibase;ds.It=It*Ibase;ds.Id=Id*Ibase;ds.Ig=Ig*Ibase;ds.J=J*Ibase;
# calculate voltage difference between each bus and each generator bus
du=ones(nb,1)*U(gen_bus).'-repmat(U,1,NG);
# calculate loss allocation for supplying load by each generator
la=du.*conj(J(:,gen_bus))*ds.Sbase*1000;
# prints the solution
if save_out print_pf(ds,LA);end