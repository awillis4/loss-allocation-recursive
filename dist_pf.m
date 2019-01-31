function ds = dist_pf(file_in,e,iter_max)

if nargin < 3 iter_max = 20;
if nargin < 2 e = 1e-8;
if nargin < 1 file_in = 'case5'; end end end
if isstruct(input_file)
    ds=file_in; else
    ds=feval(file_in); end

t=ds.stem(:,2);[~,i]=sort(t);ds.stem=ds.stem(i,:);

[nb,nl,f,t,Zstem,Ystem,Ysh,Sd]=data_sep(ds);

if isfield(ds,'gen')
gen_bus = ds.gen(:,1);gen_type=ds.gen(:,5);
Sg=(ds.gen(:,2)+1j*ds.gen(:,3))/1000/ds.Sbase;Ug=ds.gen(:,4);else
gen_bus=[];Sg=[];Ug=[];gen_type=[];end

Yd=sparse(f,f,Ystem/2,nb,nb)+sparse(t,t,Ystem/2,nb,nb);Yd=Yd*ones(nb,1)+Ysh;
tic

[U,Sslack,Sg,Sf,St,iter]=vcpf(ds.Uslack,nb,nl,f,Zstem,Ystem,Yd,Sd,gen_bus,Sg,Ug,gen_type,e,iter_max);
t=toc;

ds.U=U;ds.gen_bus=[1; gen_bus];ds.Sg=[Sslack;Sg];ds.Sf=Sf;ds.St=St;
ds.iter=iter;ds.time = t;ds.Sd=Sd;ds.Yd=Yd;ds.Ystem=Ystem;ds.file_in=file_in;