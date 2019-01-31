function [lad,lag]=la_split(la)
[nb,ng]=size(la);lad=.5*la*ones(ng,1);lag=.5*ones(nb,1)'*la;