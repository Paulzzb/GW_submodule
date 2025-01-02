
%%%%%%%%%%%%%%%%%%%%%%%%
% Check mtxel.f90 all right ? 
% Currently, aqstmp = ksinfo.vol * aqs_prod = n123.^2 * aqs_mtxel
clear all;
close all;
dbstop if error
cd ../../../;
randn('state', 0);
rand('state', 0);
KSSOLV_startup;
cd silicon
si2
cd ../src/GW/testSi2

options = [];
options.isISDF        = false;
options.isCauchy      = false;
options.fileName      = 'hello.mat';
options.frequency_dependence      = 1;
options.nv_ener       = 4;
options.nc_ener       = 4;
options.nv_oper       = mol.nel / 2;
options.nc_oper       = mol.nel / 2;

options = GWOptions(mol, options);
ksinfor = ksinfo(mol, options.Groundstate);
ksinfor2 = load('Siksinfo.mat');
ksinfor2 = ksinfor2.ksinfo;
ksinfor.Z = ksinfor2.Z * sqrt(ksinfor.vol);

nv = options.Groundstate.nv;
nc = options.Groundstate.nc;
sum_range = 1:nv+nc;
for nn = 1:nv+nc
  aqstemp = mtxel_sigma(nn, ksinfor, options.Groundstate, sum_range);
  aqstemp_ = ksinfor2.aqs{nn};
  norm(aqstemp - aqstemp_)
end


