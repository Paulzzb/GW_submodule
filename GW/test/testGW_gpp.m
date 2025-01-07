clear all;
close all;
dbstop if error
cd ../../../
KSSOLV_startup;
cd silicon
si2
cd ../src/GW/testSi2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We gives a most simple choice for options, check $home/src/GW/README.md 
% for details.
options = [];
options.isISDF        = false;
options.isCauchy      = false;
options.fileName      = 'GWoutputSi2gpp.mat';
options.frequency_dependence      = 1;
options.nv_ener       = 4;
options.nc_ener       = 4;
options.nv_oper       = mol.nel / 2;
options.nc_oper       = mol.nel / 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = GWOptions(mol, options);
% frequency_dependence = 0
ksinfor = ksinfo(mol, options.Groundstate);

gw_gpp(ksinfor, options);
