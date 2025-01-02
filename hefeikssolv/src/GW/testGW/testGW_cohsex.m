clear all;
close all;
dbstop if error
cd ../../../
KSSOLV_startup;
cd silicon
si16
cd ../src/GW/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We gives a most simple choice for options, check $home/src/GW/README.md 
% for details.
options = [];
options.isISDF        = false;
options.isCauchy      = false;
options.fileName      = 'GWoutputSi16_NoISDF.mat';
options.nv_ener       = 4;
options.nc_ener       = 4;
options.nv_oper       = mol.nel / 2;
options.nc_oper       = mol.nel / 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gwCalculation(mol, options);

options.isISDF        = true;
options.fileName      = 'GWoutputSi16_ISDF.mat';
gwCalculation(mol, options);

options.isCauchy      = true;
options.fileName      = 'GWoutputSi16_ISDFCauchy.mat';
gwCalculation(mol, options);
