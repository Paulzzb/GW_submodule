clear all;
close all;
dbstop if error
cd ../../../
KSSOLV_startup;
cd silicon
si2
cd ../src/GW/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We gives a most simple choice for options, check $home/src/GW/README.md 
% for details.
options = [];
options.isISDF        = false;
options.isCauchy      = false;
options.fileName      = 'GWoutputSi16_NoISDF.mat';
options.frequency_dependence      = 1;
options.nv_ener       = 4;
options.nc_ener       = 4;
options.nv_oper       = mol.nel / 2;
options.nc_oper       = mol.nel / 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = GWOptions(mol, options);
% frequency_dependence = 0
ksinfor = ksinfo(mol, options.Groundstate);

for igcol = 1:10
  out = wpeff(ksinfor, options, igcol);
end

