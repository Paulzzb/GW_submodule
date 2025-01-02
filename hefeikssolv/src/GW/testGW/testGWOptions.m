cd ../../../
KSSOLV_startup
cd silicon
si2

options.isGW          = true; 
% options.frequency_dependence      = 0; 
options.frequency_dependence      = 0; 
options.amin          = 5.0; % Currently, this options with unit Ha to consist with the main part of GW method 
options.nv            = mol.nel / 2;
% options.nc            = mol.nel/2;
options.nc            = mol.nel / 2;
options.input         = 'kssolv';
% options.input         = 'qe-gw-bgw';
options.isBSE         = false;
options.isISDF        = true;
options.isCauchy      = true;
% options.optionsISDF
options.vcrank_ratio  = 2; 
options.vsrank_ratio  = 2; 
options.ssrank_ratio  = 2; 
options.fileName      = 'GWoutput.mat';
options.nv_ener       = mol.nel / 2;
options.nc_ener       = mol.nel / 2;
options.nv_oper       = mol.nel / 2;
options.nc_oper       = mol.nel / 2;


cd ../src/GW
options_out = GWOptions(mol, options);
