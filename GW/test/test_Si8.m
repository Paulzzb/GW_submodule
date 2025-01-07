flag = 0;

testfolder = '/public/home/jlyang/mahh/hefeikssolv/src/GW/Si8ex';

dbstop if error
cd ../../../
KSSOLV_startup;
cd ./src/GW/
cd(testfolder);
load("molinfo.mat");
load("scinfo.mat");
options_input.inputfile = "scinfo.mat";

options_input.isGW          = true; 
options_input.frequency_dependence      = 2; 
options_input.amin          = 5.0; % Currently, this options_input with unit Ha to consist with the main part of GW method 
options_input.nv            = mol.nel/2;
options_input.nc            = mol.nbnd-options_input.nv-5;
options_input.input         = 'kssolv';
options_input.isbse         = false;
options_input.isISDF        = true;
options_input.iscauchy      = true;
options_input.vcrank_ratio  = 16; 
options_input.vsrank_ratio  = 16; 
options_input.ssrank_ratio  = 16; 
options.optionsISDF = [];
options.optionsISDF.exxmethod = 'kmeans';
options.optionsISDF.weight = 'prod';
options.inputfile = "scinfo.mat";

options_input.fileName      = 'GWoutput.mat';
options_input.nv_ener       = mol.nel/2;
options_input.nc_ener       = mol.nbnd-options_input.nv-5;
options_input.nv_oper       = mol.nel/2;
options_input.nc_oper       = mol.nbnd-options_input.nv-5;

options_input.frequency_dependence_method= 2;
options_input.setting_method = 'kssolv';

options_input.delta_frequency = 2.0;
options_input.frequency_low_cutoff = 40;
options_input.number_imaginary_freqs = 20;   %default number of frequency imaginary part



options = GWOptions(mol, options_input);
ksinfor = ksinfo(mol, options.Groundstate);

save("kstest.mat", "ksinfor", '-v7.3');
save("optest.mat", "options");

gw_x(ksinfor, options);


