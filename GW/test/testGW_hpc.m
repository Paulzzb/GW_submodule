flag = 0;

testfolder = '/public/home/jlyang/mahh/hefeikssolv/src/GW/Si8';
if flag == 1
  % clear all;
  % close all;
  dbstop if error
  cd ../../../
  KSSOLV_startup;
  % cd silicon
  % si96
  % si32
  
  cd ./src/GW/
  cd(testfolder);

  load("molinfo.mat");
  load("scinfo.mat");
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % We gives a most simple choice for options, check $home/src/GW/README.md 
  % for details.
  options = [];
  % options.isISDF        = false;
  options.isISDF        = true;
  options.isCauchy      = true;
  options.fileName      = 'GWoutput.mat';
  options.optionsISDF = [];
  options.optionsISDF.exxmethod = 'kmeans';
  options.optionsISDF.weight = 'prod';
  options.inputfile = "scinfo.mat";
  mol.nel = mol.nel / 2;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  options = GWOptions(mol, options);
  ksinfor = ksinfo(mol, options.Groundstate);

  save("kstest.mat", "ksinfor", '-v7.3');
  save("optest.mat", "options");
else
  cd ../../../
  KSSOLV_startup;
  cd src/GW
  cd(testfolder);
  load("kstest.mat");
  load("optest.mat");
end
flag = 0;

if 0
profile on
gw_cohsex_gpu(ksinfor, options);
profile off
profsave;
load GWoutput.mat
GWenergynew = GWenergy;
gw_cohsex(ksinfor, options);
load GWoutput.mat
else
  % do the test to ISDF_kernelg
  gw_testISDF(ksinfor, options);
end

% options.isISDF        = true;
% options.fileName      = 'GWoutputSi8.mat';
% gwCalculation(mol, options);
% 
% options.isCauchy      = true;
% options.fileName      = 'GWoutputSi8.mat';
% gwCalculation(mol, options);
