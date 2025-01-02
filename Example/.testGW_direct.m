
cd LiH/222
LiH_scf
cd ../../

flag = 1;
testfolder = '/public/home/zzb/hefeikssolvToPub/test/LiH/222';
if flag == 1
  % clear all;
  % close all;
  dbstop if error
  cd ../
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
  % options.nv_ener       = 4;
  % options.nc_ener       = 4;
  % options.nv_oper       = mol.nel / 2;
  % options.nc_oper       = mol.nel / 2;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  fprintf("Start initialize data for GW calculation.\n");
  options = GWOptions(mol, options);
  % options.ISDFCauchy.flagsunit2super = 1
  ksinfor = ksinfo(mol, options.Groundstate);
  fprintf('Start saving prepared data for GW calculation.\n')

  save("kstest.mat", "ksinfor", '-v7.3');
  save("optest.mat", "options");
  fprintf('Data for GW calculation is saved.\n');
else
  cd ../../
  KSSOLV_startup;
  cd src/GW
  cd(testfolder);
  load("kstest.mat");
  load("optest.mat");
end
