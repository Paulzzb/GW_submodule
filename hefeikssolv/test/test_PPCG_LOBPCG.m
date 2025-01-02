clear all;
close all;
randn('state', 0);
rand('state', 0);
global version;
version = switchver('CPU');
cd ..
KSSOLV_startup;
cd example;
sih4_setup;
%si8_setup;
%si64_setup;
%options=setksopt('eigmethod','lobpcg');
options=setksopt('eigmethod','ppcg');
[mol,H,X,info] = scf(mol);
cd ..
