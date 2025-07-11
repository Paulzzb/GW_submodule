% Cleaning up the workspace and command line window:
clc
clear all;
close all;
randn('state', 0);
rand('state', 0);
% Initializing the environment path for KSSOLV:
KSSOLV_startup;
% Initializing the structure of a crystal or a molecule:
cd examples;
sih4_setup;

% If use hybrid functional to calculate:
%ppvar.funct = 'HSE06';
%mol.ppvar = ppvar;
%mol.funct = ppvar.funct;

% Setting the calculation parameters by default:
options = setksopt();
%options.useace = 1;

% Opening the profile function of MATLAB:
profile on;
[mol,H,X,info] = scf(mol,options);
profile viewer;

% Plotting the SCF convergence curve: 
semilogy(info.SCFerrvec,'-.o');
xlabel('SCF iteration number');                                           
ylabel('SCF error')

% Plotting the electron density:
view(3);
isosurface(fftshift(H.rho));
view(3);
vol3d('cdata',fftshift(H.rho));
cd ..
