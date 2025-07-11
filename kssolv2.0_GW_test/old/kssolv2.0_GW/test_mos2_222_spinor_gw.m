% Cleaning up the workspace and command line window:
clc
clear all;
close all;
randn('state', 0);
rand('state', 0);
% Initializing the environment path for KSSOLV:
KSSOLV_startup;

[sys, options, syms] = read_qe_gw_bgw('.\gw\qe_data\mos2_222_spinor');
[sys, options] = gw_setup(sys, options);

eps.nbnd = 29;
eps.nv = options.nv;
eps.nc = eps.nbnd - eps.nv;
omega = 0;
eta = 0;
eps.cutoff = 10; %unit is Ry
eps.coul_cutoff = 10; %coulomb truncation radius in epsilon
eps = epsilon(sys, options, syms, eps);

ndiag_max = 29;
coul_cutoff = 10; %coulomb truncation radius in sigma
sigma_inp.no_symmetries_q_grid = 1;
[cor, sig] = sigma(ndiag_max, sys.nbnd, epsmat, coul_cutoff, sys, options, syms, sigma_inp);
save('test')