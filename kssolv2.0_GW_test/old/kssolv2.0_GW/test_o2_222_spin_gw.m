tic
% Cleaning up the workspace and command line window:
clc
clear all;
close all;
randn('state', 0);
rand('state', 0);
% Initializing the environment path for KSSOLV:
KSSOLV_startup;

[sys, options, syms] = read_qe_gw_bgw('.\gw\qe_data\o2_222_spin');
[sys, options] = gw_setup(sys, options);

eps.nbnd = 15;
eps.nv = options.nv;
eps.nc = eps.nbnd - eps.nv;
omega = 0;
eta = 0;
eps.cutoff = 2;
eps.coul_cutoff = 2; %coulomb truncation radius in epsilon
eps = epsilon(sys, options, syms, eps);

sig.nbnd = 15;
sig.ndiag_max = 15;
sig.coul_cutoff = 2; %coulomb truncation radius in sigma
sig.no_symmetries_q_grid = 0;
sig.exact_static_ch = 1;
sig = sigma(eps, sig, sys, options, syms);

toc
save('test_02.mat')