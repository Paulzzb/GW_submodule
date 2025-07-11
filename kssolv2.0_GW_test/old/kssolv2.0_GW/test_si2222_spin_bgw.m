% Cleaning up the workspace and command line window:
clc
clear all;
close all;
randn('state', 0);
rand('state', 0);
% Initializing the environment path for KSSOLV:
KSSOLV_startup;

[sys, options, syms] = read_qe_gw_bgw('.\gw\qe_data\spin_si2_222');
[sys, options] = gw_setup(sys, options);

eps.nbnd = 10;
eps.nv = options.nv;
eps.nc = eps.nbnd - eps.nv;
omega = 0;
eta = 0;
eps.cutoff = 10;
eps.coul_cutoff = 5; %coulomb truncation radius in epsilon
eps = epsilon(sys, options, syms, eps);

sig.nbnd = 10;
sig.ndiag_max = 10;
sig.coul_cutoff = 5; %coulomb truncation radius in sigma
sig.no_symmetries_q_grid = 0;
sig.exact_static_ch = 1;
sig = sigma(eps, sig, sys, options, syms);