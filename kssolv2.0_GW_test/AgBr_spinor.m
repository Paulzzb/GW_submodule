tic
% Cleaning up the workspace and command line window:
clc
clear all;
close all;
randn('state', 0);
rand('state', 0);

% 启动 Profiler
profile on

% Initializing the environment path for KSSOLV:
KSSOLV_startup;

[sys, options, syms] = read_qe_gw_bgw('./gw/qe_data/AgBr');
[sys, options] = gw_setup(sys, options);

eps.nbnd = 99;
eps.nv = options.nv;
eps.nc = eps.nbnd - eps.nv;
omega = 0;
eta = 0;
eps.cutoff = 5; % Ry
eps.coul_cutoff = 5; %coulomb truncation radius in epsilon
eps = epsilon(sys, options, syms, eps);

sig.nbnd = 99;
sig.ndiag_max = 90;
sig.coul_cutoff = 5; %coulomb truncation radius in sigma
sig.no_symmetries_q_grid = 0;
sig.exact_static_ch = 1;
sig = sigma(eps, sig, sys, options, syms);

sig_nv = -sig.eqp0(options.nv, :)
toc

% 停止 Profiler 并生成报告
profile off
profile report

% 保存 Profiler 数据
profinfo = profile('info');
save('profiler_results.mat', 'profinfo');

save('AgBr.mat','-v7.3')