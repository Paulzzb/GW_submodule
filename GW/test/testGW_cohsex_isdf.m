clear all;
close all;
%
INFO_DIR = 'data/Si8/';

CPATH = mfilename('fullpath');
CPATH = fileparts(CPATH);
[GWPATH, ~, ~] = fileparts(CPATH);
% cfile = mfilename('fullpath');
% CPATH = fileparts(cfile);
CPATH = [CPATH, '/'];
KS_DIR = 'D:/kssolvGW/hefeikssolv/';
GWsetfile = [CPATH, INFO_DIR, 'setGW_cohsex_isdf.m'];
INTER_DIR = [CPATH, INFO_DIR, 'IntermediateFiles/'];
scoutput = [INTER_DIR, 'scinfo.mat'];
 

% 1. Managing paths.
cd(KS_DIR);
KSSOLV_startup;
cd(GWPATH);
gw_startup();
cd(CPATH);

% 2. Load
load(scoutput);

% 3. run
run(GWsetfile);
optionsGW.inputfile = scoutput;
[GWinput, optionsGW] = kssolv2GW(mol, optionsGW);
gwCalculation(GWinput, optionsGW);


