clear;
%
INFO_DIR = 'Si8/';

cfile = mfilename('fullpath');
CPATH = fileparts(cfile);
CPATH = [CPATH, '/'];
KS_DIR = ['D:/kssolvGW/hefeikssolv/'];
GW_DIR = ['D:/kssolvGW/GW/'];
GWsetfile = [CPATH, INFO_DIR, 'setGW.m'];
INTER_DIR = [CPATH, INFO_DIR, 'IntermediateFiles/'];
GWinput = [INTER_DIR, 'scinfo.mat'];
 

% 1. Managing paths.
cd(KS_DIR);
KSSOLV_startup;
cd(CPATH);
addpath(GW_DIR);

% 2. Load
load(GWinput);

% 3. run
run(GWsetfile);
options_GW.inputfile = GWinput;
gwCalculation(mol, options_GW);


