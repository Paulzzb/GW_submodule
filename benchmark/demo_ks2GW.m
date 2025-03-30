clear all;
close all;
%
INFO_DIR = 'Si8/';

cfile = mfilename('fullpath');
CPATH = fileparts(cfile);
CPATH = [CPATH, '/'];
KS_DIR = ['D:/kssolvGW/GW_submodule/hefeikssolv/'];
GW_DIR = ['D:/kssolvGW/GW_submodule/GW/'];
GWsetfile = [CPATH, INFO_DIR, 'setGW_fullfreq_cd.m'];
INTER_DIR = [CPATH, INFO_DIR, 'IntermediateFiles/'];
GWinput = [INTER_DIR, 'scinfo.mat'];
 

% 1. Managing paths.
cd(KS_DIR);
KSSOLV_startup;
cd(GW_DIR);
gw_startup;
cd(CPATH);

% 2. Load
load(GWinput);

% 3. run
run(GWsetfile);
options_GW.inputfile = GWinput;
[GWinput, optionsGW] = kssolv2GW(options_GW, mol);

% 4. Save
save([INTER_DIR, 'GWinput.mat'], 'GWinput', 'optionsGW');



