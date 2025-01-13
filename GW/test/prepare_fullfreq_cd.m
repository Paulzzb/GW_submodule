clear all;
close all;
restoredefaultpath;

INFO_DIR = 'data/Si8/';
KS_DIR = 'D:/kssolvGW/hefeikssolv/';

CPATH = mfilename('fullpath');
CPATH = fileparts(CPATH);
[GWPATH, ~, ~] = fileparts(CPATH);
CPATH = [CPATH, '/'];
INTER_DIR = [CPATH, INFO_DIR, 'IntermediateFiles/'];
 

% Output from kssolv--scf
scoutput = [INTER_DIR, 'scinfo.mat'];
% Some parameters for GW_cohsex
GWsetfile = [CPATH, INFO_DIR, 'setGW_fullfreq_cd.m'];
% Save GW input information on this file
GWinputfile = [INTER_DIR, 'GWinput_fullfreq_cd.mat'];

cd(KS_DIR); KSSOLV_startup; cd(CPATH);
cd(GWPATH); gw_startup; cd(CPATH);

load(scoutput)
run(GWsetfile)
optionsGW.inputfile = scoutput;
[GWinput, optionsGW] = kssolv2GW(optionsGW, mol);
save(GWinputfile, 'GWinput', 'optionsGW');


