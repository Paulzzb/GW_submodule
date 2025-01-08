clear all;
close all;
%
INFO_DIR = 'data/Si8/';

CPATH = mfilename('fullpath');
CPATH = fileparts(CPATH);
[GWPATH, ~, ~] = fileparts(CPATH);
CPATH = [CPATH, '/'];
INTER_DIR = [CPATH, INFO_DIR, 'IntermediateFiles/'];
GWinputfile = [INTER_DIR, 'GWinput_fullfreq_cd_isdf.mat'];
 
cd(GWPATH); gw_startup(); cd(CPATH);

% 2. Load
load(GWinputfile);

% 3. run
gwCalculation(GWinput, optionsGW);
