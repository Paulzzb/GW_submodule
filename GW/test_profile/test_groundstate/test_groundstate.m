clear;
%
INFO_DIR = './';
% Please change the kspath to yours

cfile = mfilename('fullpath');
CPATH = fileparts(cfile);
CPATH = [CPATH, '/'];
% KS_DIR    = ['D:/kssolvGW/GW_submodule/hefeikssolv/'];
KS_DIR = '..\..\..\kssolv\';
QP_DIR = '..\..\';
molinfile = [CPATH, INFO_DIR, 'setmol.m'];
scfinfile = [CPATH, INFO_DIR, 'setscf.m'];
addpath('../') % Add savefunc into search path

%
% 1. Managing paths.

cd(KS_DIR);
KSSOLV_startup;
cd(CPATH);





run(molinfile);
run(scfinfile);
[mol,H,X0,info] = scf(mol, options_scf);
save_groundstate_to_GWformat(mol, H, X0, info, './');


