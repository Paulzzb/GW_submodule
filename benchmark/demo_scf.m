clear;
%
INFO_DIR = 'Si8/';
% Please change the kspath to yours
cfile = mfilename('fullpath');
CPATH = fileparts(cfile);
CPATH = [CPATH, '/'];
KS_DIR    = ['D:/kssolvGW/GW_submodule/hefeikssolv/'];
molinfile = [CPATH, INFO_DIR, 'setmol.m'];
scfinfile = [CPATH, INFO_DIR, 'setscf.m'];
INTER_DIR = [CPATH, INFO_DIR, 'IntermediateFiles/'];
%
scfoutfile = [INTER_DIR, 'scinfo.mat'];
moloutfile = [INTER_DIR, 'molinfo.mat'];
%
% 1. Managing paths.

cd(KS_DIR);
KSSOLV_startup;
cd(CPATH);


run(molinfile);
run(scfinfile);
[mol,H,X0,info] = scf(mol, options_scf);
save(scfoutfile, 'mol', 'H', 'X0', 'info', '-v7.3');
save(moloutfile, 'mol');


