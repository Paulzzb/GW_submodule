clc;
clear all;
close all;
%
% 0. Choose pseudopotential
%
kssolvpptype('default');
%
% 1. Construct atoms
%
a1 = Atom('C');
a2 = Atom('H');
atomlist = [a1, a1, a2, a2, a2, a2];
%
% 2. Set up supercell
%
cell_A = 10; Bohr = 0.529177210905667; cell_au = cell_A/Bohr;
C = cell_au*eye(3);
%
% 3. Define the coordinates of atoms 
%
coefs = [
     0.429980000         0.526960000         0.478800000
     0.583350000         0.531610000         0.491880000
     0.381600000         0.525080000         0.375590000
     0.364820000         0.525400000         0.572330000
     0.648520000         0.533200000         0.398350000
     0.631720000         0.533460000         0.595090000
];
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist',xyzlist, ...
    'ecut',20, 'name','C2H4', 'funct','PBE', 'n1',80, 'n2',80, 'n3',80, ...
    'extranbnd',0, 'nspin',1);
%
% 5. Configure the options
%
opt = setksopt();

opt.eigmethod = 'davidson_qe';
opt.scftol    = 1e-9;
opt.phitol    = 1e-9;
opt.maxscfiter= 300;
opt.maxphiiter= 300;
opt.mixtype   = 'broyden_modify';
opt.what2mix  = 'rho';
opt.betamix   = 0.7;
opt.mixdim    = 8;
opt.brank     = 8;

opt.exxmethod = 'ace';
opt.x_gamma_extrapolation = false;

opt.relaxmethod = 5; % 1√2△3√4√5√
opt.relaxtol = 1e-4; % Tolerance of 2-Norm Force. In Gaussian, 3e-4*sqrt(3*natoms) (default) and 1e-5*sqrt(3*natoms) (tight) are used.
opt.maxrelaxiter = 200;

mass = 4;                           % In QE, mass = 1822.88848684726
dtStart = 0.05 * 41.3413733349185; 	% 1 Femtoseconds [fs] = 41.3413733349185 atomic time unit
FIREpars.nMin       = 5;
FIREpars.fDec       = 0.5;
FIREpars.fInc       = 1.1;
FIREpars.alphaStart = 0.1;
FIREpars.fAlpha     = 0.99;
FIREpars.dtMax      = 10;           % Note: dtMax = dtStart * FIREpars.dtMax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[molopt,Hopt,Xopt,infoopt] = relaxatoms(mol,opt,[],[],mass,dtStart,FIREpars);
relaxHistory = infoopt.relaxHistory;
[molopt_cont,Hopt_cont,Xopt_cont,infoopt_cont] = relaxatoms(molopt,opt,Hopt,Xopt,mass,dtStart,FIREpars);
relaxHistory_cont = infoopt_cont.relaxHistory;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%