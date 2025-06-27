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
a1 = Atom('H');
atomlist = [a1, a1];
%
% 2. Set up supercell
%
cell_au = 10; % cell_A is too small, only for test
C = cell_au*eye(3);
%
% 3. Define the coordinates of atoms 
%
bond = 1.5/cell_au;
coefs = [
   0.5000000000000000    0.5000000000000000    0.5-bond/2
   0.5000000000000000    0.5000000000000000    0.5+bond/2
];
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist',xyzlist, ...
    'ecut',20, 'name','H2', 'funct','PZ', 'n1',45, 'n2',45, 'n3',45, ...
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
dtStart = 0.05 * 41.3413733349185;  % 1 Femtoseconds [fs] = 41.3413733349185 atomic time unit
FIREpars.nMin       = 5;
FIREpars.fDec       = 0.5;
FIREpars.fInc       = 1.1;
FIREpars.alphaStart = 0.1;
FIREpars.fAlpha     = 0.99;
FIREpars.dtMax      = 10;           % Note: dtMax = dtStart * FIREpars.dtMax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% measure the initial bond length between two H atoms
bondlength = norm(mol.xyzlist(1,:)-mol.xyzlist(2,:));
fprintf('starting bond length = %11.8e (Bohr)\n', bondlength);

% perform geometry optimization of the bond length between two H atoms
[molopt,Hopt,Xopt,infoopt] = relaxatoms(mol,opt,[],[],mass,dtStart,FIREpars);
relaxHistory = infoopt.relaxHistory; % The information of each optimization step can be found in infoopt.relaxHistory

% measure the optimized bond length between two H atoms
bondlength = norm(molopt.xyzlist(1,:)-molopt.xyzlist(2,:));
fprintf('optimized bond length = %11.8e (Bohr)\n', bondlength);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%