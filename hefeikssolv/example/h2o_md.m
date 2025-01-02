clc;
clear all;
close all;
%
% 1. construct atoms
%
a1 = Atom('O');
a2 = Atom('H');

atomlist = [a1 a2 a2];
%
% 2. set up the supercell
%
BoxSize = 20;
C = BoxSize*eye(3);
%
% 3. define the coordinates the atoms in a.u.
%
coefs = [
+1.03346393e+01 +8.64870918e+00 +1.08235485e+01
+9.94418316e+00 +1.15205130e+01 +1.05705583e+01
+9.90257706e+00 +9.95417846e+00 +9.56849266e+00
];
coefs = coefs/BoxSize; % Scale the atom coordinate relative to BoxSize.
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',60, 'name','H2O' );

opt = setksopt();
opt.maxscfiter = 200;
opt.betamix = 0.08;
opt.mixdim = 6;
opt.scftol = 1e-6;

opt = setksopt();
opt = setksopt(opt,'mixtype','broyden','eigmethod','eigs');
[mol,H0,X0,info] = scf(mol, opt);

nMdSteps = 2;

[X, V, F] = md(mol, H0, X0, nMdSteps)
