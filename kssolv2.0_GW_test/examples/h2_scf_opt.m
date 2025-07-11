%
% This script shows how to use the relaxatom function to optimize the bond
% length of the H2 molecule

clc;
clear all;
close all;
%
% choose the ONCV pseudopotential
%
kssolvpptype('oncv');
%
% construct the H2 molecule 
%
% 1. construct atoms
%
a1 = Atom('H');
atomlist = [a1 a1];
%
% 2. set up supercell
%
C = 10*eye(3);
%
% 3. define the coordinates the atoms 
%
xyzlist = [
 1.5     0.0      0.0
 0.0     0.0      0.0
];
%
% 4. Configure the molecule and set the kinetic energy cutoff
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
               'ecutwfc',50.0, 'name','h2' );
%
% call scf to compute the ground state of the inital atomic configuration
%
[mol,H,X,info] = scf(mol);

% reset some of the scf options for subsequent calls
ksopts = setksopt;
ksopts = setksopt(ksopts,'scftol',1e-7,...
                 'cgtol',1e-8,...
                 'maxscfiter',100, ...
                 'maxcgiter',100, ...
                 'relaxmethod','fminunc', ...
                 'relaxtol', 1e-4);
%
% measure the initial bond length between two H atoms
%
bondlength = norm(mol.xyzlist(1,:)-mol.xyzlist(2,:));
fprintf('starting bond length = %11.3e (Bohr)\n', bondlength);

%
% perform geometry optimization of the bond length between 2 H atoms
%
[molopt,Hopt,Xopt,infoopt] = relaxatoms(mol,ksopts,H,X);
%
% measure the optimized bond length
%
bondlength = norm(molopt.xyzlist(1,:)-molopt.xyzlist(2,:));
fprintf('optimized bond length = %11.3e (Bohr)\n', bondlength);
