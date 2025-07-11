%
% This script shows how to use the scf function to compute the ground state 
% of the SiH4 molecule

%
% construct the SiH4 (Silane) molecule 
%
%kssolvpptype('default');
%
% 1. construct atoms
%
a1 = Atom('Si');
a2 = Atom('H');
atomlist = [a1 a2 a2 a2 a2];
%
% 2. set up supercell
%
C = 10*eye(3);
%
% 3. define the coordinates the atoms 
%
redxyz = [
 0.0     0.0      0.0
 0.161   0.161    0.161
-0.161  -0.161    0.161
 0.161  -0.161   -0.161
-0.161   0.161   -0.161
];
xyzlist = redxyz*C';
%
% 4. Configure the molecule and set the kinetic energy cutoff
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
               'ecutwfc',25, 'name','SiH4' );

%
% call the scf function to perform the SCF calcuation
[mol,H,X,info] = scf(mol);
