% This script simply shows how to set up a SiH4 molecule 
% that can be passed to scf or trdcm for ground state calculations
%

%
% construct the SiH4 (Silane) molecule 
%
%kssolvpptype('pz-hgh', 'UPF');
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
