% This script shows how to use trust-region enabled direct constrained
% minimization (TRDCM) to compute the ground state of the SiH4 molecule
%
% select default pseudopotentials
%
kssolvpptype('default');

%
% construct the SiH4 (Silane) molecule 
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
coefs = [
 0.0     0.0      0.0
 0.161   0.161    0.161
-0.161  -0.161    0.161
 0.161  -0.161   -0.161
-0.161   0.161   -0.161
];
xyzlist = coefs*C';
%
% 4. Configure the molecule and set the kinetic energy cutoff
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecutwfc',12.5, 'name','SiH4' );
%
% reset some options for trdcm and call trdcm
%
opts = setksopt();
opts.verbose = 'on';
opts.dcmtol = 1e-8;
opts.maxdcmiter = 30;
[mol,H,X,info] = trdcm(mol,opts);
