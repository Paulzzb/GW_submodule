%
% construct the CO2 molecule
%

%
% 1. construct atoms
%
a1 = Atom('C');
a2 = Atom('O');
atomlist = [a1 a2 a2];
%
% 2. set up the supercell
%
BoxSize = 10;
C = BoxSize*eye(3);
%
% 3. define the coordinates the atoms 
%
coefs = [
0.000000000000    -0.026381925226     0.000000000000
2.078381991574     0.009896367383     0.000000000000
-2.078381991574     0.009896367383     0.000000000000
];
coefs = coefs/BoxSize;
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',12.5, 'name','CO2' );