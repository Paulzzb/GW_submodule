%
% construct Al4 
%
clear;
kssolvpptype('default');
%
% 1. construct atoms
%
a1 = Atom('Al');
atomlist = [a1, a1, a1, a1];
%
% 2. set up supercell
%
C = 7.632*eye(3);
%
% 3. define the coordinates the atoms 
%
coefs = [
    0.000000000         0.000000000         0.000000000
    0.000000000         0.500000000         0.500000000
    0.500000000         0.000000000         0.500000000
    0.500000000         0.500000000         0.000000000
];
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
degauss = 0.0038;
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',100, 'name','Al4','n1',72,'n2',72,'n3',72,'funct','PBE','temperature',300,'nbnd', 32,'smear','gaussian');
