%
% construct Si8 
%
kssolvpptype('pz-hgh', 'UPF');
%
% 1. construct atoms
%
a1 = Atom('Si');
atomlist = [a1, a1, a1, a1, a1, a1, a1, a1];
%
% 2. set up supercell
%
C = 10.216*eye(3);
%
% 3. define the coordinates the atoms 
%
coefs = [
    0.000000000         0.000000000         0.000000000
    0.000000000         0.500000000         0.500000000
    0.500000000         0.000000000         0.500000000
    0.500000000         0.500000000         0.000000000
    0.750000000         0.250000000         0.750000000
    0.250000000         0.250000000         0.250000000
    0.250000000         0.750000000         0.750000000
    0.750000000         0.750000000         0.250000000
];
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',12.5, 'name','Si8' );
