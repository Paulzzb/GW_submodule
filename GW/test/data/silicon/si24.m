%
% construct Si2 
%
kssolvpptype('pz-hgh', 'UPF');
%
% 1. construct atoms
%
a1 = Atom('Si');
atomlist = [
a1, a1, a1, a1, a1, a1, a1, a1, ...
a1, a1, a1, a1, a1, a1, a1, a1, ...
a1, a1, a1, a1, a1, a1, a1, a1, ...
];
%
% 2. set up supercell
%
C = [
     31.0032  0.0000  0.0000
      0.0000 10.3344  0.0000
      0.0000  0.0000 10.3344
];
%
% 3. define the coordinates the atoms 
%
coefs = [
    0.0833333333333334    0.7500000000000000    0.2500000000000000
    0.4166666666666667    0.7500000000000000    0.2500000000000000
    0.7500000000000000    0.7500000000000000    0.2500000000000000
    0.0000000000000000    0.0000000000000000    0.5000000000000000
    0.3333333333333334    0.0000000000000000    0.5000000000000000
    0.6666666666666667    0.0000000000000000    0.5000000000000000
    0.0833333333333334    0.2500000000000000    0.7500000000000000
    0.4166666666666667    0.2500000000000000    0.7500000000000000
    0.7500000000000000    0.2500000000000000    0.7500000000000000
    0.0000000000000000    0.5000000000000000    0.0000000000000000
    0.3333333333333334    0.5000000000000000    0.0000000000000000
    0.6666666666666667    0.5000000000000000    0.0000000000000000
    0.2500000000000000    0.7500000000000000    0.7500000000000000
    0.5833333333333334    0.7500000000000000    0.7500000000000000
    0.9166666666666666    0.7500000000000000    0.7500000000000000
    0.1666666666666667    0.0000000000000000    0.0000000000000000
    0.5000000000000000    0.0000000000000000    0.0000000000000000
    0.8333333333333334    0.0000000000000000    0.0000000000000000
    0.2500000000000000    0.2500000000000000    0.2500000000000000
    0.5833333333333334    0.2500000000000000    0.2500000000000000
    0.9166666666666666    0.2500000000000000    0.2500000000000000
    0.1666666666666667    0.5000000000000000    0.5000000000000000
    0.5000000000000000    0.5000000000000000    0.5000000000000000
    0.8333333333333334    0.5000000000000000    0.5000000000000000
    ];
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',10, 'name','Si24', 'n1', 90, 'n2', 30, 'n3', 30, 'nbnd', 96);