%
% construct the H2O molecule
%

%
% 1. construct atoms
%
a1 = Atom('O');
a2 = Atom('H');
atomlist = [a1 a2 a2];
%
% 2. set up the supercell
%
BoxSize = 10;
C = BoxSize*eye(3);
%
% 3. define the coordinates the atoms in a.u.
%
coefs = [
0.000000000000    -0.123909374404     0.000000000000
1.429936611037     0.983265845431     0.000000000000
-1.429936611037     0.983265845431     0.000000000000
];
coefs = coefs/BoxSize; % Scale the atom coordinate relative to BoxSize.
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',10.0, 'name','H2O' );
