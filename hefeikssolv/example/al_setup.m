%
% construct an Aluminum cluster
%

%
% 1. construct atoms
%
a = Atom('Al');
atomlist = [a a a a a a a a];
%
% 2. set up supercell
%
C = 20*eye(3);
%
% 3. define the coordinates the atoms
%
xyzlist = [
0 0 0
10 0 0
0 10 0
0 0 10
10 10 10
10 10 0
0 10 10
10 0 10
];
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',12.5, 'name','Aluminum cluster' );
