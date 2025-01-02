%
% construct the CO molecule
%
%
% 1. construct atoms
%
a1 = Atom('C');
a2 = Atom('O');
atomlist = [a1 a2];
%
% 2. set up the supercell
%
BoxSize = 13.228;
C = BoxSize*eye(3);
%
% 3. define the coordinates the atoms 
%
xyzlist = [
 -1.066 0.0 0.0
  1.066 0.0 0.0
];
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',30, 'name','CO' );
