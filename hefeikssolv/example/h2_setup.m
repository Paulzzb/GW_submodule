%
% construct the H2 molecule 
%
kssolvpptype('default');
%
% 1. construct atoms
%
a1 = Atom('H');
atomlist = [a1 a1];
%
% 2. set up supercell
%
C = 10*eye(3);
%
% 3. define the coordinates the atoms 
%
xyzlist = [
  0.75     0.0      0.0
 -0.75     0.0      0.0
];

%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',5.0, 'name','h2' );


