%
% construct the HNCO (isocyanic acid) molecule
%

%
% 1. construct atoms
%
a1 = Atom('H');
a2 = Atom('N');
a3 = Atom('C');
a4 = Atom('O');
atomlist = [a1 a2 a3 a4];
%
% 2. set up the supercell
%
BoxSize = 10;
C = BoxSize*eye(3);
%
% 3. define the coordinates the atoms 
%
coefs = [
 3.284938618980    -1.461359490637     0.001371230530
 2.258342591219     0.125197975455     0.001932017162
-0.006404171175     0.010271400781    -0.004893948816
-2.179288390874    -0.025234185924     0.001893804798
];
coefs = coefs/BoxSize;
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',12.5, 'name','HNCO' );
