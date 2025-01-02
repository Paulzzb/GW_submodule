%
% construct the NiC molecule
%

%
% 1. construct atoms
%
a1 = Atom('Ni');
a2 = Atom('C');
atomlist = [a1 a2];
%
% 2. set up the supercell
%
C = 5*eye(3);
%
% 3. define the coordinates the atoms in a.u.
%
xyzlist = [
0.0 0.0 0.0
1.8 0.0 0.0
]/0.529177;
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',12.5, 'name','NiC' );