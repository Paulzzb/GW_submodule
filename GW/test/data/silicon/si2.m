%
% construct Si2 
%
kssolvpptype('pz-hgh', 'UPF');
%
% 1. construct atoms
%
a1 = Atom('Si');
atomlist = [a1, a1];
%
% 2. set up supercell
%
C = 5.108*eye(3);
%
% 3. define the coordinates the atoms 
%
coefs = [
    0.000000000         0.000000000         0.000000000
%    0.500000000         0.500000000         0.500000000
    0.600000000         0.500000000         0.500000000
    ];
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut', 5, 'name', 'Si2', 'n1', 12, 'n2', 12, 'n3', 12, 'nbnd', 10);
mol_rho = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut', 20, 'name', 'Si2', 'n1', 12, 'n2', 12, 'n3', 12, 'nbnd', 10);
