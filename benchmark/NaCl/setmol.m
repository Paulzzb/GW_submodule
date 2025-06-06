%
% construct Si2 
%
% kssolvpptype('pz-hgh', 'UPF');
%
% 1. construct atoms
%
a1 = Atom('Na');
a2 = Atom('Cl');
atomlist = [a1, a2];
%
% 2. set up supercell
%
C = [10.6380  0.0000  0.0000
      0.0000 10.6380  0.0000
      0.0000  0.0000 10.6380];
%
% 3. define the coordinates the atoms 
%
coefs = [
   0.000000 0.000000 0.500000
   0.500000 0.500000 0.500000
    ];
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut', 50, 'name','NaCl', 'nbnd', 40);
