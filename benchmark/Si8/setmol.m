%
% construct Si2 
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
C = [10.3344  0.0000  0.0000
      0.0000 10.3344  0.0000
      0.0000  0.0000 10.3344];
%
% 3. define the coordinates the atoms 
%
coefs = [
   0.250000 0.750000 0.250000
   0.000000 0.000000 0.500000
   0.250000 0.250000 0.750000
   0.000000 0.500000 0.000000
   0.750000 0.750000 0.750000
   0.500000 0.000000 0.000000
   0.750000 0.250000 0.250000
   0.500000 0.500000 0.500000
    ];
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut', 10, 'name','Si8', 'nbnd', 36);
