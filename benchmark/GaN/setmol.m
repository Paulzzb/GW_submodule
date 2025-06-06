%
% construct Si2 
%
kssolvpptype('pz-hgh', 'UPF');
%
% 1. construct atoms
%
a1 = Atom('Ga');
a2 = Atom('N');
atomlist = [a1, a2];
%
% 2. set up supercell
%
C = [
  5.5065   0.0000   0.0000       % a ≈ 3.189 Å = 6.028 Bohr
 -2.7533   4.7692   0.0000       % a * cos(120°), a * sin(120°)
  0.0000   0.0000   9.8106       % c ≈ 5.185 Å = 9.790 Bohr
];
%
% 3. define the coordinates the atoms 
%
coefs = [
   0.000000 0.000000 0.500000
   0.333333  0.666667  0.375000
    ];
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut', 50, 'name','GaN', 'nbnd', 40);
