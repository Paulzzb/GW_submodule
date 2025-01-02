%
% construct the Benzene molecule (ethane)
%
%kssolvpptype('pz-hgh', 'UPF');
%
% 1. construct atoms
%
a1 = Atom('H');
a2 = Atom('C');
atomlist = [a1 a1 a1 a1 a1 a1 a2 a2 a2 a2 a2 a2];
%
% 2. set up supercell
%
C = [22.676722 0  0
     0 22.676722  0
     0  0 22.676722];
%
% 3. define the coordinates the atoms 
%
coefs = [
     0.292420000         0.500000000         0.500000000
     0.396239996         0.320270002         0.500000000
     0.603760004         0.320270002         0.500000000
     0.707579970         0.500000000         0.500000000
     0.603760004         0.679729998         0.500000000
     0.396239996         0.679729998         0.500000000
     0.383379996         0.500000000         0.500000000
     0.441709995         0.399040014         0.500000000
     0.558290005         0.399040014         0.500000000
     0.616620004         0.500000000         0.500000000
     0.558290005         0.600960016         0.500000000
     0.441709995         0.600960016         0.500000000
];
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',20.0, 'name','Benzene' );
