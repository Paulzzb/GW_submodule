%
% construct the SiH4 (Silane) molecule 
%
kssolvpptype('default');
%
% 1. construct atoms
%
a1 = Atom('Si');
a2 = Atom('H');
atomlist = [a1 a2 a2 a2 a2];
%
% 2. set up supercell
%
C = 10*eye(3);
%
% 3. define the coordinates the atoms 
%
coefs = [
 0.0     0.0      0.0
 0.161   0.161    0.161
-0.161  -0.161    0.161
 0.161  -0.161   -0.161
-0.161   0.161   -0.161
];
xyzlist = coefs*C';
autokpts = [ 2 2 2 0 0 0 ];
%
% 4. Configure the molecule (crystal)
%
cry = Crystal('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',12.5, 'name','SiH4', 'autokpts',autokpts );

[cry,H,X,info] = scf(cry);
