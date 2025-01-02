%
% construct the LiH dimer 
%
kssolvpptype('ONCV_PBE')
%
% 1. construct atoms
%
a1 = Atom('H');
a2 = Atom('Li');
atomlist = [a1 a2];
%
% 2. set up the supercell
%
BoxSize = 10;
C = BoxSize*eye(3);
%
% 3. define the coordinates the atoms in a.u.
%
xyzlist = [
-1.682 0.000 0.000
 1.682 0.000 0.000
];
%%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',30, 'name','LiH' );

[molupf,~,~,~]=scf(mol);



%%
% 5. Configure the molecule (crystal)
%
kssolvpptype('LDA','psp8')

mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',30, 'name','LiH' );

[molpsp,~,~,~]=scf(mol);