%
% construct the LiH dimer 
%
kssolvpptype('LDA','psp8')
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
bondVec = [(0:10)*0.05+1.482];
energyVec = zeros(size(bondVec));

for iList =  1 : length(bondVec)
  bond = bondVec(iList);
  xyzlist = [
  -bond 0.000 0.000
   bond 0.000 0.000
  ];
  mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',80, 'name','LiH' );
  [molpsp,~,~,info]=scf(mol);
  energyVec(iList) = info.Etotvec(end);
end

plot(bondVec, energyVec-min(energyVec),'r-p');
