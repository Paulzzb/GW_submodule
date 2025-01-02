%
% construct an electron quantum dot confined by 
% a harmonic external potential 
%

%
% 2. set up the primitive cell 
%
C = 10*eye(3);
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, ...
    'ecut',12.5, 'name','Quantum dot' );