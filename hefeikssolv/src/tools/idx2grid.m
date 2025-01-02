function ga = idx2grid(mol,idx)
%
% ga = idx2grid(mol,idx)
%
% returns integer reciprocal (G) space lattice points (igx,igy,igz) 
% contained in the Ggrid object associated with mol indexed by idx.
%
% ga  --- an nga by 3 array of integer reciprocal (G) space grid points. 
%         Each row of ga is in the form of [igx igy igz]. 
%         The actual grid point is 
%              [igx*2*pi/Lx igy*2*pi/Ly igz*2*pi/Lz]
% idx --- an integer that gives the position of each row in ga
%         in the list of all reciprocal (G) space grid points defined
%         by the Ggrid object associated with mol.
%         
%
C = get(mol,'supercell');
Lx = norm(C(:,1));
Ly = norm(C(:,2));
Lz = norm(C(:,3));

ggrid = Ggrid(mol);

igkx = ggrid.gkx*Lx/(2*pi);
igky = ggrid.gky*Ly/(2*pi);
igkz = ggrid.gkz*Lz/(2*pi);
igkk = igkx.^2 + igky.^2 + igkz.^2;

if (nargin < 2)
   ga = [igkx igky igkz igkk]
else
   ga = [igkx(idx) igky(idx) igkz(idx) igkk(idx)];
end
