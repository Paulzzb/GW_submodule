function idx = grid2idx(mol,ga)
%
% idx = grid2idx(mol,ga)
%
% converts each reciprocal (G) space grid point (gx,gy,gz) in ga 
% to an index that gives the position of (gx,gy,gz) in the list
% of reciprocal (G) lattice grid points within the Ecut radius
% contained in the Ggrid object associated with mol.
%
% ga  --- an nga by 3 array of integer reciprocal (G) space grid points 
%         each row of ga is in the form of [igx igy igz]. 
%         The actual grid point is 
%              [igx*2*pi/Lx igy*2*pi/Ly igz*2*pi/Lz]
% idx --- an integer that gives the position of each row in ga
%         in the list of all reciprocal (G) space grid points defined
%         by the Ggrid object associated with mol.
%         

C = get(mol,'supercell');
Lx = norm(C(:,1));
Ly = norm(C(:,2));
Lz = norm(C(:,3));

ggrid = Ggrid(mol);

nga = size(ga,1);
idx = -1*ones(nga,1);

for j = 1:nga
   gx = ga(j,1)*2*pi/Lx;
   gy = ga(j,2)*2*pi/Ly;
   gz = ga(j,3)*2*pi/Lz;

   % should check to see if [gx gy gz] is valid

   hits = (gx == ggrid.gkx).*(gy == ggrid.gky).*(gz == ggrid.gkz);
   idx(j) = find(hits);
end
