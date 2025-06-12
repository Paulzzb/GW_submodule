function gvec_in = mill(gvec_in, mill_in)
% This is a subfunction of gvec
% This function use mill_in (as same as gvec_in) to add the
% 'components' and 'idxnz' properties of gvec_in.
% Input: 
%       mill_in: as 'gvec_in' in gvec.m, check gvec.m for details

n1 = mill_in.n1;
n2 = mill_in.n2;
n3 = mill_in.n3;
supercell = mill_in.supercell;
ecut = mill_in.ecut;


[gkxind, gkyind, gkzind] = ...
  ndgrid((0:n1-1)-((0:n1-1) >= n1/2)*n1, ...
    (0:n2-1)-((0:n2-1) >= n2/2)*n2, ...
    (0:n3-1)-((0:n3-1) >= n3/2)*n3);
gkxind = gkxind(:);
gkyind = gkyind(:);
gkzind = gkzind(:);


gkind = [gkxind, gkyind, gkzind];
rcplat = 2*pi*inv(supercell);
gkvec = gkind*rcplat';
gkabs = sum(gkvec.^2, 2);
idxnz = find(gkabs <= ecut);

gvec_in.idxnz = idxnz;
gvec_in.ng = length(gvec_in.idxnz);
gvec_in.components = [gkxind(idxnz), gkyind(idxnz), gkzind(idxnz)];


end % EOF
