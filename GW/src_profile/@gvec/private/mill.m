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
qpoint = mill_in.qpoint;


[gkxind, gkyind, gkzind] = ...
  ndgrid((0:n1-1)-((0:n1-1) >= n1/2)*n1, ...
    (0:n2-1)-((0:n2-1) >= n2/2)*n2, ...
    (0:n3-1)-((0:n3-1) >= n3/2)*n3);
gkxind = gkxind(:);
gkyind = gkyind(:);
gkzind = gkzind(:);
idlist = 1:length(gkxind);


gkind = [gkxind, gkyind, gkzind];
rcplat = 2*pi*inv(supercell);
gkvec = gkind*rcplat';
qgkvec = gkvec + qpoint;
qgkabs = sum(qgkvec.^2, 2);
idxnz = find(qgkabs <= ecut);

% Sort according to |G|.^2
gkxind = gkxind(idxnz);
gkyind = gkyind(idxnz);
gkzind = gkzind(idxnz);
idlist = idlist(idxnz);

qgkabs = qgkabs(idxnz);
[~, sort_ind] = sort(qgkabs);





gvec_in.idxnz = idxnz(sort_ind);
gvec_in.ng = length(idxnz);
gvec_in.components = [gkxind(sort_ind), gkyind(sort_ind), gkzind(sort_ind)];


end % EOF
