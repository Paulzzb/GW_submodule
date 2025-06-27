function Floc = getFloc(mol,rho)
% GETFLOC Calculate the local potential contribution to force.
%    Floc = getFloc(mol,rho) calculate the local potential contribution 
%    to force which is determined only by total electron density.
%
%    See also getFewald, getFnl, getFscc.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

alist = mol.alist;
xyzlist = mol.xyzlist;
na = sum(mol.natoms);
ecut2    = mol.ecut2;
vol = mol.vol;

grid2 = Ggrid(mol,ecut2);
gkx   = grid2.gkx;
gky   = grid2.gky;
gkz   = grid2.gkz;
idxnz = grid2.idxnz;

rhog = fft3(rho)/(mol.n1*mol.n2*mol.n3);
rhog = rhog(idxnz);

Floc = zeros(na,3);
for it = 1:na
    ia = alist(it);
    arg = [gkx gky gkz]*xyzlist(it,:)';
    vlocg = mol.ppvar.vlocg(:,ia)*4/vol*pi;
    argrhog = sin(arg).*real(rhog)+cos(arg).*imag(rhog);
    Floc(it,1) = sum(gkx.*vlocg.*argrhog);
    Floc(it,2) = sum(gky.*vlocg.*argrhog);
    Floc(it,3) = sum(gkz.*vlocg.*argrhog);
end

Floc = Floc*vol;

end
