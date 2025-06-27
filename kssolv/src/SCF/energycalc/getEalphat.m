function Ealphat = getEalphat(mol)
% GETEALPHAT Calculate the Ealpha energy.
%    Ealphat = getealphat(mol) get the Ealpha energy
%    where mol is a Molecule or Crystal object.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

natoms = mol.natoms;
vol    = mol.vol;
ntypes = length(natoms);

ealpha   = zeros(1,ntypes);

for j = 1:ntypes
    % TODO: the evaluation of the vloc here is not correct.
%    anum = inz(j);
%    pp = PpData(anum);
%    vloc = pp.vloc;
%    r = pp.r;
%    ch = pp.venum;
%    nrr = size(vloc,1);
%    i15 = find(r(2:nrr-1)<15.0)+1;
%    s = sum( (ch*r(i15)+vloc(i15).*r(i15).^2).*(r(i15+1)-r(i15-1))/2 );
%    ealpha(inz(j)) = s*4*pi;
end

Ealphat = sum(ealpha.*natoms);
Ealphat = Ealphat*mol.nel/vol;

end
