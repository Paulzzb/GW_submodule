function Fscc = getFscc(mol,dv)
% GETFSCC Calculate the force correction term.
%    Fscc = getFscc(mol,dv) calculate the force correction term 
%    owing to the difference of density in last two scf iterations
%    dv = v(out) - v(in).
%
%    This term is very small when SCF convergence is reached.
%
%    See also getFewald, getFloc, getFnl.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

alist = mol.alist;
atoms = mol.atoms;
ntypes = length(atoms);
xyzlist = mol.xyzlist/mol.alat';
na = sum(mol.natoms);
ecut2    = mol.ecut2;

grid2 = Ggrid(mol,ecut2);
gkk   = grid2.gkk;
gkx   = grid2.gkx;
gky   = grid2.gky;
gkz   = grid2.gkz;

if abs(gkk(1)) <= 1e-8
    idxg = 2:length(gkk);
else
    idxg = 1:length(gkk);
end

if mol.nspin==1||(mol.nspin==4 && ~mol.domag)
    psic = dv;
elseif mol.nspin == 2
    psic = (dv{1} + dv{2})/2;
elseif mol.nspin == 4 && mol.domag
    psic = dv{1};
end

psicg = fft3(psic)/mol.n1/mol.n2/mol.n3;
psicg = psicg(grid2.idxnz);
Fscc = zeros(na,3);
for itype = 1:ntypes
    rhog = mol.ppvar.rhog(:,itype);
    for ia = 1:na
        if alist(ia) == itype
            g = [gkx(idxg) gky(idxg) gkz(idxg)];
            arg = mol.alat*g*xyzlist(ia,:)';
            Fscc(ia,:) = real(sum(repmat(rhog(idxg).*complex(sin(arg),cos(arg)).*conj(psicg(idxg)),1,3).*g,1));
        end
    end
end
end
    
    
