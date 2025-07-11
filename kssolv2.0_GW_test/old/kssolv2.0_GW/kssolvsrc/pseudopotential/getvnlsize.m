function [siz,offset] = getvnlsize(mol)
% GETVNLSIZE returns the size of non-local pseudo potential.
%    [vnlmat,vnlsign] = GETVNLSIZE(mol) returns the size of the non-local
%    pseudo potential for each atom of the molecule.
%
%    See also getvloc, getvnl.


ppvar = mol.ppvar;
alist = mol.alist;
ntypes = length(mol.atoms);
na = length(alist);

siz = zeros(na,1);
for itype = 1:ntypes
    totall = 0;
    for it = 1:length(ppvar.lll{itype})
        ll = ppvar.lll{itype}(it);
        totall = totall + 2*ll+1;
    end
    index  = find(alist == itype);
    for ita = index'
        siz(ita) = totall;
    end
end

offset = cumsum(siz);
offset = [0; offset(1:end-1)];

end
