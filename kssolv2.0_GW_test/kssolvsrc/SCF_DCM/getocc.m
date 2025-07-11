function [occ,efermi] = getocc(ev,nocc,Tbeta)
% GETOCC calculates occupation ratio for each eigenvalues.
%    [occ,efermi] = GETOCC(ev,nocc,Tbeta) returns the occupation ration for
%    each eigenvalues and the fermi parameter efermi such that the sum of
%    Fermi-Dirac function with efermi and ev equals nocc.
%   
%    See also fermidirac.

tol = 1e-15;
maxiter = 200;

if Tbeta < eps
    [ev,idx] = sort(ev);
    occ = zeros(size(ev));
    occ(1:nocc) = 1;
    occ(idx) = occ;
    efermi = ev(nocc);
    return;
end

lev = min(ev);
uev = max(ev);
efermi = (lev*uev)/2;

for iter = 1:maxiter
    occ = fermidirac(ev,efermi,Tbeta);
    occsum = sum(occ);
    if abs(occsum - nocc) < tol*nocc
        return;
    elseif occsum < nocc
        lev = efermi;
    else
        uev = efermi;
    end
    efermi = (lev+uev)/2;
end

end
