function [occ,efermi] = getocc(ev,nocc,Tbeta,smear)
% GETOCC calculates occupation ratio for each eigenvalues.
%    [occ,efermi] = GETOCC(ev,nocc,Tbeta,smear) returns the occupation ration for
%    each eigenvalues and the fermi parameter efermi such that the sum of
%    smearing function with efermi and ev equals nocc.
%   
%    See also fermidirac.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

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
efermi = (lev+uev)/2;

for iter = 1:maxiter
    occ = wgauss(ev,efermi,Tbeta,smear);
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
