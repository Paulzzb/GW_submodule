function ind = end(X,k,n)
% WAVEFUN/END End function for wave function class
%    Y = X(2:end,:) returns the sub wave function.
%
%    See also Wavefun.

if X.trans
    szd = size(X.psi');
else
    szd = size(X.psi);
end

if k < n
    ind = szd(k);
else
    ind = prod(szd(k:end));
end

end