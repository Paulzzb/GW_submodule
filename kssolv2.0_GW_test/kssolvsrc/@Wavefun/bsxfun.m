function Z = bsxfun(fun,X,Y)
% WAVEFUN/BSXFUN Bsxfun function for wave function class
%    Z = BSXFUN(fun,X,Y) returns a wave function as the fun operation of
%    two wave functions.
%
%    See also Wavefun, bsxfun.

if isa(X,'Wavefun') && isa(Y,'Wavefun')
    Z = X;
    Z.psi = bsxfun(fun, X.psi, Y.psi);
    return;
end

if isa(X,'Wavefun')
    Z = X;
    Z.psi = bsxfun(fun, X.psi, Y);
    return;
end

if isa(Y,'Wavefun')
    Z = Y;
    Z.psi = bsxfun(fun, X, Y.psi);
    return;
end

end
