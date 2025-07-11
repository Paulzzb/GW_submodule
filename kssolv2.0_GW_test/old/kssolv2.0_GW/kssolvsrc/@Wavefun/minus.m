function Z = minus(X,Y)
% WAVEFUN/MINUS Minus function for wave function class
%    Z = MINUS(X,Y) returns the wave function of X - Y.
%
%    See also Wavefun.

if isa(X,'Wavefun') && isa(Y,'Wavefun')
    Z = X;
    Z.psi = X.psi - Y.psi;
    return;
end

if isa(X,'Wavefun')
    Z = X;
    Z.psi = X.psi - Y;
    return;
end

if isa(Y,'Wavefun')
    Z = Y;
    Z.psi = X - Y.psi;
    return;
end

end
