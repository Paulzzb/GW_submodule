function VexxACE = calculateACE(Vexx, X)
%
% Usage: VexxACE = calculateACE(Vexx, X)
%
% Purpose:
%    Computes the Adaptively Compressed Exchange Operator
%
% Input:
%    Vexx --- Exact Exchange Operator
%    X --- Wavefunction
%
% Ouptut:
%    VexxACE --- ACE Exchange Operator
%
    
    W = Vexx(X);
    M = X' * W;
    M = (M + M')/2;
    R = chol(-M);
    Xi = W / R;
    VexxACE = @(x) -Xi * (Xi' * x);
    
end