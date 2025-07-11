function exxgkk = getExxgkk(mol,ecut)
%
% Usage:
%    exxgkk = getExxgkk(mol);
%
% Purpose:
%    Compute eigenvalues of the exact exchange operator
%
% Input:
%    mol  --- Molecule information
%    ecut --- ecut for reciprocal grids
%
% Output:
%    Exxgkk --- Eigenvalues of the exact change operator for solving Poisson
%               equations
%
n1 = mol.n1;
n2 = mol.n2;
n3 = mol.n3;

if (nargin < 2) 
  ecut = mol.ecutwfc;
end

grid = Ggrid(mol,ecut);

% Hard-coded constants
Mu = 0.106;     % Screening Parameter
epsDiv = 1e-8;  % Divergence threshold

% Gygi-Baldereschi regularization
exxAlpha = 10 / (grid.ecut * 2);
gkki = exp(-exxAlpha * grid.gkk) ./ grid.gkk .* (1 - exp(-grid.gkk / (4*Mu^2)));
idx = grid.gkk < epsDiv;
gkki(idx) = 0;
exxDiv = 4 * pi * (sum(gkki) + 1/(4*Mu^2));
nqq  = 1e5; % What is this?
dq = 5 / sqrt(exxAlpha) / nqq;
qt2 = (((1:nqq+1)' - 0.5) * dq).^2;
aa = -sum(exp(-(exxAlpha + 1/(4*Mu^2)) * qt2)) * dq * 2 / pi + 1/sqrt(exxAlpha * pi);
exxDiv = exxDiv - mol.vol * aa;

exxgkk = 4 * pi ./ grid.gkk .* (1 - exp(-grid.gkk / (4*Mu^2)));
exxgkk(idx) = -exxDiv + pi / Mu^2;

end





