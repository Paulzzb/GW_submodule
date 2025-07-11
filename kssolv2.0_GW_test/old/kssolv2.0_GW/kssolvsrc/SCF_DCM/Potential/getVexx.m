function Vexxf = getVexx(X, mol, dfrank, exxgkk, F)
%
% Usage: Vexx = getVexx(X, F, exxgkk, mol)
%
% Purpose:
%    Computes exact exchange operator
%
% Input:
%    X  --- Wavefunction
%    F  --- Fourier Transform
%    dfrank --- Rank for density fitting
%    exxgkk --- Eigenvalue of exchagne operator
%    mol --- Molecule information
%
% Ouptut:
%    Vexx --- Exact exchange operator
%
if nargin < 5, F = KSFFT(mol);  end
if nargin < 4, exxgkk = getExxgkk(mol); end

% n123 = mol.n1 * mol.n2 * mol.n3;
nocc = mol.nel/2*mol.nspin;
Phi = F' * X.psi(:,1:nocc);
if dfrank
    Vexxf = @(Psi)Vexxdf(Psi, Phi, dfrank, F, exxgkk, mol);
else
    Vexxf = @(Psi)Vexx(Psi, Phi, F, exxgkk, mol);
end

end
