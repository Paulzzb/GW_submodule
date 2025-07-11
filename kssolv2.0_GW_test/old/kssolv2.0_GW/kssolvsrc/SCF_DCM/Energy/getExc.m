function Exc = getExc(mol,rho,uxc2,uxcsr)
%
% Usage: Exc = getExc(mol,rho,vexc);
%
% Purpose:
%    Compute the exchange correlation energy
%
% Input:
%  mol  --- Molecule object
%  rho  --- input charge density (3D array)
%  vexc --- exchange correlation potential (3D array)
%
% Output:
%  Exc  --- Exchange correlation energy (scalar)
%
hhh = mol.vol/mol.n1/mol.n2/mol.n3;
if nargin < 4
    Exc = sumel(uxc2.*rho)*hhh;
else
    Exc = sumel(uxc2.*rho-0.25*uxcsr)*hhh;
end

end
