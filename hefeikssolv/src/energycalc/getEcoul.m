function Ecoul = getEcoul(mol,rho,vhart)
%
% Usage: Ecoul = getEcoul(mol, rho, vhart);
%
% Purpose:
%    Compute the Coulomb potential energy induced by vhart
%
% Input:
%  mol    --- a Molecule object
%  rho    --- input charge density (3D array or cell array)
%  vhart  --- Hartee (Coulomb) potential (3D array)
%
% Output:
%  Ecoul  --- Coulomb potential energy induced by vhart.
%
nspin = mol.nspin;
hhh = mol.vol/mol.n1/mol.n2/mol.n3;

if nspin == 1
    arho = abs(rho);
elseif nspin == 2
    arho = abs(rho{1}+rho{2});
elseif nspin == 4
    arho = abs(rho{1});
end

Ecoul      = sumel(arho.*vhart)/2;
Ecoul      = real(Ecoul)*hhh;
