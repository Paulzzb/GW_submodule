function Exc = getExc(mol,rho,uxc2,uxcsr)
%
% Usage: Exc = getExc(mol,rho,vexc);
%
% Purpose:
%    Compute the exchange correlation energy
%
% Input:
%  mol  --- Molecule object
%  rho  --- input charge density (3D array or cell array)
%  vexc --- exchange correlation potential (3D array or cell array)
%
% Output:
%  Exc  --- Exchange correlation energy (scalar)
%
nspin = mol.nspin;
hhh = mol.vol/mol.n1/mol.n2/mol.n3;

if nspin == 1
    arho = abs(rho);
elseif nspin == 2
    arho = abs(rho{1}+rho{2});
elseif nspin==4 && ~mol.domag
    arho = abs(rho{1});
end
% For spin-restricted, spin-unrestricted and spin-noncollinear
% calculations without magnetization, the Exc is calculated by 
% Vxc and rho directly. However, the calculation procedures are
% complicated for noncollinear spin density so we get Exc with
% Vxc together in function getVxc
if ~mol.noncolin||(mol.noncolin && ~mol.domag)
    if nargin < 4
        Exc = sumel(uxc2.*arho)*hhh;
    else
        Exc = sumel(uxc2.*arho-0.25*uxcsr)*hhh;
    end
else
    Exc = uxc2*hhh;
end

