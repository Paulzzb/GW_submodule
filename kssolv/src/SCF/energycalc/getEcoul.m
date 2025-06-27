function Ecoul = getEcoul(mol,rho,vhart)
% GETECOUL Compute the Coulomb potential energy induced by vhart
%    Ecoul = getEcoul(mol, rho, vhart) calculate the Coulomb repulsion
%    energy by density and Hartree potential
%
%    mol    --- a Molecule object
%    rho    --- input charge density (3D array or cell array)
%    vhart  --- Hartee (Coulomb) potential (3D array)

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

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

end
