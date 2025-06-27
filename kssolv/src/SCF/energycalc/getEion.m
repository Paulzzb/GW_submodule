function Eion = getEion(mol,rho,vion)
% GETEION Compute the ionic potential energy.
%    Eion = getEion(rho, vion) get the ionic potential energy
%    by ionic potential which is described by pseudopotential.
%
%    mol  --- a Molecule object
%    rho  --- input charge density (3D array)
%    vion --- ionic potential (3D array)

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

hhh = mol.vol/mol.n1/mol.n2/mol.n3;
Eion      = sumel(rho.*vion);
Eion      = real(Eion)*hhh;

end
