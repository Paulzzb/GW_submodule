function Eext = getEext(mol,rho,vext)
% GETEEXT Compute the potential energy induced by an external potential vext.
%    Eext = getEext(rho, vext) calculate the potential energy
%    by the external scalar potential and electron density.
%
%    mol  --- Molecule object 
%    rho  --- input charge density (3D array or cell array)
%    vext --- external potential (3D array or cell array)

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

hhh = mol.vol/mol.n1/mol.n2/mol.n3;
Eext      = sumel(rho.*vext);
Eext      = real(Eext)*hhh;

end
