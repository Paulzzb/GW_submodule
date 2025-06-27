function Exc = getExc(mol,exc)
% GETEXC Compute the exchange correlation energy.
%    Exc = getExc(mol,exc) calculates the exchange correlation energy
%    by local energy per unit volume (exc).

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

hhh = mol.vol/mol.n1/mol.n2/mol.n3;

if ~mol.noncolin||(mol.noncolin && ~mol.domag)
    Exc = sumel(exc)*hhh;
else
    Exc = exc*hhh;
end

end
