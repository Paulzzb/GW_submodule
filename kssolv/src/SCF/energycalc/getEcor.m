function Ecor = getEcor(mol, rho, vtot, vion, vext)
% GETECOR Computes an energy correction term.
%    Ecor = getEcor(mol, rho, vtot, vion, vext) calculate the
%    correction energy which originates from the density and
%    potential difference between two SCF iteration
%
%    mol  --- a Molecule object
%    rho  --- charge density (3D array or cell array)
%    vtot --- total local potential (3D array or cell array)
%    vion --- local ionic potential (3D array or cell array)
%    vext --- external potential (3D array or cell array)

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

nspin = mol.nspin;
hhh = mol.vol/mol.n1/mol.n2/mol.n3;

if nspin == 1
    Ecor = sumel((vion+vext-vtot).*rho);
elseif nspin == 2
    Ecor = sumel((vion+vext-vtot{1}).*rho{1})...
         + sumel((vion+vext-vtot{2}).*rho{2});
elseif mol.nspin == 4
    if mol.domag
        Ecor = sumel((vion+vext-vtot{1}).*rho{1});
        for i = 2:4
            Ecor = Ecor + sumel((vext-vtot{i}).*rho{i});
        end
    else
        Ecor = sumel((vion+vext-vtot).*rho{1});
    end
end
Ecor = real(Ecor)*hhh;

end
