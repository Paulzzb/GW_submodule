function [vhart,vxc,exc,rho] = getVhxc(mol,rho)
% GETVHXC Computes the potential depending on electron density.
%    [vhart,vxc,exc,rho] = getVhxc(mol,rho) get the Hartree
%    potential, exchange correlation potential and exchange 
%    correlation energy density by electron density.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

% computing the potential from the Hartree approximation
vhart = getVhart(mol,rho);

% computing the exchange correlation, for rho and the molecule
[vxc,exc,rho] = getVxc(mol,rho);

end
