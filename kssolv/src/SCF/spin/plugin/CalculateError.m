function rhoerr = CalculateError(mol,rho,rhoin)
% CALCULATEERROR Calculate the SCF convergence residual/error.
%    rhoerr = CalculateError(mol,rho,rhoin) calculate the density/potential
%    residual/error between two SCF iterations. 
%    the defination of error is the integral of difference of rho/pot with
%    coulomb operator 1/r and calculated in G space.
%
%    See also transform_rho, rho_ddot.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if mol.nspin == 1
    rho = transform_rho(mol,rho,'real2rec');
    rhoin = transform_rho(mol,rhoin,'real2rec');
    rho = rho - rhoin;
    rhoerr = rho_ddot(mol,rho,rho);
elseif mol.nspin == 2
    rho = transform_rho(mol,rho,'updw2sum');
    rhoin = transform_rho(mol,rhoin,'updw2sum');
    rho = transform_rho(mol,rho,'real2rec');
    rhoin = transform_rho(mol,rhoin,'real2rec');
    rho{1} = rho{1} - rhoin{1};
    rho{2} = rho{2} - rhoin{2};
    rhoerr = rho_ddot(mol,rho,rho);
elseif mol.nspin == 4
    rho = transform_rho(mol,rho,'real2rec');
    rhoin = transform_rho(mol,rhoin,'real2rec');
    for i = 1:4
        rho{i} = rho{i} - rhoin{i};
    end
    rhoerr = rho_ddot(mol,rho,rho);
end
