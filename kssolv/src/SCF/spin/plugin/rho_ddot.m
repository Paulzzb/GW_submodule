function rhoerr = rho_ddot(mol,rho1,rho2)
% RHO_DDOT Calculate function integral with Coulomb operator in G space.
%    rhoerr = rho_ddot(mol,rho1,rho2) calculate integral 
%    \int rho1(r)rho2(r)/r dr in G space 
%
%    this function is used to estimate the self-consistency error on the energy.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if isempty(rho1)||isempty(rho2)
    rhoerr = 0;
else
    grid2 = Ggrid(mol,mol.ecut2);
    gkk = grid2.gkk;
    if mol.nspin ~= 1
        rhoerr = 4*pi*sum(real(conj(rho1{1}(2:end)).*rho2{1}(2:end)./gkk(2:end)),'all');
    else
        rhoerr = 4*pi*sum(real(conj(rho1(2:end)).*rho2(2:end)./gkk(2:end)),'all');
    end

    if mol.nspin ~= 1
        for i = 2:mol.nspin
            fac = 4*pi/(2*pi)^2;
            rhoerr = rhoerr + fac*sum(real(conj(rho1{i}).*rho2{i}),'all');
        end
    end
    rhoerr = rhoerr*mol.vol/2;
end
    
    
    
    
