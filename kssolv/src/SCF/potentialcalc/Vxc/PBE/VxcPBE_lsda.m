function [vxc,exc,rho] = VxcPBE_lsda(mol,rho)
% VXCPBE_LSDA PBE exchange correlation for spin-polarized systems.
%    [vxc,exc,rho] = VXCPBE_LSDA(mol,rho) returns the PBE exchange correlation
%    of the two-component rho for spin-polarized systems.
%
%   See also exRef.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

% PBE functional for LSDA   
    e2 = e2Def();
    
    vxc = cell(2,1);
    vxc{1} = zeros(mol.n1,mol.n2,mol.n3);
    vxc{2} = zeros(mol.n1,mol.n2,mol.n3);
    exc = zeros(mol.n1,mol.n2,mol.n3);
% Do LDSA calculation first      
    arho = abs(rho{1} + rho{2});
    idxxc = arho > 1e-10;
    arho = arho(idxxc);
    rs   = ((3/(4*pi))./arho).^(1/3);
    zeta = (rho{1}(idxxc) - rho{2}(idxxc))./arho;
    id = abs(zeta) > 1;
    zeta(id) = sign(zeta(id));
    
    [vx,ux] = VExchange_sla_spin(rs,zeta);
    [vc,uc] = VCorrelation_pw_spin(rs,zeta);

    vxc{1}(idxxc) = e2*(vx{1} + vc{1});
    vxc{2}(idxxc) = e2*(vx{2} + vc{2});
    exc(idxxc) = e2*(ux + uc).* arho;
% GGA correction by gradient of rho
    [grho_up,~] = getgradrho(mol,rho{1});
    [grho_dw,~] = getgradrho(mol,rho{2});
    grho = cell(2,1);
    grho{1} = grho_up;
    grho{2} = grho_dw;
    
    [ex,ec,v1x,v2x,v1c,v2c,v2c_ud] = Vxc_gcxc(rho,grho);
    vxc{1} = vxc{1} + e2*(v1x{1}+v1c{1});
    vxc{2} = vxc{2} + e2*(v1x{2}+v1c{2});
    
    h_up = cell(3,1);
    h_dw = cell(3,1);
    for ipol=1:3
        grup = grho_up(:,:,:,ipol);
        grdw = grho_dw(:,:,:,ipol);
        h_up{ipol} = e2*((v2x{1}+v2c{1}).*grup+v2c_ud.*grdw);
        h_dw{ipol} = e2*((v2x{2}+v2c{2}).*grdw+v2c_ud.*grup);
    end
    
    dh_up = graddot(mol,h_up);
    dh_dw = graddot(mol,h_dw);
    vxc{1} = vxc{1} - dh_up;
    vxc{2} = vxc{2} - dh_dw;
    exc = exc + ex + ec;
end

