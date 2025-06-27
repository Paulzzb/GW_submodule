function [vxc,exc,rho] = VxcPW_lsda(rho)
% VXCPW_LSDA PW exchange correlation for spin-polarized systems.
%    [vxc,exc,rho] = VXCPW_LSDA(rho) returns the PW exchange correlation of the
%    two-component rho for spin-polarized systems.
%
%   See also exRef.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

% PW functional for LSDA
% LDA functional with spin
    e2 = e2Def();
    
    vxc = cell(2,1);
    vxc{1} = zeros(size(rho{1}));
    vxc{2} = zeros(size(rho{1}));
    exc = zeros(size(rho{1}));
       
    arho = abs(rho{1} + rho{2});
    idxxc = arho > 1e-10;
    arho = arho(idxxc);
    rs   = ((3/(4*pi))./arho).^(1/3);
    zeta = (rho{1}(idxxc) - rho{2}(idxxc))./arho;
    id = abs(zeta)>1;
    zeta(id) = sign(zeta(id));
    
    [vx,ux] = VExchange_sla_spin(rs,zeta);
    [vc,uc] = VCorrelation_pw_spin(rs,zeta);

    vxc{1}(idxxc) = e2*(vx{1} + vc{1});
    vxc{2}(idxxc) = e2*(vx{2} + vc{2});
    exc(idxxc)    = e2*(ux+uc).*arho;
           
end
