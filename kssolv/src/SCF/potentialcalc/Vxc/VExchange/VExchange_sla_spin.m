function [vx,ux] = VExchange_sla_spin(rs,zeta)
% VEXCHANGE_SLA_SPIN Slater exchange with alpha=2/3 for spin-polarized systems.
%    [vx,ux] = VEXCHANGE_SLA_SPIN(rs,zeta) returns the Slater exchange with alpha=2/3 of
%    the rs = (3/4/pi/rho)^(1/3) and zeta (spin polarization) for spin-polarized systems.
%
%   See also exRef.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

% alpha = 2/3;
% f = -9/8*(3/2/pi)^(2/3);
% falpha = f*alpha;
falpha = -0.458165293283143;
third = 1.0/3.0;
p43 = 4.0/3.0;
vx = cell(2,1);

rhoup = (1.0+zeta).^third./rs;
uxup = falpha.*rhoup;
vx{1} = p43*falpha.*rhoup;
  
rhodw = (1.0-zeta).^third./rs;
uxdw = falpha.*rhodw;
vx{2} = p43*falpha.*rhodw;
  
ux = 0.5*((1.0+zeta).*uxup+(1.0-zeta).*uxdw);
  
end

