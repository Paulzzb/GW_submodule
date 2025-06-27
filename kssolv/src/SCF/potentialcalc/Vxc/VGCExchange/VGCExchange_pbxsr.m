function [v1gcxsr,v2gcxsr,egcxsr] = VGCExchange_pbxsr(rho,grho2,omega)
% VGCEXCHANGE_PBXSR PBE short-range exchange gradient correction for spin-restricted systems.
%    [v1gcxsr,v2gcxsr,egcxsr] = VGCEXCHANGE_PBXSR(rho,grho2,omega) returns the PBE short-range
%    exchange gradient correction of the rho and the gradient square of rho, which already 
%    contains the Slater short-range exchange.
%
%   See also exRef.
%
%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

US = 0.161620459673995492;
AX = -0.738558766382022406;
f1 = -1.10783814957303361;
alpha = 2/3;

rs = rho.^(1/3);
vx = 4/3*f1*alpha.*rs;        % Slater exchange potential

AA    = grho2;
RR    = 1./(rho.*rs);
ex    = AX./RR;               % Slater exchange energy per unit volume (= ux.*rho)
s2    = US*US.*AA.*RR.*RR;

s = sqrt(s2);
s(s>8.3) = 8.572844 - 18.796223./s2(s>8.3);   % Gradient scaling for Lieb-Oxford bound

[Fx,D1X,D2X] = wpbe_analy_erfc_approx_grad(rho,s,omega);

egcxsr = ex.*Fx;
DSDN = -4/3.*s./rho;
v1gcxsr = vx.*Fx + (DSDN.*D2X+D1X).*ex;
DSDG = US.*RR;
v2gcxsr = ex./sqrt(AA).*DSDG.*D2X;

end
