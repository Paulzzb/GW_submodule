function [v1gcc,v2gcc,egcc] = VGCCorrelation_pbc(rho,grho2)
% VGCCORRELATION_PBC PBE correlation gradient correction for spin-restricted systems.
%    [v1gcc,v2gcc,egcc] = VGCCORRELATION_PBC(rho,grho2) returns the PBE
%    correlation gradient correction of the rho and the gradient square of rho.
%
%   See also exRef.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

% ga = 0.031091;
ga = 0.03109069086965; % ga=(1-log(2))/pi^2
% be = 0.066725;
be = 0.06672455060314922;

pi34=0.6203504908994; % pi34=(3/(4*pi))^(1/3)
xkf=1.919158292677513; xks=1.128379167095513; % xkf=(9*pi/4)^(1/3); xks=sqrt(4/pi)
rs = pi34 ./ rho.^(1/3);

[vc,uc] = VCorrelation_pw(rs);
ks = xks.*sqrt(xkf./rs);

t2 = grho2./(2*ks.*rho).^2;
expe = exp(-uc/ga);
af = be/ga*(1./(expe-1));
bf = expe.*(vc-uc);
y = af.*t2;
xy = (1+y)./(1+y+y.*y);
qy = y.*y.*(2+y)./(1+y+y.*y).^2;
s1 = 1+be/ga*t2.*xy;
h0 = ga*log(s1);
dh0 = be*t2./s1.*(-7/3*xy-qy.*(af.*bf/be-7/3));
ddh0 = be./(2*ks.*ks.*rho).*(xy-qy)./s1;

egcc = rho.*h0;
v1gcc = h0 + dh0;
v2gcc = ddh0;

end