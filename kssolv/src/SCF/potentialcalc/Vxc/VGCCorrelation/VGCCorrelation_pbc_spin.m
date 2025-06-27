function [v1gcc,v2gcc,egcc] = VGCCorrelation_pbc_spin(rho,zeta,grho2)
% VGCCORRELATION_PBC_SPIN PBE correlation gradient correction for spin-polarized systems.
%    [v1gcc,v2gcc,egcc] = VGCCORRELATION_PBC_SPIN(rho,zeta,grho2) returns the PBE
%    correlation gradient correction of the rho and the gradient square of rho.
%
%   See also exRef.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

ga = 0.031091;
% ga = 0.03109069086965; % ga=(1-log(2))/pi^2
% be = 0.066725;
be = 0.06672455060314922;

pi34=0.6203504908994; % pi34=(3/(4*pi))^(1/3)
xkf=1.919158292677513; xks=1.128379167095513; % xkf=(9*pi/4)^(1/3); xks=sqrt(4/pi)
rs = pi34 ./ rho.^(1/3);

[vc,uc] = VCorrelation_pw_spin(rs,zeta);
ks = xks.*sqrt(xkf./rs);

fz = 0.5*( (1+zeta).^(2/3) + (1-zeta).^(2/3) );
fz2 = fz .* fz;
fz3 = fz2 .* fz;

dfz = ( (1+zeta).^(-1/3) - (1 - zeta).^(-1/3) ) / 3;

t  = sqrt(grho2) ./ (2 * fz .* ks .* rho);
expe = exp( - uc ./ (fz3 * ga) );
af   = be / ga * (1 ./ (expe-1) );
bfup = expe .* (vc{1} - uc) ./ fz3;
bfdw = expe .* (vc{2} - uc) ./ fz3;

y  = af .* t .* t;
xy = (1 + y) ./ (1 + y + y .* y);
qy = y .* y .* (2 + y) ./ (1 + y + y .* y).^2;
s1 = 1 + be / ga .* t .* t .* xy;

h0 = fz3 .* ga .* log(s1);

dh0up = be * t .* t .* fz3 ./ s1 .* ( -7/3 * xy - qy .* ...
      (af .* bfup / be-7/3) );

dh0dw = be * t .* t .* fz3 ./ s1 .* ( -7/3 * xy - qy .* ...
      (af .* bfdw / be-7/3) );

dh0zup =   (3 .* h0 ./ fz - be * t .* t .* fz2 ./ s1 .*  ...
         (2 .* xy - qy .* (3 * af .* expe .* uc ./ fz3 ./ ...
         be+2) ) ) .* dfz .* (1 - zeta);

dh0zdw = - (3 .* h0 ./ fz - be * t .* t .* fz2 ./ s1 .*  ...
         (2 .* xy - qy .* (3 * af .* expe .* uc ./ fz3 ./ ...
         be+2) ) ) .* dfz .* (1 + zeta);

ddh0 = be .* fz ./ (2 * ks .* ks .* rho) .* (xy - qy) ./ s1;

egcc     = rho .* h0;

v1gcc = cell(2,1);
v1gcc{1} = h0 + dh0up + dh0zup;
v1gcc{2} = h0 + dh0dw + dh0zdw;

v2gcc    = ddh0; 
end