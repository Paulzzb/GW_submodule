function [vc,uc] = VCorrelation_pz_spin(rs,zeta)
% VCORRELATION_PZ_SPIN Perdew-Zunger correlation for spin-polarized systems.
%    [vc,uc] = VCORRELATION_PZ_SPIN(rs,zeta) returns the Perdew-Zunger correlation of
%    the rs = (3/4/pi/rho)^(1/3) and zeta (spin polarization) for spin-polarized systems.
%
%   See also exRef.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

p43 = 4.0/3.0;
third = 1.0/3.0;
vc = cell(2,1);

[vcu,ucu] = pz(rs,1); % unpolarized part
[vcp,ucp] = pz_polarized(rs); % polarized part

fz = ((1.0+zeta).^p43 + (1.0-zeta).^p43 - 2.0)/(2.0^p43-2.0);
dfz = p43*((1.0+zeta).^third-(1.0-zeta).^third)/(2.0^p43-2.0);
  
uc = ucu + fz .* (ucp - ucu);
vc{1} = vcu + fz .* (vcp - vcu) + (ucp - ucu) .* dfz .* ( 1.0 - zeta);
vc{2} = vcu + fz .* (vcp - vcu) + (ucp - ucu) .* dfz .* (-1.0 - zeta);

end

function [vc,uc] = pz(rs,iflag)
% PZ Perdew-Zunger/Ortiz-Ballone correlation for spin-restricted systems
%    [vc,uc] = PZ(rs,iflag) calculate the PZ or OB local correlation
%    which only depends on the total electron density.
%
%    LDA parametrization from Monte Carlo Data:
%    iflag = 1: J.P. Perdew and A. Zunger  
%    iflag = 2: G. Ortiz and P. Ballone

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

a = [0.0311, 0.031091];
b = [-0.048, -0.046644];
c = [0.0020, 0.00419];
d = [-0.0116,-0.00983];
gc = [-0.1423, -0.103756];
b1 = [1.0529, 0.56371]; 
b2 = [0.3334,  0.27358];  

uc = zeros(size(rs));
vc = zeros(size(rs));

idxl = rs < 1; 
rsl  = rs(idxl);
lnrs = log(rsl);
uc(idxl) = a(iflag)*lnrs + b(iflag) + c(iflag)*rsl.*lnrs + d(iflag)*rsl;
vc(idxl) = a(iflag)*lnrs + ( b(iflag) - a(iflag)/3.0 ) + 2.0/3.0 * ...
  c(iflag)*rsl.*lnrs + ( 2.0*d(iflag) - c(iflag) )/3.0*rsl;

idxg = rs >= 1;
rsg  = rs(idxg);
rs12 = sqrt(rsg);  
ox  = 1.0 + b1(iflag)*rs12 + b2(iflag)*rsg;
dox = 1.0 + 7.0/6.0*b1(iflag)*rs12 + 4.0/3.0*b2(iflag)*rsg;    
uc(idxg) = gc(iflag)./ox;
vc(idxg) = uc(idxg).*dox./ox;

end

function [vc,uc] = pz_polarized(rs)
% PZ_POLARIZED Unpolarized part of Perdew-Zunger correlation for spin-polarized systems.
%    [vc,uc] = PZ_POLARIZED(rs,iflag) calculate the PZ local correlation
%    which depends on the spin-up and spin-down electron density.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

% J.P. Perdew and A. Zunger
% spin_polarized energy and potential. 
a = 0.01555;
b = -0.0269;
c = 0.0007;
d = -0.0048;
gc = -0.0843;
b1 = 1.3981;
b2 = 0.2611;

uc = zeros(size(rs));
vc = zeros(size(rs));

idxl = rs < 1;
rsl  = rs(idxl);
lnrs = log(rsl);
uc(idxl) = a*lnrs + b + c*rsl.*lnrs + d*rsl;
vc(idxl) = a*lnrs + (b - a/3) + 2/3*c*rsl.*lnrs + (2*d - c)/3*rsl;

idxg = rs >= 1;
rsg  = rs(idxg);
rs12 = sqrt(rsg);
ox = 1 + b1*rs12 + b2*rsg;
dox = 1 + 7/6*b1*rs12 + 4/3*b2*rsg;
uc(idxg) = gc./ox;
vc(idxg) = uc(idxg).*dox./ox;

end