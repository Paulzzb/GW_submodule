function [vc,uc] = VCorrelation_pw_spin(rs,zeta)
% VCORRELATION_PW_SPIN Perdew-Wang correlation for spin-polarized systems.
%    [vc,uc] = VCORRELATION_PW_SPIN(rs,zeta) returns the Perdew-Wang correlation of
%    the rs = (3/4/pi/rho)^(1/3) and zeta (spin polarization) for spin-polarized systems.
%
%   See also exRef.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.
 
% xc parameters, unpolarised
a = 0.031091; a1 = 0.21370;
b1 = 7.5957; b2 = 3.5876; b3 = 1.6382;b4 = 0.49294;

% xc parameters, polarised
ap = 0.015545; a1p = 0.20548; 
b1p = 14.1189; b2p = 6.1977; b3p = 3.3662;b4p = 0.62517;

% xc PARAMETERs, antiferro
aa = 0.016887; a1a = 0.11125;
b1a = 10.357; b2a = 3.6231; b3a = 0.88026; b4a = 0.49671;

fz0 = 1.709921;

zeta2 = zeta .* zeta;
zeta3 = zeta2 .* zeta;
zeta4 = zeta3 .* zeta;
rs12 = sqrt(rs);
rs32 = rs .* rs12;
rs2 = rs.^2;

% unpolarised
om = 2.0 * a * (b1 * rs12 + b2 * rs + b3 * rs32 + b4 * rs2);
dom = 2.0 * a * (0.50 * b1 * rs12 + b2 * rs + 1.5 * b3 * rs32 ...
   + 2.0 * b4 * rs2);
olog = log(1.0 + 1.0 ./ om);
upwc = - 2.0 * a * (1.0 + a1 * rs) .* olog;
vpwc = - 2.0 * a * (1.0 + 2.0 / 3.0 * a1 * rs) .* olog - 2.0/ ...
     3.0 * a * (1.0 + a1 * rs) .* dom ./ (om .* (om + 1.0) );

% polarized
omp  = 2.0 * ap * (b1p * rs12 + b2p * rs + b3p * rs32 + b4p * rs2);
domp = 2.0 * ap * (0.5 * b1p * rs12 + b2p * rs + 1.50 * b3p * ...
     rs32 + 2.0 * b4p * rs2);
ologp = log(1.0 + 1.0 ./ omp);
upwcp = - 2.0 * ap * (1.0 + a1p * rs) .* ologp;
vpwcp = - 2.0 * ap * (1.0 + 2.0 / 3.0 * a1p * rs) .* ologp - ...
      2.0 / 3.0 * ap * (1.0 + a1p * rs) .* domp ./ (omp .* (omp + 1.0));

% antiferro
oma = 2.0 * aa * (b1a * rs12 + b2a * rs + b3a * rs32 + b4a * rs2);
doma = 2.0 * aa * ( 0.5 * b1a * rs12 + b2a * rs + 1.5 * b3a * ...
     rs32 + 2.0 * b4a * rs2 );
ologa = log( 1.0 + 1.0./oma );
alpha = 2.0 * aa * (1.0 + a1a*rs) .* ologa;
vpwca = 2.0 * aa * (1.0 + 2.0/3.0 * a1a .* rs) .* ologa + ...
      2.0 / 3.0 * aa .* (1.0 + a1a*rs) .* doma ./ (oma .* (oma + 1.0));

fz = ( (1.0 + zeta).^(4.0 / 3.0) + (1.0 - zeta).^(4.0 / ...
      3.0) - 2.0) ./ (2.0^ (4.0 / 3.0) - 2.0);
dfz = ( (1.0 + zeta).^(1.0 / 3.0) - (1.0 - zeta).^(1.0 / ...
      3.0) ) * 4.0 / (3.0 * (2.0^ (4.0 / 3.0) - 2.0) );

uc = upwc + alpha .* fz .* (1.0 - zeta4) ./ fz0 + (upwcp - upwc) ...
          .* fz .* zeta4;

vc = cell(2,1);

vc{1} = vpwc + vpwca .* fz .* (1.0 - zeta4) ./ fz0 + (vpwcp - vpwc) ...
             .* fz .* zeta4 + (alpha ./ fz0 .* (dfz .* (1.0 - zeta4) ...
             - 4.0 * fz .* zeta3) + (upwcp - upwc) .* (dfz .* zeta4 +... 
             4.0 .* fz .* zeta3) ) .* (1.0 - zeta); % vc_up

vc{2} = vpwc + vpwca .* fz .* (1.0 - zeta4) ./ fz0 + (vpwcp - vpwc) ...
             .* fz .* zeta4 - (alpha ./ fz0 .* (dfz .* (1.0 - zeta4)... 
             - 4.0 * fz .* zeta3) + (upwcp - upwc) .* (dfz .* zeta4 +... 
             4.0 .* fz .* zeta3) ) .* (1.0 + zeta); % vc_dw
end
