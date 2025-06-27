function occ = wgauss(ev,efermi,Tbeta,smear)
% WGAUSS returns the occupation number of electrons
%    occ = wgauss(ev,efermi,Tbeta,smear) calculate the fractional occupation
%    number according to eigenvalues ev, fermi energy efermi, temperature Tbeta
%    and smearing method smear
%
%    See also getocc.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

% FIXME The results obtained with MV and MP method are not consistant with QE for Cu4.
% Check it and fix this problem !

maxarg = 200;
c = sqrt(1/2);
x = (efermi-ev)/Tbeta;
if strcmp(smear,'fd')||strcmp(smear,'fermi-dirac')
    id1 = (x < -maxarg);
    id2 = (x > maxarg);
    occ = 1.0 ./ (1.0 + exp ( - x) );
    occ(id1) = 0;
    occ(id2) = 1;
elseif strcmp(smear,'cold')||strcmp(smear,'mv')   
    xp = x - 1.0 / sqrt (2.0);
    arg = min (maxarg, xp.^2);
    occ = 0.5 * erf(xp) + 1.0 / sqrt (2.0 * pi) * exp ( - ...
          arg) + 0.5;
elseif strcmp(smear,'gauss')||strcmp(smear,'gaussian')
    occ = 0.5 * erfc( - x * sqrt (2.0) * c);
elseif strcmp(smear,'mp')
    occ = 0.5 * erfc( - x * sqrt (2.0) * c);
    hd = 0.0;
    arg = min (maxarg, x.*2);
    hp = exp(-arg);
    ni = 0;
    a = 1.0 / sqrt (pi);
    for k = 1
        hd = 2.0 * x .* hp - 2.0 * ni .* hd;
        ni = ni + 1;
        a = - a / (k * 4.0);
        occ = occ - a * hd;
        hp = 2.0 * x .* hd-2.0 * ni .* hp;
        ni = ni + 1;
    end
else
    fprintf('%s method has not been achieved...\n',smear);
end

