function entropy = getEntropy(ev,efermi,Tbeta,smear)
% GETENTROPY calculates entropy of the statues.
%    entropy = GETENTROPY(ev,efermi,Tbeta,smear) calculates the entropy of the
%    occupation status and bete temperature. The calculation is very
%    standard, i.e., entropy = sum_i occ(i)log occ(i).
%
%   See also scf4m.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

x = (efermi-ev)/Tbeta;
%The entropy calculation method is from QE, w1gauss.f90.
%Tbeta is not included in entropy, it should be multiplied outside.
if Tbeta < eps
    entropy = 0;
elseif strcmp(smear,'fd')||strcmp(smear,'fermi-dirac')
    %Fermi-Dirac smearing.
    idx = (abs(x) <= 36);
    occ = 1./(1+exp(-x(idx)));
    entropy = sum(occ.*log(occ)+(1-occ).*log(1-occ));
elseif strcmp(smear,'cold')||strcmp(smear,'mv')
    % 
    xp = x - 1/sqrt(2.0);
    arg = min(200,x.^2);
    occ = 1/sqrt(2*pi)*xp.*exp(-arg);
    entropy = sum(occ);
elseif strcmp(smear,'gauss')||strcmp(smear,'gaussian')
    % Gaussian smearing
    arg = min(200,x.^2);
    occ = -0.5*exp(-arg)/sqrt(pi);
    entropy = sum(occ);
elseif strcmp(smear,'mp')
    % 
    arg = min(200,x.^2);
    occ = -0.5*exp(-arg)/sqrt(pi);
    hd = 0;
    hp = exp(-arg);
    ni = 0;
    a = 1/sqrt(pi);
    for k = 1
        hd = 2*x.*hp - 2*ni*hd;
        ni = ni + 1;
        hpm1 = hp;
        hp = 2*x.*hd - 2*ni*hp;
        ni = ni + 1;
        a = -a/(k*4);
        occ = occ - a*(0.5*hp + ni*hpm1);
    end
    entropy = sum(occ);
end
