function occ = wgauss(ev,efermi,Tbeta,smear)
% gauss smearing method to calculate occupation rate

maxarg = 200;
c = sqrt(1/2);
x = (efermi-ev) / Tbeta;
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
