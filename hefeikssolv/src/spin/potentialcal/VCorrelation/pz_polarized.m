function [vc,uc] = pz_polarized(rs)
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

