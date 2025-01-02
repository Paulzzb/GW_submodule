function fxcPZ = getfxcPZ(rho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Exchange functional: slater 
%         Ex = -3/4 * (3/pi*rho)^(4/3) 
%   Correlation functional: PZ
%         Ec = A*ln(rs) + B + rs*(C*ln(rs)+D)   for rs < 1
%         Ec = gc / ( 1 + b1/rs^1/2 + b2/rs)    otherwise 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e2 = e2Def();
falpha = -0.458165293283143;

idxxc = rho > 1e-15;
rhoxc = rho(idxxc);

rs = ((3/(4*pi))./rhoxc).^(1/3);
drsdrho = -rs.^4 * 4*pi/9;

% 1. Exchange part
dvxdrho = 16 * pi / 27 * falpha * rs.^2;

% 2. Correlation part
a  =  0.0311;
b  = -0.048;
c  =  0.0020;
d  = -0.0116;
gc = -0.1423;
b1 =  1.0529;
b2 =  0.3334;
dvcdrho = zeros(size(rs));

% high density formula
idxl = rs < 1;
rsl  = rs(idxl);
lnrs = log (rsl);
dvcdrho(idxl) = a ./rsl + 2/3*c*(lnrs + 1) + (2*d-c)/3;

% interpolation formula
idxg = rs >= 1;
rsg  = rs(idxg);
rs12 = sqrt(rsg);
ox   = 1 + b1*rs12 + b2*rsg;
dox  = 1 + 7/6*b1*rs12 + 4/3*b2*rsg;
dvcdrho(idxg) =  gc ./(ox.^2) ...
    .* ( 7/12*b1 ./ rs12 + 4/3*b2 - dox./ox .*(b1 ./rs12 + 2*b2) );

% 3. Sum up exchange and correlation parts
fxcPZtmp = e2*(dvcdrho .*drsdrho + dvxdrho);
%fxcPZtmp = e2*(dvxdrho);
%fxcPZtmp = e2*(dvcdrho .*drsdrho);
fxcPZ = zeros(size(rho));
fxcPZ(idxxc) = fxcPZtmp;
