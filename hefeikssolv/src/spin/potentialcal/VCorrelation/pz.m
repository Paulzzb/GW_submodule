function [vc,uc] = pz(rs,iflag)
% LDA parametrization from Monte Carlo Data:
% iflag = 1: J.P. Perdew and A. Zunger  
% iflag = 2: G. Ortiz and P. Ballone
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

