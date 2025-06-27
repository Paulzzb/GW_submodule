function [vxc,exc,rho,h2xc] = VxcHSE06(mol,rho,grho2)
% VXCHSE06 HSE06 exchange correlation for spin-restricted systems.
%    [vxc,exc,rho,h2xc] = VXCHSE06(mol,rho,grho2) returns the HSE06 
%    exchange correlation of the one-component rho for spin-restricted systems.
%
%   See also exRef.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

alpha = 1.0 - mol.exx_lr;
beta = mol.exx_lr - mol.exx_sr;
gamma = mol.dft_c;

e2 = e2Def();

if nargin < 3 % Skipped when using getFxc_FD.m, where grho2 is given.
    [grho,grho2] = getgradrho(mol,rho);
end

vxc = zeros(size(rho));
exc = zeros(size(rho));
h2xc = zeros(size(rho));

idxxc = abs(rho) > 1e-10;
rhoxc = abs(rho(idxxc));

rs = ((3/(4*pi))./rhoxc).^(1/3);
[vx,ux] = VExchange_sla(rs);
[vc,uc] = VCorrelation_pw(rs);

vxc(idxxc) = e2*(vx*alpha+vc*gamma);
exc(idxxc) = e2*(ux*alpha+uc*gamma).*rhoxc;

idxcxc = abs(rho) > 1e-6 & grho2 > 1e-10;
rhocxc = abs(rho(idxcxc));
grho2cxc = grho2(idxcxc);

[v1gcx,v2gcx,egcx] = VGCExchange_pbx(rhocxc,grho2cxc);

if beta ~= 0
    [v1gcxsr,v2gcxsr,egcxsr] = VGCExchange_pbxsr(rhocxc,grho2cxc,mol.omega);
    v1gcx = v1gcx*alpha + v1gcxsr*beta;
    v2gcx = v2gcx*alpha + v2gcxsr*beta;
    egcx = egcx*alpha + egcxsr*beta;
else
    v1gcx = v1gcx*alpha;
    v2gcx = v2gcx*alpha;
    egcx = egcx*alpha;
end

[v1gcc,v2gcc,egcc] = VGCCorrelation_pbc(rhocxc,grho2cxc);

vxc(idxcxc) = vxc(idxcxc) + e2*(v1gcx+v1gcc*gamma);
exc(idxcxc) = exc(idxcxc) + e2*(egcx+egcc*gamma);
h2xc(idxcxc) = e2*(v2gcx+v2gcc*gamma);

if nargin < 3 % Skipped when using getFxc_FD.m, where 'vxc without v2xc' and 'h2xc' are demanded.
    h2xcgrho2 = repmat((h2xc),1,1,1,3).*grho;
    
    n1 = mol.n1;
    n2 = mol.n2;
    n3 = mol.n3;
    ecut2 = mol.ecut2;
    grid2 = Ggrid(mol,ecut2);
    idxnz2 = grid2.idxnz;
    gkx2 = grid2.gkx;
    gky2 = grid2.gky;
    gkz2 = grid2.gkz;
    hxc2g = reshape(fft3(h2xcgrho2),n1*n2*n3,3);
    gaux = 1i*(hxc2g(idxnz2,1).*gkx2 ...
        + hxc2g(idxnz2,2).*gky2 + hxc2g(idxnz2,3).*gkz2);
    gauxfull = zeros(n1*n2*n3,1);
    gauxfull(idxnz2) = gaux;
    v2xc = real(ifft3(reshape(gauxfull,n1,n2,n3)));
    
    vxc = vxc-v2xc;
end

end
