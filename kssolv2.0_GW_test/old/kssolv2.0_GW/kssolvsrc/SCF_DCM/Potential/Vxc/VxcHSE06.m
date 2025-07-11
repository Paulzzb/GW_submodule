function [vxc,uxc,rho,uxcsr] = VxcHSE06(mol,rho)

e2 = e2Def();

[grho,grho2] = getgradrho(mol,rho);

vxc = zeros(size(rho));
uxc = zeros(size(rho));
uxcsr = zeros(size(rho));
h2xc = zeros(size(rho));
    
idxxc = abs(rho) > 1e-10;
rhoxc = abs(rho(idxxc));

rs = ((3/(4*pi))./rhoxc).^(1/3);
[vx,ux] = VExchange_sla(rs);
[vc,uc] = VCorrelation_pw(rs);

vxc(idxxc) = e2*(vx+vc);
uxc(idxxc) = e2*(ux+uc);

idxcxc = abs(rho) > 1e-6 & grho2 > 1e-10;
if nnz(idxcxc)
    rhocxc = abs(rho(idxcxc));
    grho2cxc = grho2(idxcxc);
    
    [v1gcx,v2gcx,ugcx] = VGCExchange_pbx(rhocxc,grho2cxc);
    [uxcsr(idxcxc),v1gcxsr,v2gcxsr] = VGCExchange_pbxsr(rhocxc,grho2cxc);
    v1gcx = v1gcx - 0.25*v1gcxsr;
    v2gcx = v2gcx - 0.25*v2gcxsr;
    [v1gcc,v2gcc,ugcc] = VGCCorrelation_pbc(rhocxc,grho2cxc);

    vxc(idxcxc) = vxc(idxcxc) + e2*(v1gcx+v1gcc);
    uxc(idxcxc) = uxc(idxcxc) + e2*(ugcx+ugcc);
    h2xc(idxcxc) = e2*(v2gcx+v2gcc);
end

h2xcgrho2 = repmat((h2xc),1,1,1,3).*grho;

n1 = mol.n1;
n2 = mol.n2;
n3 = mol.n3;
ecut2 = mol.ecutrho;
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
