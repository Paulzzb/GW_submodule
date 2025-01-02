function [vxc,exc,rho] = VxcPZ_nc(rho)
% PZ functional for noncolinear spin 
% note : the output exc is the exchange-correlation energy
% LDA functional with spin
    e2 = e2Def();
    
    vxc = cell(4,1);
    for i = 1:4
        vxc{i} = zeros(size(rho{1}));
    end
    uxc = zeros(size(rho{1}));
       
    arho = abs(rho{1});
    idxxc = arho > 1e-10;
    arho = arho(idxxc);
    rs   = ((3/(4*pi))./arho).^(1/3);
    amag = sqrt(rho{2}.^2+rho{3}.^2+rho{4}.^2);
    zeta = amag(idxxc)./arho;
    id = abs(zeta)>1;
    zeta(id) = sign(zeta(id));
    
    [vx,ux] = VExchange_sla_spin(arho,zeta);
    [vc,uc] = VCorrelation_pz_spin(rs,zeta);
    
    idmag = amag > 1e-20;
    idxmag = amag(idxxc) > 1e-20;
    vxc{1}(idxxc) = e2*(vx{1} + vc{1} + vx{2} + vc{2})/2;
    vs = (vx{1} + vc{1} - vx{2} - vc{2})/2;
    for i = 2:4
        vxc{i}(idmag&idxxc) = e2*vs(idxmag).*rho{i}(idmag&idxxc)./amag(idmag&idxxc);
    end
    uxc(idxxc)    = e2*(ux+uc); 
    exc = sum(uxc.*abs(rho{1}),'all');
end