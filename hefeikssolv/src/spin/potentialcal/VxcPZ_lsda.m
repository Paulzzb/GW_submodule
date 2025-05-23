function [vxc,uxc,rho] = VxcPZ_lsda(rho)
% PZ functional for LSDA
% LDA functional with spin
    e2 = e2Def();
    
    vxc = cell(2,1);
    vxc{1} = zeros(size(rho{1}));
    vxc{2} = zeros(size(rho{1}));
    uxc = zeros(size(rho{1}));
       
    arho = abs(rho{1} + rho{2});
    idxxc = arho > 1e-10;
    arho = arho(idxxc);
    rs   = ((3/(4*pi))./arho).^(1/3);
    zeta = (rho{1}(idxxc) - rho{2}(idxxc))./arho;
    id = abs(zeta)>1;
    zeta(id) = sign(zeta(id));
    
    [vx,ux] = VExchange_sla_spin(arho,zeta);
    [vc,uc] = VCorrelation_pz_spin(rs,zeta);

    vxc{1}(idxxc) = e2*(vx{1} + vc{1});
    vxc{2}(idxxc) = e2*(vx{2} + vc{2});
    uxc(idxxc)    = e2*(ux+uc);
           
end