function [vxc,exc,rho] = VxcHSE06_nc(mol,rho)
% HSE06 functional for noncolinear spin
% note : the output exc is the exchange-correlation energy
% LDA part
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
    vc = cell(2,1);
    [uc,vc{1},vc{2}] = VCorrelation_pw_spin(rs,zeta);
    
    idmag = amag > 1e-20;
    idxmag = amag(idxxc) > 1e-20;
    vxc{1}(idxxc) = e2*(vx{1} + vc{1} + vx{2} + vc{2})/2;
    vs = (vx{1} + vc{1} - vx{2} - vc{2})/2;
    for i = 2:4
        vxc{i}(idmag&idxxc) = e2*vs(idxmag).*rho{i}(idmag&idxxc)./amag(idmag&idxxc);
    end
    uxc(idxxc)    = e2*(ux+uc);
    exc = sum(uxc.*rho{1},'all');
% GGA part
    [rhoaux,segni] = rho_compute(rho,mol.lsign,mol.ux);
    [grho_up,~] = getgradrho(mol,rhoaux{1});
    [grho_dw,~] = getgradrho(mol,rhoaux{2});
    grho = cell(2,1);
    grho{1} = grho_up;
    grho{2} = grho_dw;
    
    [ex,ec,v1x,v2x,v1c,v2c,v2c_ud] = Vxc_gcx(rhoaux,grho,1);
    vgxc = cell(2,1);
    vgxc{1} = e2*(v1x{1}+v1c{1});
    vgxc{2} = e2*(v1x{2}+v1c{2});
    
    h_up = cell(3,1);
    h_dw = cell(3,1);
    for ipol=1:3
        grup = grho_up(:,:,:,ipol);
        grdw = grho_dw(:,:,:,ipol);
        h_up{ipol} = e2*((v2x{1}+v2c{1}).*grup+v2c_ud.*grdw);
        h_dw{ipol} = e2*((v2x{2}+v2c{2}).*grdw+v2c_ud.*grup);
    end
    
    dh_up = graddot(mol,h_up);
    dh_dw = graddot(mol,h_dw);
    vgxc{1} = vgxc{1} - dh_up;
    vgxc{2} = vgxc{2} - dh_dw;
    exc = exc + sum(ex,'all') + sum(ec,'all');
    
    vxc{1} = vxc{1} + 0.5*(vgxc{1}+vgxc{2});
    idmag = amag > 1e-12;
    dvgxc = vgxc{1}(idmag) - vgxc{2}(idmag);
    for i = 2:4      
        vxc{i}(idmag) = vxc{i}(idmag) + segni(idmag).*dvgxc.*rho{i}(idmag)./amag(idmag)/2;
    end   
end