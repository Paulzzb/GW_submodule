function [ux,uc,v1x,v2x,v1c,v2c,v2c_ud] = Vxc_gcx(rho,grho,ishybrid)
% Graident part of GGA functional
% ishybrid = 0 PBE case: pbx + pbc
% ishybrid = 1 HSE06 case: pbx + pbc - 0.25*pbxsr

    grho2 = cell(2,1);
    v1c   = cell(2,1);
    v2c   = cell(2,1);
    for i = 1 : 2
        grho2{i} = sum(grho{i}.^2,4);
    end
    % PBX
    [ux,v1x,v2x] = VGCExchange_pbx_spin(rho,grho2,ishybrid);
    
    arho = rho{1} + rho{2};
    grho_tmp = grho{1} + grho{2};
    grho2{1} = sum(grho_tmp.^2,4);
    
    uc = zeros(size(rho{1}));
    for i = 1:2
       v1c{i} = zeros(size(rho{1})); 
       v2c{i} = zeros(size(rho{1}));
    end

    id1 = arho>1e-6;
    zeta = zeros(size(rho{1}));
    zeta(id1) = (rho{1}(id1)-rho{2}(id1))./arho(id1);
    
    id2 = abs(zeta)<=1;
    zeta(id2) = sign(zeta(id2)).*min(abs(zeta(id2)),1-1e-6);
    
    idxnz = arho>1e-6&sqrt(abs(grho2{1}))>1e-6&abs(zeta)<=1.0;
    % PBC
    [uc(idxnz),v1c{1}(idxnz),v1c{2}(idxnz),v2c{1}(idxnz)] = VGCCorrelation_pbc_spin...
    (arho(idxnz),zeta(idxnz),grho2{1}(idxnz),1);
    v2c{2} = v2c{1};
    v2c_ud = v2c{1};    
end


