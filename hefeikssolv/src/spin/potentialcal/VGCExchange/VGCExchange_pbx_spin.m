function [ugcx,v1gcx,v2gcx] = VGCExchange_pbx_spin(rho,grho2,ishybrid)
% Gradient corrections for exchange  

% Filter the points where the density or gradient of density is very 
% small to avoid the overflow of  machine number
idxnz_up = rho{1}>1e-10 & sqrt(abs(grho2{1}))>1e-10;
idxnz_dw = rho{2}>1e-10 & sqrt(abs(grho2{2}))>1e-10;    
for i = 1 : 2
    rho{i} = 2*rho{i};
    grho2{i} = 4*grho2{i};
end
v1gcx = cell(2,1);
v2gcx = cell(2,1);
v1gcxsr = cell(2,1);
v2gcxsr = cell(2,1);

for i = 1:2
    v1gcx{i} = zeros(size(rho{1}));
    v2gcx{i} = zeros(size(rho{1}));
    v1gcxsr{i} = zeros(size(rho{1}));
    v2gcxsr{i} = zeros(size(rho{1}));
end
ugcx_up = zeros(size(rho{1}));
ugcx_dw = zeros(size(rho{1}));
exsr_up = zeros(size(rho{1}));
exsr_dw = zeros(size(rho{1}));
% PBX density functional
[v1gcx{1}(idxnz_up),v2gcx{1}(idxnz_up),ugcx_up(idxnz_up)] = ...
VGCExchange_pbx(rho{1}(idxnz_up),grho2{1}(idxnz_up));
[v1gcx{2}(idxnz_dw),v2gcx{2}(idxnz_dw),ugcx_dw(idxnz_dw)] = ...
VGCExchange_pbx(rho{2}(idxnz_dw),grho2{2}(idxnz_dw));

ugcx = (ugcx_up + ugcx_dw)/2;
for i = 1 : 2
    v2gcx{i} = 2*v2gcx{i};  
end

if ishybrid==1
    % For HSE06 case, subtract  1/4 of the short range part of exchange energe
    % Instead , add 1/4 of the hartree-fock exchange energy
    [exsr_up(idxnz_up),v1gcxsr{1}(idxnz_up),v2gcxsr{1}(idxnz_up)] = ...
        VGCExchange_pbxsr(rho{1}(idxnz_up),grho2{1}(idxnz_up));
    [exsr_dw(idxnz_dw),v1gcxsr{2}(idxnz_dw),v2gcxsr{2}(idxnz_dw)] = ...
        VGCExchange_pbxsr(rho{2}(idxnz_dw),grho2{2}(idxnz_dw));
    ugcx = ugcx - 0.25*(exsr_up + exsr_dw)/2;
    for i = 1:2
        v1gcx{i} = v1gcx{i} - 0.25*v1gcxsr{i};
        v2gcx{i} = v2gcx{i} - 0.25*v2gcxsr{i}*2;
    end
end
end

