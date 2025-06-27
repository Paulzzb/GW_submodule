function [egcx,egcc,v1gcx,v2gcx,v1gcc,v2gcc,v2gcc_ud] = Vxc_gcxc(rho,grho,mol)
% VXC_GCXC Graident part of GGA functional for spin systems.
%    [egcx,egcc,v1gcx,v2gcx,v1gcc,v2gcc,v2gcc_ud] = VXC_GCXC(rho,grho,mol) calculate
%    the gradient-dependent part in GGA functional for spin-polarized
%    systems and spin-noncollinear systems.
%
%    PBE case: pbx + pbc
%    Hybrid case: pbx*alpha + pbxsr*beta + pbc*gamma

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

grho2 = cell(2,1);
for i = 1 : 2
	grho2{i} = sum(grho{i}.^2,4);
end

% PBX
[v1gcx,v2gcx,egcx] = VGCExchange_pbx_spin(rho,grho2);

if nargin > 2
    alpha = 1.0 - mol.exx_lr;
    beta = mol.exx_lr - mol.exx_sr;
    if beta ~= 0
        [v1gcxsr,v2gcxsr,egcxsr] = VGCExchange_pbxsr_spin(rho,grho2,mol.omega);
        egcx = egcx*alpha + egcxsr*beta;
        for i = 1:2
            v1gcx{i} = v1gcx{i}*alpha + v1gcxsr{i}*beta;
            v2gcx{i} = v2gcx{i}*alpha + v2gcxsr{i}*beta;
        end
    else
        egcx = egcx*alpha;
        for i = 1:2
            v1gcx{i} = v1gcx{i}*alpha;
            v2gcx{i} = v2gcx{i}*alpha;
        end
    end
end

for i = 1 : 2
    v2gcx{i} = 2*v2gcx{i};  
end

% PBC
arho = abs(rho{1} + rho{2});
grho_tmp = grho{1} + grho{2};
grho2{1} = sum(grho_tmp.^2,4);

egcc = zeros(size(rho{1}));
v1gcc   = cell(2,1);
v2gcc   = cell(2,1);
for i = 1:2
	v1gcc{i} = zeros(size(rho{1})); 
    v2gcc{i} = zeros(size(rho{1}));
end

id1 = arho > 1e-6;
zeta = zeros(size(rho{1}));
zeta(id1) = (rho{1}(id1)-rho{2}(id1))./arho(id1);
    
id2 = abs(zeta) <= 1;
zeta(id2) = sign(zeta(id2)).*min(abs(zeta(id2)),1-1e-6);
    
idxnz = id1 & grho2{1} > 1e-12 & id2;

[v1gccnz,v2gcc{1}(idxnz),egcc(idxnz)] = VGCCorrelation_pbc_spin(arho(idxnz),zeta(idxnz),grho2{1}(idxnz));
v1gcc{1}(idxnz) = v1gccnz{1};
v1gcc{2}(idxnz) = v1gccnz{2};

if nargin > 2
    gamma = mol.dft_c;
    egcc = egcc*gamma;
    v1gcc{1} = v1gcc{1}*gamma;
    v1gcc{2} = v1gcc{2}*gamma;
    v2gcc{1} = v2gcc{1}*gamma;
end

v2gcc{2} = v2gcc{1};
v2gcc_ud = v2gcc{1};

end
