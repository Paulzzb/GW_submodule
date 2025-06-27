function [v1gcxsr,v2gcxsr,egcxsr] = VGCExchange_pbxsr_spin(rho,grho2,omega)
% VGCEXCHANGE_PBXER_SPIN PBE short-range exchange gradient correction for spin-polarized systems.
%    [v1gcxsr,v2gcxsr,egcxsr] = VGCEXCHANGE_PBXSR_SPIN(rho,grho2,omega) returns the PBE short-range
%    exchange gradient correction of the rho and the gradient square of rho, which already 
%    contains the Slater short-range exchange.
%
%   See also exRef.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

% Filter the points where the density or gradient of density is very 
% small to avoid the overflow of machine number
rho{1} = abs(rho{1});
rho{2} = abs(rho{2});
idxnz_up = rho{1}>1e-10 & grho2{1}>1e-20;
idxnz_dw = rho{2}>1e-10 & grho2{2}>1e-20;
for i = 1 : 2
    rho{i} = 2*rho{i};
    grho2{i} = 4*grho2{i};
end
v1gcxsr = cell(2,1);
v2gcxsr = cell(2,1);
for i = 1:2
    v1gcxsr{i} = zeros(size(rho{1}));
    v2gcxsr{i} = zeros(size(rho{1}));
end
egcxsr_up = zeros(size(rho{1}));
egcxsr_dw = zeros(size(rho{1}));

[v1gcxsr{1}(idxnz_up),v2gcxsr{1}(idxnz_up),egcxsr_up(idxnz_up)] = ...
    VGCExchange_pbxsr(rho{1}(idxnz_up),grho2{1}(idxnz_up),omega);
[v1gcxsr{2}(idxnz_dw),v2gcxsr{2}(idxnz_dw),egcxsr_dw(idxnz_dw)] = ...
    VGCExchange_pbxsr(rho{2}(idxnz_dw),grho2{2}(idxnz_dw),omega);

egcxsr = (egcxsr_up + egcxsr_dw)/2;

end