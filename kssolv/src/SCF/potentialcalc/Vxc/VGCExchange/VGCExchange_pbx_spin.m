function [v1gcx,v2gcx,egcx] = VGCExchange_pbx_spin(rho,grho2)
% VGCEXCHANGE_PBX_SPIN PBE exchange gradient correction for spin-polarized systems.
%    [v1gcx,v2gcx,egcx] = VGCEXCHANGE_PBX_SPIN(rho,grho2) returns the PBE
%    exchange gradient correction of the rho and the gradient square of rho.
%
%   See also exRef.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

% Filter the points where the density or gradient of density is very 
% small to avoid the overflow of machine number
idxnz_up = rho{1}>1e-10 & grho2{1}>1e-20;
idxnz_dw = rho{2}>1e-10 & grho2{2}>1e-20;
for i = 1 : 2
    rho{i} = 2*rho{i};
    grho2{i} = 4*grho2{i};
end
v1gcx = cell(2,1);
v2gcx = cell(2,1);
for i = 1:2
    v1gcx{i} = zeros(size(rho{1}));
    v2gcx{i} = zeros(size(rho{1}));
end
egcx_up = zeros(size(rho{1}));
egcx_dw = zeros(size(rho{1}));

[v1gcx{1}(idxnz_up),v2gcx{1}(idxnz_up),egcx_up(idxnz_up)] = ...
    VGCExchange_pbx(rho{1}(idxnz_up),grho2{1}(idxnz_up));
[v1gcx{2}(idxnz_dw),v2gcx{2}(idxnz_dw),egcx_dw(idxnz_dw)] = ...
    VGCExchange_pbx(rho{2}(idxnz_dw),grho2{2}(idxnz_dw));

egcx = (egcx_up + egcx_dw)/2;

end