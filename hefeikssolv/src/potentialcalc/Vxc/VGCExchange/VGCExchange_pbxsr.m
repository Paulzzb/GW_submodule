function [ugcxsr,v1gcxsr,v2gcxsr] = VGCExchange_pbxsr(rho,grho2,omega)
% VGCEXCHANGE_PBXSR PBE short-range exchange gradient correction.
%    [v1gcxsr,v2gcxsr,ugcxsr] = VGCEXCHANGE_PBXSR(rho,grho2) returns the short
%    range PBE exchange gradient correction of the rho and the gradient square 
%    of rho.
%
%   See also exRef.

if nargin < 3 
    omega = 0.106;
end

[ugcxsr,v1gcxsr,v2gcxsr] = pbexsrF(rho(:),grho2(:),omega);
%[ugcxsr,v1gcxsr,v2gcxsr] = pbexsr_wrap(rho(:),grho2(:),omega);
if length(rho) > 1
    v1gcxsr = reshape(v1gcxsr, size(rho));
    v2gcxsr = reshape(v2gcxsr, size(rho));
    ugcxsr  = reshape(ugcxsr,  size(rho));
end

end
