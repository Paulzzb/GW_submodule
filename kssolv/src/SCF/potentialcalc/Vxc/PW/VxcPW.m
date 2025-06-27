function [vxc,exc,rho] = VxcPW(rho)
% VXCPW PW exchange correlation for spin-restricted systems.
%    [vxc,exc,rho] = VXCPW(rho) returns the PW exchange correlation of the
%    one-component rho for spin-restricted systems.
%
%   See also exRef.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

e2 = e2Def();

vxc = zeros(size(rho));
exc = zeros(size(rho));

idxxc = abs(rho) > 1e-10;
rhoxc = abs(rho(idxxc));

rs = ((3/(4*pi))./rhoxc).^(1/3);
[vx,ux] = VExchange_sla(rs);
[vc,uc] = VCorrelation_pw(rs);

vxc(idxxc) = e2*(vx+vc);
exc(idxxc) = e2*(ux+uc).*rhoxc;

end
