function fxc = getfxc(mol,rho)
% GETVXC exchange correlation.
%    [vxc,uxc,rho] = GETVXC(mol,rho) returns the exchange correlation of
%    the rho. The type of the exchange correlation is determined by the
%    pseudopotential.
%
%   See also exRef.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

switch upper(mol.ppvar.funct)
    case 'PBE'
        fprintf('PBE functional NYI');
        pause;
        %fxc = getfxcPBE(mol,rho);
    case 'SLA-PW-PBX-PBC'
        fprintf('SLA-PW-PBX-PBC functional NYI');
        pause;
        %fxc = getfxcPBE(mol,rho);
    case 'PZ'
        fxc = getfxcPZ(rho);
    otherwise
        fxc = getfxcPZ(rho);
end

end

