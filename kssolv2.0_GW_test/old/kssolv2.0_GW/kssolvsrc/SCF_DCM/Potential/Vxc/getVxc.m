function [vxc,uxc,rho,uxcsr] = getVxc(mol,rho)
% GETVXC exchange correlation.
%    [vxc,uxc,rho] = GETVXC(mol,rho) returns the exchange correlation of
%    the rho. The type of the exchange correlation is determined by the
%    pseudopotential.
%
%   See also exRef.

uxcsr = zeros(size(rho));
switch upper(mol.ppvar.funct)
    case 'HSE06'
        [vxc,uxc,rho,uxcsr] = VxcHSE06(mol,rho);
    case 'PBE'
        [vxc,uxc,rho] = VxcPBE(mol,rho);
    case 'SLA-PW-PBX-PBC'
        [vxc,uxc,rho] = VxcPBE(mol,rho);
    case 'PZ'
        [vxc,uxc,rho] = VxcPZ(rho);
    otherwise
        [vxc,uxc,rho] = VxcPZ(rho);
end

vxc = bandlimV(mol,vxc);

end
