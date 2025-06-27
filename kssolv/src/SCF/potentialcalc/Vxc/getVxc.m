function [vxc,exc,rho] = getVxc(mol,rho)
% GETVXC computes the exchange correlation potential (vxc) and energy per unit volume (exc).
%    [vxc,exc,rho] = GETVXC(mol,rho) returns the exchange correlation of
%    the rho. The type of the exchange correlation is determined by the
%    pseudopotential.
%
%   See also exRef.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

nspin = mol.nspin;

switch upper(mol.funct)
    case {'PZ','PW','PBE','SLA-PW-PBX-PBC'}
        if mol.exx_sr ~= 0 || mol.exx_lr ~= 0 || mol.dft_c ~= 1
            error('Hybrid parameters are not available for pure DFT.')
        end
end

switch upper(mol.funct)
    case 'PZ'
        if nspin == 1
            [vxc,exc,rho] = VxcPZ(rho);
        elseif nspin == 2
            [vxc,exc,rho] = VxcPZ_lsda(rho);
        elseif nspin == 4
            if mol.domag
                [vxc,exc,rho] = VxcPZ_nc(rho);
            else
                [vxc,exc,~] = VxcPZ(rho{1});
            end
        end
    case 'PW'
        if nspin == 1
            [vxc,exc,rho] = VxcPW(rho);
        elseif nspin == 2
            [vxc,exc,rho] = VxcPW_lsda(rho);
        elseif nspin == 4
            if mol.domag
                [vxc,exc,rho] = VxcPW_nc(rho);
            else
                [vxc,exc,~] = VxcPW(rho{1});
            end
        end
    case {'PBE','SLA-PW-PBX-PBC'}
        if nspin == 1
            [vxc,exc,rho] = VxcPBE(mol,rho);
        elseif nspin == 2
            [vxc,exc,rho] = VxcPBE_lsda(mol,rho);
        elseif nspin == 4
            if mol.domag
                [vxc,exc,rho] = VxcPBE_nc(mol,rho);
            else
                [vxc,exc,~] = VxcPBE(mol,rho{1});
            end
        end
    case {'HSE06','HSE','PBE0','LC-WPBE','LC-PBE0','HF','PUREHF'}
        if nspin == 1
            [vxc,exc,rho] = VxcHSE06(mol,rho);
        elseif nspin == 2
            [vxc,exc,rho] = VxcHSE06_lsda(mol,rho);
        elseif nspin == 4
            if mol.domag
                [vxc,exc,rho] = VxcHSE06_nc(mol,rho);
            else
                [vxc,exc,~] = VxcHSE06(mol,rho{1});
            end
        end
    otherwise
        error('Unknown functional is used.')
end

end
