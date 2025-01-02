function [vxc,uxc,rho,uxcsr] = getVxc(mol,rho)
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
nspin = mol.nspin;
if nspin == 1
    uxcsr = zeros(size(rho));
elseif nspin == 2||nspin == 4
    uxcsr = zeros(size(rho{1}));
end
switch upper(mol.funct)
    case {'HSE06','HSE'}
	%fprintf('The type of the exchange correlation is HSE\n');        
        if nspin == 1
            [vxc,uxc,rho,uxcsr] = VxcHSE06(mol,rho);
        elseif nspin == 2
            [vxc,uxc,rho] = VxcHSE06_lsda(mol,rho);
        elseif nspin == 4
            if mol.domag
                [vxc,uxc,rho] = VxcHSE06_nc(mol,rho);
            else
                [vxc,uxc,~,uxcsr] = VxcHSE06(mol,rho{1});
            end
        end         
    case 'PBE'
	%fprintf('The type of the exchange correlation is PBE\n');
        if nspin == 1
            [vxc,uxc,rho] = VxcPBE(mol,rho);
        elseif nspin == 2
            [vxc,uxc,rho] = VxcPBE_lsda(mol,rho);
        elseif nspin == 4
            if mol.domag
                [vxc,uxc,rho] = VxcPBE_nc(mol,rho);
            else
                [vxc,uxc,~] = VxcPBE(mol,rho{1});
            end
        end
    case 'SLA-PW-PBX-PBC'
	%fprintf('The type of the exchange correlation is PBE\n');
        if nspin == 1
            [vxc,uxc,rho] = VxcPBE(mol,rho);
        elseif nspin == 2
            [vxc,uxc,rho] = VxcPBE_lsda(mol,rho);
        elseif nspin == 4
            if mol.domag
                [vxc,uxc,rho] = VxcPBE_nc(mol,rho);
            else
                [vxc,uxc,~] = VxcPBE(mol,rho{1});
            end
        end
    case 'PZ'
	%fprintf('The type of the exchange correlation is PZ\n');
         if nspin == 1
            [vxc,uxc,rho] = VxcPZ(rho);
        elseif nspin == 2
            [vxc,uxc,rho] = VxcPZ_lsda(rho);
        elseif nspin == 4
            if mol.domag
                [vxc,uxc,rho] = VxcPZ_nc(rho);
            else
                [vxc,uxc,~] = VxcPZ(rho{1});
            end
        end
    otherwise
	%fprintf('The type of the exchange correlation is PZ\n');
         if nspin == 1
            [vxc,uxc,rho] = VxcPZ(rho);
        elseif nspin == 2
            [vxc,uxc,rho] = VxcPZ_lsda(rho);
        elseif nspin == 4
            if mol.domag
                [vxc,uxc,rho] = VxcPZ_nc(rho);
            else
                [vxc,uxc,~] = VxcPZ(rho{1});
            end
        end
end

end
