function rho = corrcharge(mol,rho)
% CORRCHARGE correct the charge density.
%    rho = CORRCHARGE(mol,rho) correct the charge density, i.e., correct
%    negative density, renormalise density.
%
%   See also getcharge.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

n1    = mol.n1;
n2    = mol.n2;
n3    = mol.n3;
nspin = mol.nspin;
n123  = n1*n2*n3;
nel   = mol.nel;

if nspin == 1
    arho = rho;
elseif nspin == 2
    arho = rho{1}+rho{2};
elseif nspin == 4
    arho = rho{1};
end
charge = sumel(arho)*mol.vol/n123;

if abs(charge-nel)/nel > 1e-7
    global verbose;
    if verbose >=0
        warning(['Renormalize starting charge ' num2str(charge) ...
        	' to ' num2str(nel)]);
    end

    if charge > 1e-8
        if nspin == 1
            rho = rho*(nel/charge);
        elseif nspin == 2||nspin == 4
            for i = 1:nspin
                rho{i} = rho{i}*(nel/charge);
            end
        end
    else
        if nspin == 1
            rho = nel/mol.vol*ones(size(rho));
        elseif nspin == 2
            rho{1} = nel/2/mol.vol*ones(size(rho{1}));
            rho{2} = nel/2/mol.vol*ones(size(rho{2}));
        elseif nspin == 4
            rho{1} = nel/mol.vol*ones(size(rho));
        end
    end
else

end
