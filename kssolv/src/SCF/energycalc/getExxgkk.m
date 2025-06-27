function [options,idxnz] = getExxgkk(mol,options)
% GETEXXGKK Compute eigenvalues of the exact exchange operator with k points.
%    options = getExxgkk(mol,options) calculates the matrix elements by
%    V_{kk'GG'} = <kG|V|k'G'> where kk' is k-point index, GG' is G-grid
%    index, V is Coulomb kernel.
%
%    The calculated results are stored in options.exxgkk.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

options.exxcut = mol.ecut2;
grid = Ggrid(mol, options.exxcut);
idxnz = grid.idxnz;
n1 = mol.n1;n2 = mol.n2;n3 = mol.n3;
n123 = n1*n2*n3;

if isa(mol,'Crystal')
    if (~options.aceconv && ~options.store_tensors) || ...
            strcmp(options.exxmethod,'normal') || ...
            strcmp(options.exxmethod,'ace')
        nk_psi = mol.nkpts;
        kpts_psi = mol.kpts;
        if isfield(options,'scfX') && ~isempty(options.scfX)
            nk_phi = options.scfX.nkpts;
            kpts_phi = mol.scfkpts;
        else
            nk_phi = nk_psi;
            kpts_phi = kpts_psi;
        end
        exxgkk=zeros(n123,nk_psi,nk_phi);
        for i=1:nk_psi
            for j=1:nk_phi
                %In QE, k+q=kq, i is index of k, j is index of kq.
                %The input argument of the formulation of exxgkk is k-kq=-q.
                %See function g2_convolution of exx_base.f90 in QE.
                q=kpts_psi(i,:)-kpts_phi(j,:);
                [exxgkk(idxnz,i,j),options]=getExxgkk_q(mol,options,q);
            end
        end
    else
        [kpts2,exxidxnz]=kgrid(mol);
        nkpts2=size(kpts2,1);
        exxgkk=zeros(n123,nkpts2);
        for ik=1:nkpts2
            [exxgkk(idxnz,ik),options]=getExxgkk_q(mol,options,kpts2(ik,:));
        end
        options.exxidxnz=exxidxnz;
    end
else
    exxgkk = zeros(n123,1);
    exxgkk(idxnz) = getExxgkk_q(mol,options);
end
options.exxgkk=exxgkk;
end

function exxgkk = getExxgkk_q(mol,options,q)
% GETEXXGKK_Q Compute eigenvalues of the exact exchange operator.
%    options = getExxgkk(mol,options) calculates the matrix elements by
%    V_{GG'} = <G|V|G'> where GG' is G-grid index, V is Coulomb kernel.
%
%    The calculated results are stored in options.exxgkk.
%    This function is equivalent to the function "g2_convolution"
%    in exx_base.f90 of QE.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

% Gygi-Baldereschi regularization is used to treat Coumlomb singularity
Mu = mol.omega; % Screening Parameter
epsDiv = 1e-8;  % Divergence threshold
grid = Ggrid(mol, options.exxcut);

if isprop(mol,'nqs') && ~isempty(mol.nqs)
    nq=mol.nqs;
elseif isa(mol,'Crystal')
    nq=mol.nkxyz;
else
    nq=[1,1,1];
end
if options.x_gamma_extrapolation
    grid_factor=8/7;
else
    grid_factor=1;
end

if (nargin < 3)
    gkxyz=[grid.gkx,grid.gky,grid.gkz];
    gkk = grid.gkk;
else
    gkxyz=[grid.gkx,grid.gky,grid.gkz]+q;
    gkk = sum(gkxyz.^2,2);
end
grid_factor_track = ones(size(gkk));
if options.x_gamma_extrapolation
    x = gkxyz*(mol.supercell.'.*nq/pi/4);
    on_double_grid=sum(abs(x-round(x)),2)<3*eps;
    grid_factor_track(on_double_grid)=0;
    grid_factor_track(~on_double_grid)=grid_factor;
end

if Mu == 0 || mol.exx_sr == mol.exx_lr
    exxgkk = mol.exx_sr * 4 * pi ./ gkk;
else
    exxgkk = mol.exx_sr * 4 * pi ./ gkk + (mol.exx_lr - mol.exx_sr) * 4 * pi ./ gkk .* exp(-gkk / (4*Mu^2));
end
exxgkk = exxgkk .* grid_factor_track;

exxDiv = getexxDiv(mol,options);
idx = gkk < epsDiv;
exxgkk(idx) = -exxDiv;

    function exxDiv = getexxDiv(mol,options)
        % Calculate the divergent G=0 term.
        % Equivalent to the function "exx_divergence" in exx_base.f90 of QE.
        if isfield(options,'exxDiv')
            exxDiv=options.exxDiv;
        else
            % Gygi-Baldereschi regularization
            exxAlpha = 10 / options.exxcut*2;
            rcell=inv(mol.supercell)'*2*pi;
            exxDiv=0;
            for iq1=0:nq(1)-1
                for iq2=0:nq(2)-1
                    for iq3=0:nq(3)-1
                        qxyz=[iq1/nq(1),iq2/nq(2),iq3/nq(3)]*rcell;
                        gkxyz=[grid.gkx,grid.gky,grid.gkz]+qxyz;
                        if options.x_gamma_extrapolation
                            x = gkxyz*(mol.supercell.'.*nq/pi/4);
                            on_double_grid=sum(abs(x-round(x)),2)<3*eps;
                        else
                            on_double_grid = zeros(size(gkxyz,1),1);
                        end
                        gkk = sum(gkxyz.*gkxyz,2);
                        
                        if Mu == 0 || mol.exx_sr == mol.exx_lr
                            gkki = mol.exx_sr * exp(-exxAlpha * gkk) ./ gkk;
                        else
                            gkki = mol.exx_sr * exp(-exxAlpha * gkk) ./ gkk + (mol.exx_lr - mol.exx_sr) * exp(-exxAlpha * gkk) ./ gkk .* exp(-gkk / (4*Mu^2));
                        end
                        
                        idx = (gkk > epsDiv)&~on_double_grid;
                        exxDiv = exxDiv + sum(gkki(idx))*grid_factor;
                    end
                end
            end
            if ~options.x_gamma_extrapolation
                if Mu == 0 || mol.exx_sr == mol.exx_lr
                    exxDiv = exxDiv - mol.exx_sr * exxAlpha;
                else
                    exxDiv = exxDiv - mol.exx_lr * exxAlpha;
                end
            end
            exxDiv = 4 * pi * exxDiv;
            
            if Mu == 0 || mol.exx_sr == mol.exx_lr
                aa = mol.exx_sr / sqrt(pi * exxAlpha);
            else
                aa = mol.exx_sr / sqrt(pi * exxAlpha) + (mol.exx_lr - mol.exx_sr) / sqrt(pi * (exxAlpha + 1 / (4 * Mu^2)));
            end
            exxDiv = exxDiv - mol.vol * aa * prod(nq);
        end
    end
end