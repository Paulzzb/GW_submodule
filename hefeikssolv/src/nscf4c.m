function [cry,BH,BX,info] = nscf4c(cry,options)
% NSCF4C Non Self Consistent Field iteration for crystal.
%    [cry,BX,info] = NSCF4C(cry,options) find the ground state minimum
%    total energy and the corresponding wave functions for cry. The initial
%    density must be provided in options.
%
%   See also scf4m, scf4c, dcm, trdcm.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if nargin < 2 || isempty(options.rho0)
    error('The initial density is not provided');
end

% Set timer
tstart  = cputime;

% Initialize input variables
rho        = options.rho0;
BX0        = options.X0;
nkpts      = cry.nkpts;
temperature= cry.temperature;
Tbeta      = temperature*8.6173324e-5/13.6;
maxPhiIter = options.maxphiiter;

% Initialize Hamiltonian, Wavefun, and Preconditioners
[cry,BH,BX,BHprec,noccs] = iterinit(cry,rho,BX0);

nBXcols    = ncols(BX);

% Initialize output variables
ev         = zeros(sumel(nBXcols),1);

if isempty(options.scfX)
    maxPhiIter = 1;
    fprintf('Beging Non-SCF calculation for %s...\n',cry.name);
else
    fprintf('Beging Fake-SCF calculation for %s...\n',cry.name);
end

exxmethod = options.exxmethod;
for iterphi = 1:maxPhiIter
    fprintf('Phi iter %3d:\n', iterphi);
    if iterphi == 1
        fprintf('GACE\n')
        options.exxmethod = 'nscf';
    elseif iterphi == 2
        options.exxmethod = exxmethod;
        options.aceconv = false;
        %options.exxmethod = 'normal';
        options = getExxgkk(cry,options);
    end
    [BH.vexx, cry]  = getVexx(BX, cry, options);
    %a=diag(options.scfX{1}'*BH.vexx(options.scfX{1}));a(1:6)
    idx = 0;
    for ik = 1:nkpts
        idx = idx(end) + (1:nBXcols(ik));
        [X, ev(idx)] = updateX(cry, BH{ik}, BX{ik}, BHprec{ik}, options);
        BX{ik} = X;
    end
    %a=diag(BX{1}'*BH.vexx(BX{1}));a(1:6)
    if iterphi>=2
        ev_err = max(abs(last_ev(:)-ev(:)));
        fprintf('Fake-SCF max delta eigenvalues: %20.13e\n', ev_err);
        if ev_err<5e-4
            break
        end
    end
%     if iterphi>=3
%         zeta_err=norm(last_zeta(:)-cry.zeta(:));
%         fprintf('mean zeta: %20.13e\n',mean(abs(cry.zeta(:))));
%         fprintf('zeta error: %20.13e\n', zeta_err);
%     end
    last_ev=ev;
%     last_zeta=cry.zeta;
end

if iterphi>=2 && ev_err<5e-4
    fprintf('Fake-SCF convergence is reached.\n');
elseif iterphi>=2
    fprintf('######################################\n');
    fprintf('Warning: Fake-SCF not converge!\n');
end
%[occs,~] = getocc(ev,noccs,Tbeta);

%BX = assignoccs(BX,occs);

info.Eigvals = ev;

timetot = cputime - tstart;
fprintf('Total time used = %20.3e\n', timetot);

end
