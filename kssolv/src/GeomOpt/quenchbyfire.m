function [mol, POS, FORCE, Etot, SCFINFO, FIREINFO, exitFlag, MinIdx] = quenchbyfire(mol, ksopts, mass, dtStart, FIREpars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subhajit Banerjee
% June 2017
% SSG@CRD, LBNL
%
% Updated by Linhao Wang on 2025/3/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 2
    fprintf('FIRE executed with default parameter values, time step, and mass = 4 \n');
    % For appropriate mass value please review the following reference
    % http://nanosurf.fzu.cz/wiki/doku.php?id=fire_minimization
    % Also the mass is same for all atoms in FIRE
    mass = 4;
    
    dtStart = 0.05 * 41.3413733349185; % 1 Femtoseconds [fs] = 41.3413733349185 atomic time unit
    
    % Source:
    % https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.97.170201
    %
    % Qualitatively they should satisfy:
    % nMin larger than 1 (at least a few smooth steps after freezing);
    % fInc larger than but near to 1 (avoid too fast acceleration);
    % fDec smaller than 1 but much larger than 0 (avoid too heavy slowing down);
    % alphaStart larger than, but near to 0 (avoid too heavy damping);
    % fAlpha smaller than, but near to 1 (mixing is efficient some time after restart).
    FIREpars.nMin       = 5;
    FIREpars.fDec       = 0.5;
    FIREpars.fInc       = 1.1;
    FIREpars.alphaStart = 0.1;
    FIREpars.fAlpha     = 0.99;
    FIREpars.dtMax      = 10;
end

MAXITER         = ksopts.maxrelaxiter;
TOL             = ksopts.relaxtol;
FIREpars.dtMax  = dtStart * FIREpars.dtMax;

% Initialization
Etot = []; POS  = {}; FORCE = {}; SCFINFO = {}; FIREINFO = {};
x = mol.xyzlist(:);                                  	% initial position
[~, ~, mol, info] = ksefg(x, mol, ksopts);              % initial force (mol.xyzforce(:))
Etot(1)       = info.Etot;
POS{1}        = mol.xyzlist;
FORCE{1}      = mol.xyzforce;
SCFINFO{1}    = info;
checkConv = norm(mol.xyzforce(:));                      % check for convergence
if checkConv <= TOL
    exitFlag = true;
    MinIdx = 1;
    fprintf('\n');
    fprintf('The input structure has already converged.\n');
    fprintf('\n');
    return
else
    exitFlag = false;
end

v = zeros(size(x));                                 	% initial velocity
dt = dtStart;                                        	% initial time step
verletParams = struct('v0',v, 'dt',dt, 'mass',mass);  	% initial parameters for velocity Verlet
alpha = FIREpars.alphaStart;                        	% initial mixing coefficient for FIRE
FIREINFO(1,:) = {reshape(v, [sum(mol.natoms) 3]), dt, alpha, 0, 0, 0, 0};

% First velocity Verlet update with initial velocity = 0
fprintf('\n');
fprintf('First velocity Verlet update with initial velocity = 0.\n');
fprintf('\n');
[v, mol, Etot, POS, FORCE, SCFINFO] = velverletfire(mol, verletParams, ksopts, Etot, POS, FORCE, SCFINFO);

% Iteration
niter = 0; ncut = 0;
while niter <= MAXITER
    
    niter = niter + 1;
    
    checkConv = norm(mol.xyzforce(:));                  % check for convergence
    if checkConv < TOL
        exitFlag = true;
        break
    elseif niter == MAXITER + 1
        break
    else
        % Update v, dt, alpha using FIRE
        fprintf('Updating v, dt, and alpha using FIRE for iteration: %-5d\n\n', niter);
        [v, dt, alpha, ncut, FIREINFO] = fireupdater(niter, mol, v, dt, alpha, ncut, FIREpars, FIREINFO);
        verletParams = struct('v0',v, 'dt',dt, 'mass',mass);
        
        % Update x, F, v using velocity Verlet
        fprintf('Starting velocity Verlet update from FIRE output for iteration: %-5d\n\n', niter);
        [v, mol, Etot, POS, FORCE, SCFINFO] = velverletfire(mol, verletParams, ksopts, Etot, POS, FORCE, SCFINFO);
        fprintf('Completed velocity Verlet update for iteration: %-5d\n\n', niter);
    end
end

% Output
if exitFlag
    
    MinEtot = Etot(end);
    MinIdx = length(Etot);
    POSValAtMin   = [POS{end}];
    FORCEValAtMin = [FORCE{end}];
    
    fprintf('FIRE converged in %-5d iterations \n', niter-1);
    fprintf('The local minumum of energy is %-16.10f \n', MinEtot);
    fprintf('The corresponding configuration is:\n\n');
    disp(POSValAtMin);
    fprintf('The corresponding force is:\n\n');
    disp(FORCEValAtMin);
    
else
    
    [MinEtot, MinIdx] = min(Etot);
    POSValAtMin   = [POS{MinIdx}];
    mol = set(mol,'xyzlist',POSValAtMin);
    FORCEValAtMin = [FORCE{MinIdx}];
    
    fprintf('FIRE did not converge!\n');
    fprintf('Reached max number %-5d of iterations!\n', niter-1);
    fprintf('The infimum of energy is %-16.10f \n', MinEtot);
    fprintf('The infimum is achieved at the iteration # %-5d \n', MinIdx-2)
    fprintf('The corresponding configuration is:\n\n');
    disp(POSValAtMin);
    fprintf('The corresponding force is:\n\n');
    disp(FORCEValAtMin);
end

end

function [v, mol, Etot, POS, FORCE, SCFINFO] = velverletfire(mol, params, ksopts, Etot, POS, FORCE, SCFINFO)

x0 = mol.xyzlist(:);
F0 = mol.xyzforce(:);

v0 = params.v0;
dt = params.dt;
mass = params.mass;

% Position update in Velocity Verlet
x = x0 + dt*v0 + 0.5*dt*dt/mass*F0;

% Force update in SCF cycle
[~, ~, mol, info] = ksefg(x, mol, ksopts);
F = mol.xyzforce(:);

% Velocity update in Velocity Verlet
v = v0 + 0.5*dt/mass*(F0 + F);

Etot(end+1,1)    = info.Etot;
POS{end+1,1}     = mol.xyzlist;
FORCE{end+1,1}   = mol.xyzforce;
SCFINFO{end+1,1} = info;

if info.converge
    fprintf('\n');
    fprintf('SCF Converged in %d iterations.\n', length(info.Etotvec));
    fprintf('\n');
else
    fprintf('\n');
    fprintf('Warning: SCF did not converge.\n');
    fprintf('\n');
end

end

function [v, dt, alpha, ncut, FIREINFO] = fireupdater(niter, mol, v, dt, alpha, ncut, FIREpars, FIREINFO)

fDec            = FIREpars.fDec;
fInc            = FIREpars.fInc;
nMin            = FIREpars.nMin;
alphaStart      = FIREpars.alphaStart;
fAlpha          = FIREpars.fAlpha;
dtMax           = FIREpars.dtMax;

F = mol.xyzforce(:);
P = dot(F, v);      % Power
hatF = F/norm(F);   % a unit vector

% FIRE velocity update formula
if P > 0
    v = (1 - alpha)*v + alpha*norm(v)*hatF;
    if (niter - ncut) > nMin
        dt = min(dt*fInc, dtMax);   % increase dt
        alpha = alpha*fAlpha;       % decrease alpha
    end
else
    ncut = niter;                   % ncut <-- niter # (gets updated everytime P <= 0)
    v = zeros(size(F));             % reset velosity to 0
    dt = dt*fDec;                   % decrease dt
    alpha = alphaStart;             % reset alpha to alpha_start
end

FIREINFO(end+1,:) = {reshape(v, [sum(mol.natoms) 3]), dt, alpha, P, niter, ncut, niter - ncut};

end