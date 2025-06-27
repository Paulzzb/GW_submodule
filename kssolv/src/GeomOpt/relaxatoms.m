function [mol, H, X, info, varargout] = relaxatoms(mol, ksopts, H, X, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% RELAXATOMS optimize the positions of atoms to minimize total energy
% method = 1: use MATLAB optimization toolbox (quasi-newton)
%             https://ww2.mathworks.cn/help/optim/ug/fminunc.html
% method = 2: nonlinear conjugate gradient implemented by Amartya S. Banerjee
%             https://www.amartyabanerjee.com/
% method = 3: nonlinear conjugate gradient implemented by Michael L. Overton
%             https://www.cs.nyu.edu/overton/software/nlcg/
% method = 4: BFGS package by Michael L. Overton [part of the
%             HANSO: Hybrid Algorithm for Non-Smooth Optimization (V2.2)]
%             https://www.cs.nyu.edu/overton/software/hanso/
% method = 5: Fast Inertial Relaxation Engine (FIRE) by Amartya S. Banerjee
%             DOI: 10.1103/PhysRevLett.97.170201
%
% Updated by Linhao Wang on 2025/3/23
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nVarargs = length(varargin);
natoms = sum(mol.natoms);

% Ensure the parameters for optimization
if nargin == 1
    ksopts = setksopt();
    ksopts.relaxmethod = 'fminunc';
    ksopts.relaxtol = 1e-4;
    ksopts.maxrelaxiter = 200;
end
switch lower(ksopts.relaxmethod)
    case 'fminunc'
        method = 1;
    case {'nlcg1','nlcgb'}
        method = 2;
    case {'nlcg2','nlcgo'}
        method = 3;
    case 'bfgs'
        method = 4;
    case 'fire'
        method = 5;
    otherwise
        if isnumeric(ksopts.relaxmethod) && ksopts.relaxmethod == fix(ksopts.relaxmethod) && ksopts.relaxmethod <= 5
            method = ksopts.relaxmethod;
        else
            error('Error: invalid optimization method = %s\nmethod must be fminunc, nlcg1, nlcg2, bfgs or fire\n', ...
                num2str(ksopts.relaxmethod));
        end
end
ksopts.force = 1;

% Get the initial imformation
global saved_H saved_X relaxHistory;
if exist('H','var')
    saved_H = H;
else
    clear global saved_H
end
if exist('X','var')
    saved_X = X;
else
    clear global saved_X
end
relaxHistory = {};
x0 = mol.xyzlist(:);

% Select optimization method
if method == 1
    % use MATLAB optimization toolbox (quasi-newton)
    optimopts = optimoptions('fminunc');
    optimopts.Algorithm = 'quasi-newton';
    optimopts.Diagnostics = 'on';
    if verLessThan('optim', '7.4')
        optimopts.FinDiffRelStep = 0.01;
        optimopts.GradObj = 'on';
        optimopts.MaxFunEvals = ksopts.maxrelaxiter + 1;
        optimopts.MaxIter = ksopts.maxrelaxiter -1;
        optimopts.TolFun = ksopts.relaxtol/sqrt(3*natoms);              % Force
        optimopts.TolX = 0;                                             % Displacement (disabled)
    else
        optimopts.FiniteDifferenceStepSize = 0.01;
        optimopts.SpecifyObjectiveGradient = true;
        optimopts.MaxFunctionEvaluations = ksopts.maxrelaxiter + 1;
        optimopts.MaxIterations = ksopts.maxrelaxiter -1;
        optimopts.OptimalityTolerance = ksopts.relaxtol/sqrt(3*natoms); % Force
        optimopts.StepTolerance = 0;                                    % Displacement (disabled)
    end
    
    [x,FVAL,exitFlag,output,GRAD,~] = fminunc(@(x) ksefg(x,mol,ksopts), x0, optimopts);
    
    if nargout > 4
        varargout{1} = FVAL;                            % converged energy value
        varargout{2} = -reshape(GRAD,[natoms 3]);       % converged force value
    end
    
    xyzlist_new = reshape(x, [natoms 3]);
    mol = set(mol,'xyzlist',xyzlist_new);
    if exitFlag == 1
        converge = true;
        % 1  Magnitude of gradient small enough.
        % 2  Change in X too small.                     % won't happen (since TolX/StepTolerance = 0)
        % 3  Change in objective function too small.    % won't happen (for quasi-newton)
    else
        converge = false;
        % 5  Cannot decrease function along search direction.
        % 0  Too many function evaluations or iterations.
        % -1  Stopped by output/plot function.
        % -3  Problem seems unbounded.
    end
    relaxMessage = output.message;
    
elseif method == 2
    % nonlinear conjugate gradient implemented by Amartya S. Banerjee
    [mol, FORCE, Etot] = NLCG(mol, ksopts);
    
    if nargout > 4
        varargout{1} = Etot;                            % converged energy value
        varargout{2} = FORCE;                           % converged force value
%       varargout{3} = EtotRec;                         % history of energy values
%       varargout{4} = FORCERec;                        % history of force values
    end
    
%   fprintf('********************************\n');
%   fprintf('Number of Ionic Iterations: %3d \n', nIonicIter);
%   fprintf('********************************\n');
    
    if norm(FORCE) < ksopts.relaxtol
        converge = true;
        relaxMessage = 'gradient norm below tolerance\n';
    else
        converge = false;
        relaxMessage = 'error\n';
    end
    
elseif method == 3
    % nonlinear conjugate gradient implemented by Michael L. Overton
    if nVarargs == 2
        % output order: [x, f, g, frec, grec, alpharec]
        [x, f, g, frec, grec, ~] = nlcg(varargin{1}, varargin{2});
    else
        pars.nvar = length(x0);
        pars.fgname = 'ksefg';
        pars.mol = mol;
        pars.ksopts = ksopts;
        
        options.x0 = x0;
        options.nstart = 1;
        options.maxit = ksopts.maxrelaxiter;
        options.normtol = ksopts.relaxtol;
        options.version = 'C';
        % options.version (used to obtain different choices of beta):
        % 'P' for Polak-Ribiere-Polyak (not recommended: fails on hard problems)
        % 'F' for Fletcher-Reeves (not recommended: often stagnates)
        % 'C' for Polak-Ribiere-Polyak Constrained by Fletcher-Reeves
        % (recommended, combines advantages of 'P' and 'F'; default)
        % 'S' for Hestenes-Stiefel (not recommended)
        % 'Y' for Dai-Yuan (allows weak Wolfe line search, see nlcg.m)
        % 'Z' for Hager-Zhang
        % '-' for Steepest Descent (for comparison)
        
        options.strongwolfe = 1;   	% strong Wolfe line search (default)
        options.wolfe1 = 1e-4;  	% 0 < c1 < c2 < 1/2
        options.wolfe2 = 0.49;   	% c2 < 1/2 for NLCG
        options.prtlevel = 1;      	% minimal printing (default)
        
        % output order: [x, f, g, frec, grec, alpharec]
        [x, f, g, frec, grec, ~] = nlcg(pars, options);
    end
    nIonicIter = length(frec{1});
    
    if nargout > 4
        varargout{1} = f;                               % converged energy value
        varargout{2} = -reshape(g,[natoms 3]);          % converged force value
        varargout{3} = (frec{1})';                      % history of energy values
        
        % Get the history of force on each atom at every iteration:
        fAtoms = cell(nIonicIter,1);
        for n = 1:nIonicIter
            fAtoms{n} = -reshape(grec{1}(:,n),[natoms 3]);
        end
        varargout{4} = fAtoms;                          % history of force values
    end
    
    fprintf('********************************\n');
    fprintf('Number of Ionic Iterations: %3d \n', nIonicIter);
    fprintf('********************************\n');
    
    xyzlist_new = reshape(x,[natoms 3]);
    mol = set(mol,'xyzlist',xyzlist_new);
    if norm(g) < options.normtol
        converge = true;
        relaxMessage = 'gradient norm below tolerance\n';
    else
        converge = false;
        relaxMessage = 'error\n';
    end
    
elseif method == 4
    % BFGS package by Michael Overton [part of the
    % HANSO: Hybrid Algorithm for Non-Smooth Optimization (V2.2)]
    if nVarargs == 2
        % output order: [x, f, d, H, iter, info, X, G, w, fevalrec, xrec, drec, Hrec]
        [x, f, d, ~, nIonicIter, exitFlag, ~, ~, ~, frec, ~, drec, ~] = bfgs(varargin{1}, varargin{2});
    else
        pars.nvar = length(x0);
        pars.fgname = 'ksefg';
        pars.mol = mol;
        pars.ksopts = ksopts;
        
        options.x0 = x0;
        options.nstart = 1;
        options.maxit = ksopts.maxrelaxiter;
        options.nvec = 0;           % 0 => full BFGS, else => limited-memory BFGS
        options.H0 = sparse(eye(pars.nvar));
        options.scale = 1;
        options.normtol = ksopts.relaxtol;
        
        options.strongwolfe = 0;    % weak Wolfe line search (default)
        options.wolfe1 = 1e-4;      % 0 < c1 < c2 < 1
        options.wolfe2 = 0.5;       % c2 = 1/2 is also default
        options.ngrad = 1;          % saves the final gradient
        options.prtlevel = 1;       % minimal printing (default)
        
        % output order: [x, f, d, H, iter, info, X, G, w, fevalrec, xrec, drec, Hrec]
        [x, f, d, ~, nIonicIter, exitFlag, ~, ~, ~, frec, ~, drec, ~] = bfgs(pars, options);
    end
    
    if nargout > 4
        varargout{1} = f;                               % converged energy value
        varargout{2} = -reshape(d,[natoms 3]);          % converged force value
        varargout{3} = frec';                           % history of energy values
        
        % Get the history of force on each atom at every iteration:
        fAtoms = cell(nIonicIter,1);
        for n = 1:nIonicIter
            fAtoms{n} = -reshape(drec{n}(:,1),[natoms 3]);
        end
        varargout{4} = fAtoms;                          % history of force values
    end
    
    fprintf('********************************\n');
    fprintf('Number of Ionic Iterations: %3d \n', nIonicIter);
    fprintf('********************************\n');
    
    xyzlist_new = reshape(x,[natoms 3]);
    mol = set(mol,'xyzlist',xyzlist_new);
    if exitFlag == 0
        converge = true;
        relaxMessage = 'gradient norm below tolerance\n';
    else
        converge = false;
        relaxMessageLibrary = {'max number of iterations reached',...
            'energy reached target value',...           % won't happen
            'norm(position) exceeded limit',...         % won't happen
            'cpu time exceeded limit',...               % won't happen
            'energy or force is inf or nan at initial point',...
            'direction not a descent direction (because of rounding)',...
            'line search bracketed minimizer but Wolfe conditions not satisfied',...
            'line search did not bracket minimizer: energy may be unbounded below'};
        relaxMessage = ['error type ',num2str(exitFlag),': ',relaxMessageLibrary{exitFlag},'\n'];
    end
    
elseif method == 5
    % Fast Inertial Relaxation Engine (FIRE) by Amartya S. Banerjee
    if nVarargs == 3
        [mol, POS, FORCE, Etot, SCFINFO, FIREINFO, exitFlag, MinIdx] = quenchbyfire(mol, ksopts, ...
            varargin{1}, varargin{2}, varargin{3});
    else
        [mol, POS, FORCE, Etot, SCFINFO, FIREINFO, exitFlag, MinIdx] = quenchbyfire(mol, ksopts);
    end
    
    relaxHistory(:,end+2:end+8) = [{'velocity', 'dt', 'alpha', 'P', 'N_iter', 'N_cut', 'N(P>0)'}; FIREINFO; cell(1,7)];
    
    if nargout > 4
        varargout{1} = Etot(MinIdx);                    % converged energy value
        varargout{2} = FORCE{MinIdx};                   % converged force value
        varargout{3} = Etot;                            % history of energy values
        varargout{4} = FORCE;                           % history of force values
        
        varargout{5} = POS;                             % history of position values
        varargout{6} = SCFINFO;                         % history of SCF information
        varargout{7} = FIREINFO;                        % history of FIRE information
    end
    
    fprintf('******************************* \n');
    fprintf('Number of Ionic Iterations: %3d \n', max(0,length(Etot)-2));
    fprintf('******************************* \n');
    
    if exitFlag == 1
        converge = true;
        relaxMessage = 'gradient norm below tolerance\n';
    else
        converge = false;
        relaxMessage = 'error: max number of iterations reached\n';
    end
end

% Run SCF again to pass out H and X
[mol,H,X,lastIonicStepInfo] = scf(mol,ksopts);	% In some cases, the last step information in relaxHistory may be
                                            	% differemt from the true final structure.
info.converge = converge;
info.relaxHistory = relaxHistory;
info.relaxMessage = relaxMessage;
info.lastIonicStepInfo = lastIonicStepInfo;
info.lastIonicStepInfo.NormForce = norm(mol.xyzforce(:));

fprintf('norm(g) = %11.3e\n', info.lastIonicStepInfo.NormForce);
fprintf('\n');
fprintf('****** End of RelaxAtoms ****** \n');
fprintf(relaxMessage);
if converge
    fprintf('**** RelaxAtoms Converged! **** \n');
else
    fprintf('***** RelaxAtoms Failed!! ***** \n');
end

end
