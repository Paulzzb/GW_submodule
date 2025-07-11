function [mol, H, X, info, varargout] = relaxatoms(mol, ksopts, H, X, varargin)

% RELAXATOMS optimize the positions of atoms to minimize total energy
% method = 1: use MATLAB optimization toolbox
% method = 2: nonlinear conjugate gradient implemented by Overton
% method = 3: BFGS package by Overton [part of the 
%             HANSO: Hybrid Algorithm for Non-Smooth Optimization (V2.2)
%             http://www.cs.nyu.edu/overton/software/hanso/]

nVarargs = length(varargin);
natoms = sum(mol.natoms);


% Readjust the positions to fit within the supercell
% NOTE: Assuming cuboidal supercell
xxlim = mol.supercell(1,1);
yylim = mol.supercell(2,2);
zzlim = mol.supercell(3,3);

xtemp = mol.xyzlist(:,1) - (sum(mol.xyzlist(:,1))/natoms - xxlim/2.0); 
ytemp = mol.xyzlist(:,2) - (sum(mol.xyzlist(:,2))/natoms - yylim/2.0); 
ztemp = mol.xyzlist(:,3) - (sum(mol.xyzlist(:,3))/natoms - zzlim/2.0); 

xyztemp = [xtemp ytemp ztemp];

mol = set(mol,'xyzlist',xyztemp);
    
xyzforces = mol.xyzforce;


if (nargin == 1)
    ksopts = setksopt(ksopts,'scftol', 1e-6, 'cgtol',1e-8, 'maxscfiter', 100, 'maxcgiter',100);
    ksopts.relaxmethod = 'fminunc';
    ksopts.relaxtol = 1e-4;
end

% get the initial force (gradient)
if ( isempty(xyzforces) )
    % run SCF to obtain the forces on the atoms in MOL
    [mol,H,X,~] = scf4m(mol);
    ksopts = setksopt(ksopts,'rho0', H.rho, 'verbose','off');
else
    if (nargin > 2)
        ksopts = setksopt(ksopts, 'rho0',H.rho, 'verbose','off');           % USE with caution. Seems to worsen SCF convergence 
%         ksopts = setksopt('verbose','off');
    end
    if (nargin > 3)
%         ksopts = setksopt(ksopts, 'X0',X, 'verbose','off');       % USE with caution. Seems to worsen SCF convergence 
        ksopts = setksopt(ksopts, 'verbose','off');
    end
end
global saved_mol saved_H saved_X;
saved_mol = mol;
saved_H = H;
saved_X = X;

if ( isempty(ksopts.relaxmethod) )
   ksopts.relaxmethod = 'fminunc';
end;

if ( isempty(ksopts.relaxtol) )
   ksopts.relaxtol = 1e-4;
end;

xyzlist = mol.xyzlist;
x0 = xyzlist(:);
if (strcmp(lower(ksopts.relaxmethod),'fminunc'))
    % requires MATLAB optimization toolbox
    optimopts = optimoptions('fminunc');
    optimopts.Algorithm = 'quasi-newton';
    optimopts.Diagnostics = 'on';
    if verLessThan('optim', '7.4')
        optimopts.FinDiffRelStep = 0.01;
        optimopts.GradObj = 'on';
        optimopts.MaxFunEvals=500;
        optimopts.MaxIter = 100;
        optimopts.TolX = ksopts.relaxtol/sqrt(3*natoms);
    else
        optimopts.FiniteDifferenceStepSize = 0.01;
        optimopts.SpecifyObjectiveGradient = true;
        optimopts.MaxFunctionEvaluations=500;
        optimopts.MaxIterations = 100;
        optimopts.OptimalityTolerance = ksopts.relaxtol/sqrt(3*natoms);
    end
    
    [x,~] = fminunc(@(x) ksfg(x,mol,ksopts), x0, optimopts);
    xyzlist_new = reshape(x, [natoms 3]);
    mol = set(mol,'xyzlist',xyzlist_new);
elseif (strcmp(lower(ksopts.relaxmethod),'nlcg'))
    % default: nonlinear CG by M. Overton
    if nVarargs == 2
        % output order: [x, f, g, frec, grec, alpharec]
        [x, f, g, frec, grec, ~] = nlcg(varargin{1}, varargin{2});
    else
        pars.nvar = length(x0);
        pars.fgname = 'ksef';
        pars.mol = mol;
        pars.ksopts = ksopts;
        
        options.x0 = x0;
        options.nstart = 1;
        options.maxit = 100;
        options.normtol = ksopts.relaxtol/sqrt(3*natoms);
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
        
        options.strongwolfe = 1;
        options.wolfe1 = 0.1;           % 0 < c1 < c2 < 1/2
        options.wolfe2 = 0.49;          % c2 < 1/2 for NLCG
        options.prtlevel = 1;
        
        % output order: [x, f, g, frec, grec, alpharec]
        [x, f, g, frec, grec, ~] = nlcg(pars, options);
    
    end

    nIonicIter = length(frec{1}(:));
    natoms = sum(mol.natoms);
    
    if nargout == 8
        varargout{1} = f;                           % converged energy value
        varargout{2} = -reshape(g,[natoms 3]);      % converged force value
        varargout{3} = frec;        % history of energy values
        
        % Get the history of force on each atom at every iteration:
        fAtoms = cell(nIonicIter,1);
        
        for n = 1:nIonicIter
            fAtoms{n} = -reshape(grec{1}(:,n),[natoms 3]);
        end
        varargout{4} = fAtoms;
    end 
    
    fprintf('******************************************* \n');
    fprintf('Number of Ionic Iterations to Converge: %3d \n', nIonicIter);
    fprintf('******************************************* \n');
    
    xyzlist_new = reshape(x,[natoms 3]);
    mol = set(mol,'xyzlist',xyzlist_new);
    
elseif (strcmp(lower(ksopts.relaxmethod),'bfgs'))
    % The BFGS package by Michael Overton (part of the
    % HANSO: Hybrid Algorithm for Non-Smooth Optimization (V2.2)
    % http://www.cs.nyu.edu/overton/software/hanso/
    %
    if nVarargs == 2
        % output order: [x, f, d, H, iter, info, X, G, w, fevalrec, xrec, drec, Hrec]
        [x, f, d, ~, nIonicIter, info, ~, ~, ~, frec, ~, drec, ~] = bfgs(varargin{1}, varargin{2});
    else
        % pars (same as NLCG)
        pars.nvar = length(x0);
        pars.fgname = 'ksef';       % same as NLCG (accessed from nlcg1_0 dir)
        pars.mol = mol;
        pars.ksopts = ksopts;
        
        % options
        options.x0 = x0;
        options.nstart = 1;
        options.maxit = 1000;
        options.nvec = 0;           % 0 => full BFGS, else specify m \in [3, 20]
        options.H0 = sparse(eye(pars.nvar));
        options.scale = 1;
        options.normtol = ksopts.relaxtol/sqrt(3*natoms);   % default
        options.strongwolfe = 1;
        options.wolfe1 = 1.0e-4;    % 0 < c1 < c2 < 1
        options.wolfe2 = 0.5;       % c2 = 1/2 is also default
        options.ngrad = 1;          % saves the final gradient
        options.prtlevel = 1;       % default
    
        % Call BFGS
        % output order: [x, f, d, H, iter, info, X, G, w, fevalrec, xrec, drec, Hrec]
        [x, f, d, ~, nIonicIter, info, ~, ~, ~, frec, ~, drec, ~] = bfgs(pars, options);
    end
    
    if nargout == 7
        varargout{1} = f;                           % converged energy value
        varargout{2} = reshape(x,[natoms 3]);       % converged position
        varargout{3} = frec;                        % history of energy values
    end
    
    if nargout == 8
        varargout{1} = f;                           % converged energy value
        varargout{2} = -reshape(d,[natoms 3]);      % converged force value
        varargout{3} = frec;        % history of energy values
        
        % Get the history of force on each atom at every iteration:
        fAtoms = cell(nIonicIter,1);
        
        for n = 1:nIonicIter
            fAtoms{n} = -reshape(drec{n}(:,1),[natoms 3]);
        end
        varargout{4} = fAtoms;
    end 
    
    fprintf('EXIT INFO STATUS: %d\n', info);
    
    fprintf('******************************************* \n');
    fprintf('Number of Ionic Iterations to Converge: %3d \n', nIonicIter);
    fprintf('******************************************* \n');
    
    natoms = sum(mol.natoms);
    xyzlist_new = reshape(x,[natoms 3]);
    mol = set(mol,'xyzlist',xyzlist_new);
    
    
else
    fprintf('Error: invalid optimization method = %s\n', ksopts.relaxmethod);
    fprintf('method must be fminunc, NLCG or BFGS\n');
end
%
% run SCF again to pass out H and X
ksopts.maxscfiter = 100;
ksopts.maxcgiter = 100;
ksopts.cgtol = 1e-7;
ksopts.betamix = 0.1;
ksopts.mixdim = 10;
ksopts.scftol = 1e-6;

[mol,H,X,info]=scf4m(mol,ksopts);

end
%-------------------------------------------------------------------
