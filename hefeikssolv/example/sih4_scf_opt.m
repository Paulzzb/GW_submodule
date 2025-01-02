clc;
clear all;
close all;
%
% Set timer
%
tstart  = cputime;
%
%
% construct the SiH4 (Silane) molecule 
%
kssolvpptype('pz-vbc');
%
% 1. construct atoms
%
a1 = Atom('Si');
a2 = Atom('H');
atomlist = [a1 a2 a2 a2 a2];
%
% 2. set up supercell
%
C = 10*eye(3);
%
% 3. define the coordinates the atoms 
%
coefs = [
 0.02    0.0      0.0
 0.161   0.161    0.161
-0.161  -0.161    0.161
 0.161  -0.161   -0.161
-0.161   0.161   -0.161
];
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
               'ecut',12.5, 'name','SiH4' );


scfOpts = setksopt();
scfOpts.maxscfiter = 100;
scfOpts.betamix = 0.1;
scfOpts.mixdim = 10;
scfOpts.scftol = 1.0E-6;
% scfOpts = setksopt(scfOpts,'mixtype','pulay','eigmethod','eigs');
scfOpts = setksopt(scfOpts,'mixtype','pulay','eigmethod','lobpcg');

[mol, H, X, info]=scf4m(mol, scfOpts);

pxyzlist = xyzlist + randn(size(xyzlist))/5;

molper = mol;
molper = set(molper,'xyzlist',pxyzlist);

strOptMethod = 4;       % 5 => FIRE; 3 => nlcg(); 6 => hybrid; 2 => NLCG Amartya

if strOptMethod == 2
    pars.jmax = 6;
    pars.n = 5;
    pars.sigma0 = 0.02;
    
    opts.MAXITER = 50;
    opts.TOL1 = 1.0E-4;
    opts.TOL2 = 1.0E-3;

    % Calling relaxatom with NLCG_mod
    [molopt, Hopt, Xopt, infoopt, Etot, FORCE] = relaxatoms(molper, ...
                                                 strOptMethod, H, X, scfOpts, ...
                                                 pars, opts);
    
    % post processing for NLCG
    nIonicIter = length(Etot);
    nAtoms = size(molopt.xyzlist,1);
    
    maxAbsForceAtom = zeros(nIonicIter,nAtoms);
    maxAbsForce = zeros(nIonicIter,1);
    forceNorm = zeros(nIonicIter,1);
    
    for n = 1:nIonicIter
        f = [FORCE{1,n}];
        forceNorm(n) = norm(reshape(f, 1, []));
        
        for j = 1:nAtoms
            % max of the abs. val of force componenet on the jth atom at
            % nth iter
            
            maxAbsForceAtom(n,j) = max(f(j,:));
        end
        maxAbsForce(n) = max(maxAbsForceAtom(n,:));
    end

elseif strOptMethod == 3     % nlcg (Overton's version) 
    x0 = pxyzlist(:);
    pars.nvar = length(x0);
    pars.fgname = 'ksef';
    pars.mol = molper;
    pars.ksopts = scfOpts;
    
    options.x0 = x0;
    options.nstart = 1;
    options.maxit = 1000;
    options.normtol = 1.0e-4;
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
     
    % Call nlcg
    [molopt, Hopt, Xopt, infoopt, eAtOpt, fAtOpt, Etot, FORCE] = ...
                                relaxatoms(molper, strOptMethod, H, X, ...
                                pars, options);
    
    % post processing for nlcg
    nIonicIter = length(Etot{:});
    Etot = [Etot{:}];
    nAtoms = size(molopt.xyzlist,1);
    
    maxAbsForceAtom = zeros(nIonicIter,3);
    maxAbsForce = zeros(nIonicIter,1);
    forceNorm = zeros(nIonicIter,1);
    
    for n = 1:nIonicIter
        forceNorm(n) = norm(reshape([FORCE{n}], 1, []));
        
        for j = 1:nAtoms
            % max of the abs. val of force componenet on the jth atom at
            % nth iter
            maxAbsForceAtom(n,j) = max(abs(FORCE{n}(j,:)));
        end
        maxAbsForce(n) = max(maxAbsForceAtom(n,:));
    end
    
elseif strOptMethod == 4     % bfgs (Overton's code)
    %
    x0 = pxyzlist(:);
    pars.nvar = length(x0);
    pars.fgname = 'ksef';                       % same as NLCG (accessed from nlcg1_0 dir)
    pars.mol = molper;
    pars.ksopts = scfOpts;
    %
    options.x0 = x0;
    options.nstart = 1;
    options.maxit = 1000;
    options.nvec = 0;                           % 0 => full BFGS, else specify \in [3, 20] 
    options.H0 = sparse(eye(pars.nvar));        % explore
    options.scale = 1;
    options.normtol = 1.0e-4;   % 1.0e-6 is default but hard to reach
    options.strongwolfe = 0;
    options.wolfe1 = 1.0e-4;    % 0 < c1 < c2 < 1 (it's different from nlcg)
    options.wolfe2 = 0.5;       % c2 = 1/2 is also default
    options.ngrad = 1;          % saves the final gradient
    options.prtlevel = 1;       % default
   
    % Calling relaxatom with BFGS
    [molopt, Hopt, Xopt, infoopt, eAtOpt, posAtOpt, Etot, FORCE] = ...
                                relaxatoms(molper, strOptMethod, H, X, ...
                                pars, options);

    nAtoms = size(molopt.xyzlist,1);
    
    % post processing for BFGS
    for n = 1:size(Etot,2)
        Etot{n} = Etot{n}(end);
    end
    Etot = cell2mat(Etot);
    
    nIonicIter = size(Etot,2);
    natoms = size(molopt.xyzlist,1);
    
    forceNorm = zeros(nIonicIter,1);
    
    for n = 1:nIonicIter
        forceNorm(n) = norm(FORCE{n});
        maxAbsForce(n) = max(max(abs(FORCE{n})));
    end                                        
    
elseif strOptMethod == 5        % FIRE
    mass = 4.0;     % MASS = 4 unit
    
    pars.fDec = 0.5;
    pars.fInc = 1.1;
    % for SiH4 nMin = 10 seems to be performing BEST!
    % pars.nMin = 10;
    pars.nMin = 5;
    pars.alphaStart = 0.1;
    pars.fAlpha = 0.99;
    
    opts.dt = 41.3413745758/4.0;
    opts.MAXITER = 250;
    opts.TOL = 1.0E-4;

    % Calling relaxatom with FIRE
    [molopt, Hopt, Xopt, infoopt, POS, FORCE, Etot, INFO, exitFlag] = ...
                                relaxatoms(molper, strOptMethod, H, X, ...
                                           mass, pars, opts, scfOpts);
    % post processing for FIRE
    nIonicIter = length(Etot);
    nAtoms = size(molopt.xyzlist,1);
    
    maxAbsForceAtom = zeros(nIonicIter,nAtoms);
    maxAbsForce = zeros(nIonicIter,1);
    forceNorm = zeros(nIonicIter,1);
    
    for n = 1:nIonicIter
        forceNorm(n) = norm(reshape(FORCE(:,(n-1)*3+1:n*3), 1, []));
        
        for j = 1:nAtoms
            % max of the abs. val of force componenet on the jth atom at
            % nth iter
            maxAbsForceAtom(n,j) = max(abs(FORCE(j,(n-1)*3+1:n*3)));
%             max(abs(FORCE{n}(j,:)));
        end
        maxAbsForce(n) = max(maxAbsForceAtom(n,:));
    end
    
elseif strOptMethod == 6        % hybrid FIRE
    forceTol = 1.0E-4;
    mass = 10.0;     % MASS = 4 unit
    
    pars.fDec = 0.5;
    pars.fInc = 1.1;
    pars.nMin = 5;       % for SiH4 nMin = 10 seems to be performing BEST!
    pars.alphaStart = 0.1;
    pars.fAlpha = 0.99;
    
%     opts.dt = 41.3413745758/4.0;
    opts.dt = 41.3413745758;
    opts.MAXITER = 10;
    opts.TOL = forceTol;

    fireinput = struct('mass',mass, 'pars',pars, 'opts',opts);
    
    otheroptimizer = 'nlcg_overton';
    optimizeropts = struct('MAXITER',100, 'TOL',forceTol);
     
    
    % Calling relaxatom with hybrid with FIRE
    [molopt, Hopt, Xopt, infoopt, FORCE, Etot] = relaxatoms(molper, ...
                                                 strOptMethod, H, X, ...
                                                 otheroptimizer, fireinput, ...
                                                 optimizeropts, scfOpts);
    
    nIonicIter = length(Etot);
    nAtoms = size(molopt.xyzlist,1);
    
    maxAbsForceAtom = zeros(nIonicIter,3);
    maxAbsForce = zeros(nIonicIter,1);
    forceNorm = zeros(nIonicIter,1);
    
    for n = 1:nIonicIter
        forceNorm(n) = norm(reshape([FORCE{n}], 1, []));
        
        for j = 1:nAtoms
            % max of the abs. val of force componenet on the jth atom at
            % nth iter
            maxAbsForceAtom(n,j) = max(abs(FORCE{n}(j,:)));
        end
        maxAbsForce(n) = max(maxAbsForceAtom(n,:));
    end        
    
    

end
    
% 
% Plot Total energy:
%     
width = 6;     % Width in inches
height = 4;    % Height in inches
alw = 5.0;     % AxesLineWidth
fsz = 18;      % Fontsize
lw = 2.0;      % LineWidth

figure
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*1000, height*1000]);
set(gca, 'FontSize', fsz, 'LineWidth', alw);
set(get(gca,'xlabel'),'FontSize', fsz, 'FontWeight', 'Bold');
set(get(gca,'ylabel'),'FontSize', fsz, 'FontWeight', 'Bold');
set(get(gca,'title'),'FontSize', fsz, 'FontWeight', 'Bold');
hold on;
% 
plot(1:nIonicIter, Etot, 'r', 'LineWidth', 1.25)

% 
xlabel('Iteration Number','Interpreter','tex');
ylabel('Total Energy (Ha)','Interpreter','tex');
box on; 
axis square;
set(gca,'LineWidth',lw);
set(gca,'FontSize',fsz);
set(gca,'FontWeight','Bold');
set(gcf,'color','w');
% 
% Plot Force (max(max)):
%  
figure
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*1000, height*1000]);
set(gca, 'FontSize', fsz, 'LineWidth', alw);
set(get(gca,'xlabel'),'FontSize', fsz, 'FontWeight', 'Bold');
set(get(gca,'ylabel'),'FontSize', fsz, 'FontWeight', 'Bold');
set(get(gca,'title'),'FontSize', fsz, 'FontWeight', 'Bold');
hold on;
% 
plot(1:nIonicIter, log10(maxAbsForce), 'r', 'LineWidth', 1.25)
% 
xlabel('Iteration Number','Interpreter','tex');
ylabel('Max. Force (Ha/Bohr)','Interpreter','tex');
box on; 
axis square;
set(gca,'LineWidth',lw);
set(gca,'FontSize',fsz);
set(gca,'FontWeight','Bold');
set(gcf,'color','w');
% 
% Plot Force (norm):
%  
figure
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*1000, height*1000]);
set(gca, 'FontSize', fsz, 'LineWidth', alw);
set(get(gca,'xlabel'),'FontSize', fsz, 'FontWeight', 'Bold');
set(get(gca,'ylabel'),'FontSize', fsz, 'FontWeight', 'Bold');
set(get(gca,'title'),'FontSize', fsz, 'FontWeight', 'Bold');
hold on;
% 
plot(1:nIonicIter, log10(forceNorm), 'r', 'LineWidth', 1.25)
% 
xlabel('Iteration Number','Interpreter','tex');
ylabel('Norm of Force (Ha/Bohr)','Interpreter','tex');
box on; 
axis square;
set(gca,'LineWidth',lw);
set(gca,'FontSize',fsz);
set(gca,'FontWeight','Bold');
set(gcf,'color','w');

timetot = cputime - tstart;
fprintf('*************************************************************** \n');
fprintf('\n');
fprintf(' Time to complete the FPMD Simulation = %20.3e s\n', timetot);
fprintf('\n');
fprintf('*************************************************************** \n');
