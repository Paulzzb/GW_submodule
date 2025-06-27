function [mol,H,X,info] = scf(mol,options)
% SCF Self Consistent Field iteration for both hybrid and non-hybrid functionals.
% For non-hybrid functionals, it just calls scf0.
%    [mol,H,X,info] = scf(mol,options) find the ground state minimum
%    total energy and the corresponding wave functions for mol, mol can be
%    Molecule or Crystal. 
%
%    The framwork is the same as QE which is more efficient to reach convergence
%    than the original implemention.
%
%    See also nscf4c, dcm, trdcm.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if (nargin < 2)
    options = setksopt();
end

global verbose;
verbose    = options.verbose;
maxPhiIter = options.maxphiiter;
phiTol = options.phitol;

% Regular SCF Iteration
ishybrid = mol.exx_sr ~= 0 || mol.exx_lr ~= 0;
% Perform non-hybrid calculation first to get X and rho.
if ~ishybrid
    fprintf('Regular SCF for Pure DFT\n');
    [mol,H,X,info,~] = scf0(mol,options);
    return
end

Etotvec = []; SCFerrvec = [];
if (isempty(options.X0)||isempty(options.rho0))
    fprintf('Regular hybrid SCF for Initialization\n');
    omega = mol.omega; exx_sr = mol.exx_sr; exx_lr = mol.exx_lr; dft_c = mol.dft_c;
    mol.omega = 0; mol.exx_sr = 0; mol.exx_lr = 0; mol.dft_c = 1;
    [mol,H,X,info,options] = scf0(mol,options);
    Etotvec = info.Etotvec; SCFerrvec = info.SCFerrvec;
    
    options.X0 = X;
    options.rho0 = H.rho;
	mol.omega = omega; mol.exx_sr = exx_sr; mol.exx_lr = exx_lr; mol.dft_c = dft_c;
end

% Hybrid SCF Calculation
hsestart=tic;
fprintf('Beging Hybrid SCF calculation for %s...\n',mol.name);
for iterphi = 1:maxPhiIter
    fprintf('Phi iter %3d:\n', iterphi);
	hsestep=tic;

    if iterphi == 1
        options.ishybrid = true;
        options = getExxgkk(mol,options);
        [Vexx, mol, options]  = getVexx(options.X0, mol, options);
        fock0 = getExx(options.X0, Vexx, mol);
    else
        options.X0 = X;
        options.rho0 = H.rho;
    end
	options.vexx = Vexx;

    [mol,H,X,info,options] = scf0(mol, options);
    Etotvec = [Etotvec; info.Etotvec]; SCFerrvec = [SCFerrvec; info.SCFerrvec];

    fock1 = getExx(X, Vexx, mol);
	[Vexx, mol, options]  = getVexx(X, mol, options);
	fock2 = getExx(X, Vexx, mol);
    if mol.nspin == 1
        dfock =4*abs( fock1 - 0.5*(fock0+fock2));
    elseif mol.nspin == 2||mol.nspin == 4
        dfock =2*abs( fock1 - 0.5*(fock0+fock2));
    end
    info.Etot = info.Etot - fock1 + fock2;
    fock0 = fock2;
    fprintf('Etot(corrected)   = %20.13e\n', info.Etot);
    fprintf('Etot(corrected,Ry)= %20.13e\n', info.Etot*2);
    fprintf('Fock Energy       = %20.13e\n', fock2/2);
    fprintf('Fock Energy(Ry)    = %20.13e\n', fock2);
    fprintf('dfock             = %20.13e\n', dfock);
    fprintf('HSE step time = %20.3e\n', toc(hsestep));
    fprintf('======================================\n');
    if (dfock < phiTol)
        break
    end
end

if (dfock < phiTol)
	fprintf('HSE convergence is reached.\n');
else
	fprintf('######################################\n');
	fprintf('Warning: HSE not converge!\n');
end

info.Efock  = fock2;
info.Etotvec = Etotvec;
info.SCFerrvec = SCFerrvec;

timetot = toc(hsestart);
fprintf('Etot        = %20.13e\n', info.Etot);
fprintf('Etot(Ry)    = %20.13e\n', info.Etot*2);
fprintf('Efock       = %20.13e\n', fock2);
fprintf('Efock(Ry)   = %20.13e\n', fock2*2);
fprintf('Total time  = %20.13e\n', timetot);
end

function [mol,H,X,info,options] = scf0(mol,options)
% SCF0 SCF Self Consistent Field iteration, both for semiconductor and metal.
%    [mol,H,X,info] = SCF(mol,options) adopts Self Consistent Field (SCF)
%    iteration to find the ground state minimum total energy and the
%    corresponding wave functions. mol is a Molecule/Crystal object and options 
%    is the options for running the SCF. Please read setksopt for detailed
%    information about options. SCF returns the molecule mol with/without
%    force, the Hamiltonian H, the wave functions X, and the information
%    for each iteration in info.
%    This file is a merged version of the old scf, scf4m and scf4c.
% 
%    See also nscf4c, dcm, trdcm.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if (nargin < 2)
    options = setksopt();
end

scfstart  = tic;

% Initialize input variables
global verbose;
force      = options.force;
maxscfiter = options.maxscfiter;
scftol     = options.scftol;
what2mix   = options.what2mix;
mixtype    = options.mixtype;
mixdim     = options.mixdim;
betamix    = options.betamix;
brank      = options.brank;
X          = options.X0;
rho        = options.rho0;

iscryst    = isa(mol,'Crystal');
ishybrid   = options.ishybrid;
nspin      = mol.nspin;
domag      = mol.domag;
smear      = mol.smear;
Tbeta      = mol.temperature*8.6173324e-5/13.6;

% Initialize Hamiltonian, Wavefun, and Preconditioners
[mol,H,X,Hprec,nocc] =  iterinit(mol,rho,X);

% calculate Ewald and Ealphat
Eewald     = getEewald(mol);
Ealphat    = getEalphat(mol);

vion       = H.vion;
vext       = H.vext;
vtot       = H.vtot;
rho        = H.rho;

% Initialize output variables
Etotvec    = zeros(maxscfiter,1);
scferr     = zeros(maxscfiter,1);
% three entries stored for mixing
dfmat      = [];
dvmat      = [];
cdfmat     = [];

if iscryst
    nkpts      = mol.nkpts;
	wks        = mol.wks;
    nXcols     = ncols(X);
    ev         = zeros(sumel(nXcols),1);
end

if ishybrid
	Vexx       = options.vexx;
	H.vexx     = Vexx;
end

fprintf('Beging SCF calculation for %s...\n',mol.name);
info.converge = false;
for iterscf = 1:maxscfiter
    
    fprintf('SCF iter %3d:\n', iterscf);
    
    options.iterscf = iterscf;
    if iterscf > 1
        if iterscf == 2
            options.cgtol = 1e-2;
        end
        % Dynamic threshold is used to reduce time for DIAG
        options.cgtol = min(options.cgtol, 0.1*scferr(iterscf-1)/max(1,mol.nel));
        % iterative diagonalization may become unstable if threshold is too small
        options.cgtol = max(options.cgtol, 1e-13); 
    end
    rhoin  = rho;
    vtotin = vtot;
    
    first = (iterscf == 1);
    if first&&isfield(options,'ev0')
        H.eband = options.ev0;
    end
    rediag = true;
    
    while rediag 
        diagerr_min = 0;
        if first
            diagerr_min = options.cgtol*max(1,mol.nel);
        end
        
        diagt = tic;
        if iscryst       
            idx = 0;
            for ik = 1:nkpts
                idx = idx(end) + (1:nXcols(ik));
                [X{ik}, ev(idx),options] = updateX(mol, H{ik}, X{ik}, Hprec{ik}, options);
                H.eband{ik} = ev(idx); 
            end   
            if nspin == 2
                for ik = 1:nkpts
                    idx = idx(end) + (1:nXcols(ik+nkpts));
                    [X{ik+nkpts}, ev(idx),options] = updateX(mol, H{ik}, X{ik+nkpts}, Hprec{ik}, options);
                    H.eband{ik+nkpts} = ev(idx);
                end
            end 
        else
            if nspin == 1||nspin == 4
                [X, ev,options] = updateX(mol, H, X, Hprec, options);
                H.eband = ev;
            elseif nspin == 2
                [X{1}, evup,options] = updateX(mol, H, X{1}, Hprec, options);
                [X{2}, evdw,options] = updateX(mol, H, X{2}, Hprec, options);
                ev = [evup; evdw];
                H.eband = ev;
            end
        end
        fprintf('Diag time = %20.5e\n',toc(diagt));

        [occ,mol.efermi] = getocc(ev,nocc,Tbeta,smear);      
        if iscryst
            idx = 0;
            for ik = 1:nkpts*(mol.lsda+1)
                idx = idx(end) + (1:nXcols(ik));
                % the weight is multiplied after calculation of Entropy
                % ev(idx) = ev(idx)*wks(ik);
                X{ik}.occ = occ(idx);
            end
        else
            if nspin == 1||nspin == 4
                X.occ = occ;
            elseif nspin == 2
                X{1}.occ = occ(1:mol.nbnd);
                X{2}.occ = occ(mol.nbnd+1:end);
            end
         end

        % Update density function rho
        rho = getcharge(mol,X,occ);
        H.rho = rho;
        
        % ionic and external potential energy was included in Ekin
        % along with incorrect Ecoul and Exc. Need to correct them
        % later;
        if strcmpi(what2mix,'rho')
            rhoerr = 2*CalculateError(mol,rho,rhoin);
            scferr(iterscf) = rhoerr;
            fprintf('Rel Rho Err     = %20.3e\n',rhoerr);
            [cvg,resfro] = reportconverge(H,X,iterscf,maxscfiter, ...
             scferr(iterscf),scftol,verbose);
            
            % mix rho when convergence not reached
            if (~cvg)&&(rhoerr>=diagerr_min)
            	[rho,dfmat,dvmat,cdfmat] = potmixing(mol,rhoin,rho,...
                	iterscf,mixtype, betamix, ...
                    dfmat, dvmat, cdfmat, mixdim, ...
                    brank);
            end
        end
               
        if first && strcmp(what2mix,'rho')
            % only restart diagnolization when rho is mixed
            first = false;
            if rhoerr<diagerr_min
                fprintf("The threshold of diagonalization is too large, restart diagnolization\n");
                options.cgtol = 0.1*rhoerr/max(1,mol.nel);
                rediag = true;
            else
                rediag = false;
            end
        else
            rediag = false;
        end
    end
    
    Ecor = getEcor(mol, rho, vtot, vion, vext);
    
    if mol.nspin == 1
        Entropy = getEntropy(ev,mol.efermi,Tbeta,smear)*Tbeta*2;
    else
        Entropy = getEntropy(ev,mol.efermi,Tbeta,smear)*Tbeta;
    end

    if iscryst
        Entropy = Entropy/nkpts;
    end
    % Kinetic energy and some additional energy terms
    if iscryst
        ev_weight = zeros(size(ev));
        idx = 0;
        for ik = 1:nkpts*(mol.lsda+1)
            idx = idx(end) + (1:nXcols(ik));
            % the weight is multiplied after calculation of Entropy
            ev_weight(idx) = ev(idx)*wks(ik);
        end
        Ekin = sum(ev_weight.*occ);
    else
        Ekin = sum(ev.*occ);
    end
    
    if mol.nspin == 1
        Ekin = 2*Ekin;
    end
    
    % ionic and external potential energy was included in Ekin
    % along with incorrect Ecoul and Exc. Need to correct them
    % later;
    
    % Compute Hartree and exchange correlation energy and potential
    % using the new charge density; update the total potential
    [vhart,vxc,exc,rho]=getVhxc(mol,rho);
    
    % Update total potential
    vtot = getVtot(mol, vion, vext, vhart, vxc);
    if strcmpi(what2mix,'pot')
        vtoterr = CalculateError(mol,vtot,vtotin);
        scferr(iterscf) = vtoterr;
        fprintf('Rel Vtot Err    = %20.3e\n',vtoterr);
        [cvg,resfro] = reportconverge(H,X,iterscf,maxscfiter, ...
            scferr(iterscf),scftol,verbose);
        % mix potential when convergence not reached
        if ~cvg
                [vtot,dfmat,dvmat,cdfmat] = ...
                    potmixing(mol,vtotin,vtot,iterscf,mixtype,...
                    betamix,dfmat,dvmat,cdfmat,mixdim,brank);
        end
    end
    H.vtot = vtot;
    
    % Calculate the potential energy based on the new potential
    Ecoul = getEcoul(mol,rho,vhart);
    
    Exc   = getExc(mol,exc);
    Etot  =  Entropy + Ekin + Eewald + Ealphat + Ecor + Ecoul + Exc;
	if ishybrid
		% Exchange energy is double counted in Ekin, need correction.
        if mol.nspin == 1
            Exx  = getExx(X,Vexx,mol);
        else
            Exx  = getExx(X,Vexx,mol)/2;
        end
        Etot = Etot - Exx;
	end
    Etotvec(iterscf) = Etot;
    
    % Convergence check
    fprintf('Total Energy    = %20.13e\n', Etot);
    if cvg
        info.converge = true;
        break;
    end
end

if mol.nspin == 4 && ~mol.domag
    H.dv = vtot - vtotin;
else
    H.dv = rho_mix(mol,-1,vtot,vtotin); % V(out)- V(in) used to correct the forces
end
 
if iscryst
	X = assignoccs(X,occ);
else
    if nspin == 1||nspin == 4
        X.occ = occ;
    elseif nspin == 2
        X{1}.occ = occ(1:mol.nbnd);
        X{2}.occ = occ(mol.nbnd+1:end);
    end
end

if force
    mol.xyzforce = getFtot(mol,H,X,rho);
    info.Force = mol.xyzforce;
end

info.Eigvals = ev;
info.rho = rho;
options.ev0 = H.eband;

if domag
    Mag = calmag(mol,rho);
end

info.Etotvec = Etotvec(1:iterscf);
info.SCFerrvec = scferr(1:iterscf);
info.Etot = Etot;
info.Entropy = Entropy;
info.Efree = Etotvec(end) - Tbeta*Entropy;
info.Eoneelectron = Ekin + Ecor;
info.Ehart = Ecoul;
info.Exc = Exc;
info.Eewald = Eewald;
info.Efermi = mol.efermi;
if domag
    info.Mag = Mag;
end
% output
timetot = toc(scfstart);
fprintf('Etot            = %20.13e\n', Etot);
fprintf('Entropy         = %20.13e\n', Entropy);
fprintf('Ekin            = %20.13e\n', Ekin);
fprintf('Eewald          = %20.13e\n', Eewald);
fprintf('Ealphat         = %20.13e\n', Ealphat);
fprintf('Ecor            = %20.13e\n', Ecor);
fprintf('Ehart           = %20.13e\n', Ecoul);
fprintf('Exc             = %20.13e\n', Exc);
fprintf('Efermi          = %20.13e\n', mol.efermi);

if mol.temperature > 0
    nocc_max = find(occ>eps,1,'last');
    fprintf('HO energy       = %20.13e\n',ev(nocc_max));
    fprintf('LO energy       = %20.13e\n',ev(nocc_max+1));
end
 
if force
    fprintf('----------------Forces for atoms----------------\n');
    force = mol.xyzforce;
    for it = 1:numel(mol.alist)
        fprintf('%5s : [%20.13e %20.13e %20.13e]\n',mol.atoms(mol.alist(it)).symbol,...
            force(it,1),force(it,2),force(it,3));
    end
end

if mol.domag
    fprintf('----------------Magnetization----------------\n');
    if mol.lsda
        fprintf('Magtot          = %20.13e\n', Mag.totmag);
        fprintf('Magabs         = %20.13e\n', Mag.absmag);
    elseif mol.noncolin
        fprintf('Magabs          = %20.13e\n', Mag.absmag);
        fprintf('Magtot          = [%20.13e %20.13e %20.13e]\n', Mag.mx, Mag.my, Mag.mz);
    end
end

fprintf('--------------------------------------\n');
fprintf('Total time used = %20.3e\n', timetot);
fprintf('||HX-XD||_F     = %20.3e\n', resfro);
fprintf('Etot(Ry)        = %20.13e\n', Etot*2);
fprintf('Entropy(Ry)     = %20.13e\n', Entropy*2);
fprintf('Eoneelectron(Ry)= %20.13e\n', (Ekin+Ecor)*2);
fprintf('Ehart(Ry)       = %20.13e\n', Ecoul*2);
fprintf('Exc(Ry)         = %20.13e\n', Exc*2);
fprintf('Eewald(Ry)      = %20.13e\n', Eewald*2);

end
