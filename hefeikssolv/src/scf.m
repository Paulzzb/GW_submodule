function [mol,H,X,info] = scf(mol,options)
% Self Consistent Field iteration for both hybrid and non-hybrid functionals.
% For non-hybrid functionals, it just calls scf0.
% Its structure is similar to the electrons.f90 of QE.

if (nargin < 2)
    options = setksopt();
end

global verbose;
verbose    = options.verbose;
maxPhiIter = options.maxphiiter;
phiTol = options.phitol;

% Regular SCF Iteration
ishybrid = strcmp(mol.funct,'HSE06');
% Perform non-hybrid calculation first to get X and rho.
if (isempty(options.X0)&&isempty(options.rho0))
    if ishybrid
        fprintf('Regular hybrid SCF for Initialization\n');
	funct = mol.funct;
	mol.funct='PBE';
    else
	fprintf('Regular SCF for Pure DFT\n');
    end
    [mol,H,X,info,options] = scf0(mol,options);
    % The output options is used to pass adaptive cgtol and ev0
    % for davidson diagonalization
    options.X0 = X;
    options.rho0 = H.rho;

    if ishybrid
        mol.funct = funct;
    else
	return
    end
end

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
    fock1 = getExx(X, Vexx, mol);
    [Vexx, mol, options]  = getVexx(X, mol, options);
    fock2 = getExx(X, Vexx, mol);
    dfock = (abs(fock1-fock0)+abs(fock1-fock2))/2;
    info.Etot = info.Etot - fock1 + fock2;
    fock0 = fock2;
    fprintf('Etot(corrected)   = %20.13e\n', info.Etot);
    fprintf('Etot(corrected,Ry)= %20.13e\n', info.Etot*2);
    fprintf('Fock Energy       = %20.13e\n', fock2);
    fprintf('Fock Energy(Ry)    = %20.13e\n', fock2*2);
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

timetot = toc(hsestart);
fprintf('Etot        = %20.13e\n', info.Etot);
fprintf('Etot(Ry)    = %20.13e\n', info.Etot*2);
fprintf('Efock       = %20.13e\n', fock2);
fprintf('Efock(Ry)   = %20.13e\n', fock2*2);
fprintf('Total time  = %20.13e\n', timetot);
end

function [mol,H,X,info,options] = scf0(mol,options)
% SCF Self Consistent Field iteration, both for semiconductor and metal.
%    [mol,H,X,info] = SCF(mol,options) adopts Self Consistent Field (SCF)
%    iteration to find the ground state minimum total energy and the
%    corresponding wave functions. mol is a Molecule object and options is
%    the options for running the SCF. Please read setksopt for detailed
%    information about options. SCF returns the molecule mol with/without
%    force, the Hamiltonian H, the wave functions X, and the information
%    for each iteration in info.

%This file is a merged version of the old scf, scf4m and scf4c.

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
%Tbeta      = 315774.67 / mol.temperature;
%options.cgtol = 1e-9;
%options.cgtol = 1e-2;
%The original value of cgtol is 1e-2, it seems not correct, and the final result is wrong.
%Anyway, 1e-9 is safe.

% Initialize Hamiltonian, Wavefun, and Preconditioners
[mol,H,X,Hprec,nocc] = iterinit(mol,rho,X);

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

vhart = getVhart(mol,rho);
Ecoul_old = getEcoul(mol,rho,vhart);
fprintf('Beging SCF calculation for %s...\n',mol.name);
for iterscf = 1:maxscfiter
    
    fprintf('SCF iter %3d:\n', iterscf);
    options.iterscf = iterscf;
    
    rhoin  = rho;
    vtotin = vtot;
    
    first = (iterscf == 1);
    if first&&isfield(options,'ev0')
        H.eband = options.ev0;
    end

    if iscryst
    	idx = 0;
    	for ik = 1:nkpts
    	    idx = idx(end) + (1:nXcols(ik));
    	    [X{ik}, ev(idx), options] = updateX(mol, H{ik}, X{ik}, Hprec{ik}, options);
            H.eband{ik} = ev(idx);
    	end
        % Perform another Hamiltonian diagonalization for spin-down electrons in spin-unrestricted case
        if nspin == 2
            for ik = 1:nkpts
                idx = idx(end) + (1:nXcols(ik+nkpts));
                [X{ik+nkpts}, ev(idx),options] = updateX(mol, H{ik}, X{ik+nkpts}, Hprec{ik}, options);
                H.eband{ik+nkpts} = ev(idx);
            end
        end
    else
        if nspin == 1||nspin == 4
    	    [X, ev, options] = updateX(mol, H, X, Hprec, options);
            H.eband = ev;
        elseif nspin == 2
            [X{1}, evup,options] = updateX(mol, H, X{1}, Hprec, options);
            [X{2}, evdw,options] = updateX(mol, H, X{2}, Hprec, options);
            ev = [evup; evdw];
            H.eband = ev;
        end
    end
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
    
    if strcmpi(what2mix,'rho')
        if nspin == 1
            %rhoerr = norm(rho(:)-rhoin(:))/norm(rhoin(:));
            rhoerr = 2*rhoerr_qe(mol,rho,rhoin);
        else
            rhoerr = 2*rhoerr_qe(mol,rho,rhoin);
        end
        scferr(iterscf) = rhoerr;
        fprintf('Rel Rho Err     = %20.3e\n',rhoerr);
        
        [rho,dfmat,dvmat,cdfmat] =... 
            potmixing(mol,rhoin,rho,iterscf,mixtype,...
            betamix, dfmat, dvmat, cdfmat, mixdim, brank);
    end
    
    if nspin == 1
        Entropy = getEntropy(ev,mol.efermi,Tbeta,smear)*Tbeta*2;
    else
        Entropy = getEntropy(ev,mol.efermi,Tbeta,smear)*Tbeta;
    end

    if iscryst
        Entropy = Entropy/nkpts;
    end
    % Kinetic energy and some additional energy terms
    if iscryst
        idx = 0;
        for ik = 1:nkpts*(mol.lsda+1)
            idx = idx(end) + (1:nXcols(ik));
            % the weight is multiplied after calculation of Entropy
            ev(idx) = ev(idx)*wks(ik);
        end
    end

    if nspin == 1
        Ekin = 2*sum(ev.*occ);
    else
        Ekin = sum(ev.*occ);
    end

    % ionic and external potential energy was included in Ekin
    % along with incorrect Ecoul and Exc. Need to correct them
    % later;
    Ecor = getEcor(mol, rho, vtot, vion, vext);
    
    % Compute Hartree and exchange correlation energy and potential
    % using the new charge density; update the total potential
    [vhart,vxc,uxc2,rho,uxcsr]=getVhxc(mol,rho);
    
    % Update total potential
    vtot = getVtot(mol, vion, vext, vhart, vxc);
    if strcmpi(what2mix,'pot')
        if nspin == 1
            %vtoterr = norm(vtot(:)-vtotin(:))/norm(vtotin(:));
            vtoterr = rhoerr_qe(mol,vtot,vtotin);
        elseif nspin == 2||nspin == 4
            vtoterr = rhoerr_qe(mol,vtot,vtotin);
        end
        scferr(iterscf) = vtoterr;
        fprintf('Rel Vtot Err    = %20.3e\n',vtoterr);
        
        [vtot,dfmat,dvmat,cdfmat] = ...
            potmixing(mol,vtotin,vtot,iterscf,mixtype,...
            betamix,dfmat,dvmat,cdfmat,mixdim,brank);
    end
    H.vtot = vtot;
    
    % Calculate the potential energy based on the new potential
    Ecoul = getEcoul(mol,rho,vhart);
    Exc   = getExc(mol,rho,uxc2,uxcsr);
    Etot  =  Entropy + Ekin + Eewald + Ealphat + Ecor + Ecoul + Exc;

    if ishybrid
    % Exchange energy is double counted in Ekin, need correction.
    	Exx  = getExx(X,Vexx,mol);
        if nspin == 2 || nspin == 4
            Exx = Exx/2;
        end
	Etot = Etot - Exx;
    end
    Etotvec(iterscf) = Etot;
    
    % Convergence check
    fprintf('Total Energy    = %20.13e\n', Etot);
    [cvg,resfro] = reportconverge(H,X,iterscf,maxscfiter, ...
        scferr(iterscf),scftol,verbose);
    if cvg
        info.converge = true;
        break;
    end
    deltaE = abs(Ecoul - Ecoul_old);
    Ecoul_old = Ecoul;
    options.cgtol = min(options.cgtol,0.01*deltaE/max(1.0,nocc));
    options.cgtol = max(options.cgtol,1e-10);
end

if (nspin == 4 && ~domag) || nspin == 1
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


if iscryst
    info.Eigvals = reshape(ev,mol.nbnd,mol.nkpts*(mol.lsda+1));
    options.ev0 = H.eband;
else
    info.Eigvals = ev;
    info.rho=rho;
    options.ev0 = ev;
end

if domag
    Mag = calmag(mol,rho);
end

info.Etotvec = Etotvec(1:iterscf);
info.SCFerrvec = scferr(1:iterscf);
info.Etot = Etot;
info.Eentropy = Entropy;
info.Efree = Etotvec(end) - Tbeta*Entropy;
info.Eoneelectron = Ekin+Ecor;
info.Ehart = Ecoul;
info.Exc = Exc;
info.Eewald = Eewald;
info.Efermi = mol.efermi;
if domag
    info.Mag = Mag;
end

timetot = toc(scfstart);
fprintf('Etot            = %20.13e\n', Etot);
fprintf('Entropy         = %20.13e\n', Entropy);
fprintf('Ekin            = %20.13e\n', Ekin);
fprintf('Eewald          = %20.13e\n', Eewald);
fprintf('Ealphat         = %20.13e\n', Ealphat);
fprintf('Ecor            = %20.13e\n', Ecor);
fprintf('Ehart           = %20.13e\n', Ecoul);
fprintf('Exc             = %20.13e\n', Exc);

if mol.abnd
    %ha2ev = 4.3597447222071e-18/1.602176634e-19;
    %nocc_max = find(occ>eps,1,'last');
    fprintf('Efermi         = %20.13e\n',mol.efermi);
 %   fprintf('HO energy       = %20.13e\n',ev(nocc_max));
%    fprintf('LO energy       = %20.13e\n',ev(nocc_max+1));
end

if force
    fprintf('----------------Forces for atoms----------------\n');
    force = mol.xyzforce;
    for it = 1:numel(mol.alist)
        fprintf('%5s : [%20.13e %20.13e %20.13e]\n',mol.atoms(mol.alist(it)).symbol,...
            force(it,1),force(it,2),force(it,3));
    end
end

if domag
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
