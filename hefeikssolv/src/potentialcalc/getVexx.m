function [Vexxf, mol, options] = getVexx(X, mol, options)
% Purpose:
%    Computes exact exchange operator
%
% Input:
%    X  --- Wavefunction
%    mol --- Molecule information
%
% Ouptut:
%    Vexx --- Exact exchange operator
%
nspin = mol.nspin;
vexxt = tic;
if isa(mol,'Crystal')
    switch options.exxmethod
        case 'normal'
            F = KSFFT(mol);
            if isfield(options,'scfX') && ~isempty(options.scfX)
                X = options.scfX;
            end
            if nspin == 1
                realwavefunction.psi = zeros(mol.n1*mol.n2*mol.n3,ncols(X{1}),X.nkpts);
                realwavefunction.occ = zeros(size(X{1}.occ,1),X.nkpts);
                realwavefunction.nocc_max = zeros(X.nkpts,1);
                realwavefunction.wks = X.wks;
                for ik=1:X.nkpts
                    realwavefunction.psi(:,:,ik)= F'*X{ik}.psi;
                    realwavefunction.occ(:,ik) = X{ik}.occ;
                    realwavefunction.nocc_max(ik) = find(X{ik}.occ>eps,1,'last');
                end
                Vexxf = @(Psi)Vexx(Psi, realwavefunction, mol, options);
            elseif nspin == 2
                realwavefunction = cell(2,1);
                for is = 1:2
                    realwavefunction{is}.psi = zeros(mol.n1*mol.n2*mol.n3,ncols(X{1}),X.nkpts);
                    realwavefunction{is}.occ = zeros(size(X{1}.occ,1),X.nkpts);
                    realwavefunction{is}.nocc_max = zeros(X.nkpts,1);
                    realwavefunction{is}.wks = X.wks;
                    for ik=1:X.nkpts
                        realwavefunction{is}.psi(:,:,ik)= F'*X{ik+(is-1)*X.nkpts}.psi;
                        realwavefunction{is}.occ(:,ik) = X{ik+(is-1)*X.nkpts}.occ;
                        realwavefunction{is}.nocc_max(ik) = find(X{ik+(is-1)*X.nkpts}.occ>eps,1,'last');
                    end
                end
                Vexxf = @(Psi)Vexx(Psi, realwavefunction, mol, options);
            elseif nspin == 4
                npw = length(X{1}.idxnz);
                realwavefunction.psi = cell(2,1);
                realwavefunction.psi{1} = zeros(mol.n1*mol.n2*mol.n3,ncols(X{1}),X.nkpts);
                realwavefunction.psi{2} = zeros(mol.n1*mol.n2*mol.n3,ncols(X{1}),X.nkpts);
                realwavefunction.occ = zeros(size(X{1}.occ,1),X.nkpts);
                realwavefunction.nocc_max = zeros(X.nkpts,1);
                realwavefunction.wks = X.wks;
                for ik=1:X.nkpts
                    realwavefunction.psi{1}(:,:,ik)= F'*X{ik}.psi(1:npw,:);
                    realwavefunction.psi{2}(:,:,ik)= F'*X{ik}.psi(npw+1:end,:);
                    realwavefunction.occ(:,ik) = X{ik}.occ;
                    realwavefunction.nocc_max(ik) = find(X{ik}.occ>eps,1,'last');
                end
                Vexxf = @(Psi)Vexx(Psi, realwavefunction, mol, options);
            end
        case {'ACE','ace'}    
            if ~mol.lsda
                [mol.zeta, options] = ACE_k(X,[],mol,options);
            else
                mol.zeta = cell(2,1);
                X = splitX_to_updw(X);
                [mol.zeta{1}, options] = ACE_k(X{1},[],mol,options);
                [mol.zeta{2}, options] = ACE_k(X{2},[],mol,options);
            end      
            Vexxf = @(Psi)Vexx(Psi, [], mol, options);
        case {'default','qrcp','kmeans'}
            %Save zeta in mol, thus it can be exported by scf.m and used in band calculation.
	    if isfield(options,'scfX') && ~isempty(options.scfX)
                [mol.zeta, options] = ACEISDF_k(X,options.scfX,mol,options);
	    else
                if nspin == 1
                    [mol.zeta, options] = ACEISDF_k(X,[],mol,options);
                elseif nspin == 2
                    mol.zeta = cell(2,1);
                    X = splitX_to_updw(X);
                    [mol.zeta{1}, options] = ACEISDF_k(X{1},[],mol,options);
                    [mol.zeta{2}, options] = ACEISDF_k(X{2},[],mol,options);
                elseif nspin == 4
                    [mol.zeta, options] = ACEISDF_k_nlc(X,[],mol,options);
                end
	    end
            Vexxf = @(Psi)Vexx(Psi, [], mol, options);
        case 'nscf'
            if nspin == 1
                mol.zeta=GACE(mol);
            elseif nspin == 2
                zeta = cell(2,1);
                zeta{1} = GACE(mol,mol.zeta{1});
                zeta{2} = GACE(mol,mol.zeta{2});
                mol.zeta = zeta;
            elseif nspin == 4
                npw = size(mol.zeta,1)/2;
                zeta_up = GACE(mol,mol.zeta(1:npw,:,:));
                zeta_dw = GACE(mol,mol.zeta(npw+1:end,:,:));
                rank = max(size(zeta_up,2),size(zeta_dw,2));
                nk = size(zeta_up,3);
                mol.zeta = zeros(npw*2,rank,nk);
                mol.zeta(1:npw,1:size(zeta_up,2),:) = zeta_up;
                mol.zeta(npw+1:end,1:size(zeta_dw,2),:) = zeta_dw;
            end
            Vexxf = @(Psi)Vexx(Psi, [], mol, options);
    end
else
    switch options.exxmethod
        case 'normal'
            F = KSFFT(mol);
            if nspin == 1
                nocc_max =  find(X.occ>eps,1,'last');
                realwavefunction.psi = F' * X.psi(:,1:nocc_max);
                realwavefunction.occ = X.occ(1:nocc_max);
            elseif nspin == 2
                realwavefunction = cell(2,1);
                for is = 1:2
                    nocc_max =  find(X{is}.occ>eps,1,'last');
                    realwavefunction{is}.psi = F' * X{is}.psi(:,1:nocc_max);
                    realwavefunction{is}.occ = X{is}.occ(1:nocc_max);
                end
            elseif mol.nspin == 4
                F = KSFFT(mol);
                npw = length(X.idxnz);
                nocc_max =  find(X.occ>eps,1,'last');
                psir = cell(2,1);
                psir{1} = F' * X.psi(1:npw,1:nocc_max);
                psir{2} = F' * X.psi(npw+1:end,1:nocc_max);
                realwavefunction.psi = psir;
                realwavefunction.occ = X.occ(1:nocc_max);
            end
            Vexxf = @(Psi)Vexx(Psi, realwavefunction, mol, options);
        case {'ACE','ace'}
            if ~mol.lsda
                [mol.zeta, options] = ACE(X,[],mol,options);
            else
                mol.zeta = cell(2,1);
                [mol.zeta{1}, options] = ACE(X{1},[],mol,options);
                [mol.zeta{2}, options] = ACE(X{2},[],mol,options);
            end
            Vexxf = @(Psi)Vexx(Psi, [], mol, options);
        case {'default','qrcp','kmeans'}
            if ~mol.lsda
                [mol.zeta, options] = ACEISDF(X,[],mol,options);
            else
                mol.zeta = cell(2,1);
                if ~isfield(options.isdfoptions,'ind_mu')
                    ind_mu = cell(2,1);
                    options.isdfoptions.ind_mu = ind_mu;
                end
                [mol.zeta{1}, options] = ACEISDF(X{1},[],mol,options);
                [mol.zeta{2}, options] = ACEISDF(X{2},[],mol,options);
            end
            Vexxf = @(Psi)Vexx(Psi, [], mol, options);
    end
end
fprintf('Time for constructing Fock/ACE operator = %20.4e\n',toc(vexxt));
fprintf('----------------------------------\n');
end
