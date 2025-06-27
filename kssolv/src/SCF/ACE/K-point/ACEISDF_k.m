function [xi, options] = ACEISDF_k(psi_in,phi_in,mol,options)
% ACEISDF_K Adaptively Compressed Exchange Operator with Interpolative
% Separable Density Fitting algorithm with k-point sampling for hybrid 
% functional calculation
%    [xi, options] = ACEISDF_K(psi_in,phi_in,mol,options) generate xi_k 
%    for inner hybrid functional iteration. psi_in: Any orbitals to do 
%    V_x*psi_in, project V_x into the space spanned by psi_in phi_in: 
%    Occupied orbitals to construct V_x[{phi_in}]. phi_in = [] means it 
%    is same to psi_in
%
%    This function is used for Crystal with k-point sampling in case of
%    spin-restricted and spin-unrestricted calculations

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

psi_phi_same = false;
if isempty(phi_in)
	psi_phi_same = true;
	phi_in = psi_in;
end

n1 = mol.n1;n2 = mol.n2;n3 = mol.n3;
n123  = n1*n2*n3;
ispin = psi_in{1}.ispin;

F = KSFFT(mol);
exxgkk = options.exxgkk;
% neglect the influence of the phi which has very small occupying rate to
% the Fock operator
nk_psi = psi_in.nkpts;
nk_phi = phi_in.nkpts;
wks = phi_in.wks;
nb_psi = size(psi_in{1}.psi,2);
nocc_max = zeros(nk_phi,1);
occ = cell(1,nk_phi);
for ik = 1:nk_phi
    nocc_max(ik) = find(phi_in{ik}.occ>eps,1,'last');
    occ{ik} = phi_in{ik}.occ(1:nocc_max(ik));
end

Psi3 = zeros(n123,nb_psi,nk_psi);
Phi3 = cell(1,nk_phi);
for ik = 1:nk_psi
    Psi3(:,:,ik) = F'*psi_in{ik}.psi;
    if psi_phi_same
        Phi3{ik} = Psi3(:,1:nocc_max(ik),ik);
    end
end
if ~psi_phi_same
    for ik = 1:nk_phi
        Phi3{ik} = F'*phi_in{ik}.psi(:,1:nocc_max(ik));
    end
end
% Reshape Psi and Phi to 2D array to perform ISDF decomposition
Psi2  = reshape(Psi3,n123,nb_psi*nk_psi);
Phi2 = cell2mat(Phi3);
if ~options.fixip
    [zeta_mu,ind_mu,options] = isdf(mol,conj(Phi2),Psi2,[],options,ispin);
else
    [zeta_mu,ind_mu,options] = isdf(mol,conj(Phi2),Psi2,options.ind_mu,options,ispin);
end
    
if ~mol.lsda
    options.ind_mu = ind_mu;
    if ~options.dynamic_samp
       options.fixip = true;
    end
else
    options.ind_mu{ispin} = ind_mu;
    if ~options.dynamic_samp && ~isempty(options.ind_mu{1})...
       && ~isempty(options.ind_mu{2})
       options.fixip = true;
    end
end

if numel(options.rank) == 1
    rank = options.rank;
else
    rank = options.rank(ispin);
end

ng = length(psi_in{1}.idxnz);
% contraction of band index
P_ru  = zeros(n123,rank,nk_phi);
P_vu  = zeros(rank,rank,nk_phi);
for ik = 1:nk_phi
    P_ru(:,:,ik) = Phi3{ik}*(Phi3{ik}(ind_mu,:)'.*(occ{ik}*wks(ik)));
    P_vu(:,:,ik) = P_ru(ind_mu,:,ik);
end

% transform zeta_mu from real space to reciprocal space
% and store it in Fzeta, coefficients are ignored for efficiency
% 
% here is the major difference of ACE-ISDF formulation for Gamma and multi-k-points
% sampling, where the V\zeta becomes V_(k-k')\zeta
% to reduce the time for FFT of zeta, we calculate and save it in Fzeta
% rather than call Poisson solver.
Fzeta = zeros(n123,rank);
rhoc = zeros(n1,n2,n3);
for i = 1:rank
	rhoc = reshape(zeta_mu(:,i),n1,n2,n3);
    rhoc = fftn(rhoc);
    Fzeta(:,i) = rhoc(:);
end

xi = cell(nk_psi,1);
% Calculate M_ijk and W_jk according to the tensors obtained by ISDF
if options.aceconv
        % Construct ACE operator by fourier convolution in reciprocal space.
        % It is faster than normal ACE method, but its code is more complicated.
        assert(psi_phi_same,'Phi and psi must be same when use Fourier convolution!');
        exxidxnz = options.exxidxnz;
        nkxyz = mol.nkxyz;
        nkx=nkxyz(1);nky=nkxyz(2);nkz=nkxyz(3);
        nk_phi2 = prod(2*nkxyz-1);
        % Construction of kgrid for convolution between P_k and K_k
        % index mu is separated by for loop to reduce memory burden
        P2_vu = zeros(rank,2*nkx-1,2*nky-1,2*nkz-1);
        P2_ru = zeros(n123,2*nkx-1,2*nky-1,2*nkz-1);
        VM = zeros(rank,rank,nk_phi);
        VW = zeros(ng,rank,nk_phi);

        for ir=1:rank
            % Solving Poisson equations by FFT
            K2_ru = zeros(n123,nk_phi2);
            K2_vu = zeros(rank,nk_phi2);
            Fzeta_tmp = Fzeta(:,ir);   
            for ik=1:nk_phi2
                vc = Fzeta_tmp.*exxgkk(:,ik);
                vc = reshape(vc,n1,n2,n3);
                vcr = ifftn(vc);
                K2_ru(:,ik) = vcr(:);
            end
            
            K2_vu = zeta_mu'*K2_ru*mol.vol/n123;
            % Put K_k in kgrid for each mu index
            K2_ru = reshape(K2_ru,n123,2*nkx-1,2*nky-1,2*nkz-1);
            K2_vu = reshape(K2_vu,rank,2*nkx-1,2*nky-1,2*nkz-1);
            % Put P_k in kgrid for each mu index
            P2_vu(:) = 0;
            P2_ru(:) = 0;
            P2_vu(:,exxidxnz)=P_vu(:,ir,:);
            P2_ru(:,exxidxnz)=P_ru(:,ir,:);
            % Nr*Nu for loops exist here which may cause very slow speed
            % Calculation of VM_k = ifft(fft(P2_vu)*fft(K2_vu)) by convolution theorem
            K2_vu=fft(K2_vu,[],2); K2_vu=fft(K2_vu,[],3); K2_vu=fft(K2_vu,[],4);
            P2_vu=fft(P2_vu,[],2); P2_vu=fft(P2_vu,[],3); P2_vu=fft(P2_vu,[],4);
            V_vu=K2_vu.*P2_vu;
            V_vu=ifft(V_vu,[],2); V_vu=ifft(V_vu,[],3); V_vu=ifft(V_vu,[],4);
            VM(:,ir,:) = V_vu(:,exxidxnz);
            % Calculation of VW_k = ifft(fft(P2_ru)*fft(P2_ru)) by convolution theorem
            K2_ru=fft(K2_ru,[],2); K2_ru=fft(K2_ru,[],3); K2_ru=fft(K2_ru,[],4);
            P2_ru=fft(P2_ru,[],2); P2_ru=fft(P2_ru,[],3); P2_ru=fft(P2_ru,[],4);
            V_ru=reshape(K2_ru.*P2_ru,n123,nk_phi2);
            % VW is transformed to reciprocal space first
            % The FFT can be performed at last too (the speed of two methods is uncertain)
            V_gu=F*V_ru;
            V_gu=reshape(V_gu,ng,2*nkx-1,2*nky-1,2*nkz-1);
            V_gu=ifft(V_gu,[],2); V_gu=ifft(V_gu,[],3); V_gu=ifft(V_gu,[],4);
            %VW(:,ir,:) = V_gu(:,exxidxnz)*sqrt(nk_phi2);
            VW(:,ir,:) = V_gu(:,exxidxnz);
        end
        for ik=1:nk_phi
            M = Psi3(ind_mu,:,ik)'*VM(:,:,ik)*Psi3(ind_mu,:,ik);
            M = (M+M')/2;
            W = VW(:,:,ik)*Psi3(ind_mu,:,ik);
            R = chol(M);
            xi{ik} = W/R;
        end
else
    % Construct ACE operator by standard convolution.
    % The code is simpler, but the speed is slower than Fourier convolution version.
    % psi_in and phi_in can be totally different.
	if options.store_tensors
		fprintf('Begin the calculation of tensors K_ru and K_vu.\n');
        nkxyz = mol.nkxyz;
        nkx = nkxyz(1);nky = nkxyz(2);nkz = nkxyz(3);
        nk_phi2 = prod(2*nkxyz-1);
        % Calculate and store K_{k-q}_{ru} and K_{k-q}_{vu} in 3D array
        K_ruk = cell(nk_phi2,1);
        K_vuk = cell(nk_phi2,1);
        for ik = 1:nk_phi2
            K_ruk{ik} = zeros(n123,rank);
        end
        for ir = 1:rank
            Fzeta_tmp = Fzeta(:,ir);
            for ik = 1:nk_phi2
                vc = Fzeta_tmp.*exxgkk(:,ik);
                vc = reshape(vc,n1,n2,n3);
                vcr = ifftn(vc);
                K_ruk{ik}(:,ir) = vcr(:);
            end
        end
        for ik = 1:nk_phi2
            K_vuk{ik} = zeta_mu'*K_ruk{ik}*mol.vol/n123;
        end
        [I,J,K] = ndgrid((0:nkx-1)-((0:nkx-1) >= nkx/2)*nkx, ...
        	(0:nky-1)-((0:nky-1) >= nky/2)*nky, ...
            (0:nkz-1)-((0:nkz-1) >= nkz/2)*nkz);
        kpts = [I(:) J(:) K(:)];            
        mapk = zeros(nk_psi,nk_phi);
        for ik1 = 1:nk_psi
            for ik2 = 1:nk_phi
                dk = kpts(ik1,:) - kpts(ik2,:);
                ik = dk + 2*[nkx nky nkz].*(dk < 0) + (dk >= 0);
                mapk(ik1,ik2) = ik(1) + (ik(2)-1)*(2*nkx-1) + (ik(3)-1)*(2*nkx-1)*(2*nky-1);
            end
        end
        fprintf('Finish the calculation of tensors K_ru and K_vu.\n');
	end

    V_vu = zeros(rank);
    V_ru = zeros(n123,rank);
    K_ru = zeros(n123,rank);
    for ik1 = 1:nk_psi
    	V_vu(:) = 0;
        V_ru(:) = 0;
        for ik2 = 1:nk_phi
            if ~options.store_tensors
                % facb equivalent to facb in exx.f90 of QE
                facb = squeeze(exxgkk(:,ik1,ik2));
                for ir = 1:rank
                	vc = Fzeta(:,ir).*facb;
                    vc = reshape(vc,n1,n2,n3);
                    vcr = ifftn(vc);
                    K_ru(:,ir) = vcr(:);
                end
                K_vu = zeta_mu'*K_ru*mol.vol/n123;
            else
                idk = mapk(ik1,ik2);
                K_ru = K_ruk{idk};
                K_vu = K_vuk{idk};
            end
            V_vu = V_vu + P_vu(:,:,ik2).*K_vu;
            V_ru = V_ru + P_ru(:,:,ik2).*K_ru;
        end
        M = Psi3(ind_mu,:,ik1)'*V_vu*Psi3(ind_mu,:,ik1);
        M = (M+M')/2;
        W = F*(V_ru*Psi3(ind_mu,:,ik1));
        R = chol(M);
        xi{ik1} = W/R;
    end
end

end

    

















