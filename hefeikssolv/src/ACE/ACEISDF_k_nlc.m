function [xi, options] = ACEISDF_k_nlc(psi_in,phi_in,mol,options)
% ACE-ISDF method to generate xi_k for inner hybrid functional iteration.
%
% input: psi_in: Any orbitals to do V_x*psi_in, project V_x into the space spanned by psi_in
% phi_in: Occupied orbitals to construct V_x[{phi_in}]. phi_in = [] means it is same to psi_in
%
% For Crystal with k point sampling
% Spin-noncollinear case

psi_phi_same = false;
if isempty(phi_in)
        psi_phi_same = true;
        phi_in = psi_in;
end

n1=mol.n1;n2=mol.n2;n3=mol.n3;
n123  = n1*n2*n3;
ispin = psi_in{1}.ispin;

F = KSFFT(mol);
F2 = KSFFT(mol,options.exxcut);
exxgkk = options.exxgkk;
opt = options.isdfoptions;
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

ng = length(psi_in{1}.idxnz);
Psi3 = {zeros(n123,nb_psi,nk_psi);zeros(n123,nb_psi,nk_psi)};
Phi3 = {cell(1,nk_phi);cell(1,nk_phi)};
for is = 1:2
    idr = (1:ng) + (is -1)*ng;
    for ik = 1:nk_psi
        Psi3{is}(:,:,ik) = F'*psi_in{ik}.psi(idr,:);
        if psi_phi_same
            Phi3{is}{ik} = Psi3{is}(:,1:nocc_max(ik),ik);
        else
            Phi3{is}{ik} = F'*phi_in{ik}.psi(idr,1:nocc_max(ik));
        end
    end
end
% Reshape Psi and Phi to 2D array to perform ISDF decomposition
Psi2  = {reshape(Psi3{1},n123,nb_psi*nk_psi);reshape(Psi3{2},n123,nb_psi*nk_psi)};
Phi2  = {cell2mat(Phi3{1});cell2mat(Phi3{2})};
% transform zeta_mu from real space to reciprocal space
% and store it in Fzeta, coefficients of FFT are ignored for efficiency
if opt.numisdf == 1
    if ~opt.fixip
        [zeta_mu,ind_mu,options]=isdf(Phi2,Psi2,[],options,ispin);
    else
        [zeta_mu,ind_mu,options]=isdf(Phi2,Psi2,opt.ind_mu,options,ispin);
    end 
    rank = options.isdfoptions.rank;
    Fzeta = zeros(n123,rank);
    rhoc = zeros(n1,n2,n3);
    for i = 1:rank
        rhoc = reshape(zeta_mu(:,i),n1,n2,n3);
        rhoc = fftn(rhoc);
        Fzeta(:,i) = rhoc(:);
    end  
else
    zeta_mu = cell(2,1);
    ind_mu = cell(2,1);
    Fzeta = cell(2,1);
    for is = 1:2
        if ~opt.fixip
            [zeta_mu{is},ind_mu{is},options]=isdf(conj(Phi2{is}),Psi2{is},[],options,is);
        else
            [zeta_mu{is},ind_mu{is},options]=isdf(conj(Phi2{is}),Psi2{is},opt.ind_mu,options,is);
        end
        rank = options.isdfoptions.rank;
        Fzeta{is} = zeros(n123,rank(is));
        rhoc = zeros(n1,n2,n3);
        for i = 1:rank(is)
            rhoc = reshape(zeta_mu{is}(:,i),n1,n2,n3);
            rhoc = fftn(rhoc);
            Fzeta{is}(:,i) = rhoc(:);
        end
    end    
end

options.isdfoptions.ind_mu = ind_mu;
if ~opt.dynamic_sample
   options.isdfoptions.fixip = true;
end

acet = tic;
if options.aceconv
    assert(psi_phi_same,'Phi and psi must be same when use Fourier convolution!');
    fprintf('Calculate K_ru\n');
    nkxyz = mol.nkxyz;
    nkx=nkxyz(1);nky=nkxyz(2);nkz=nkxyz(3);
    nk_phi2 = prod(2*nkxyz-1);
    % Calculate and store K_{k-q}_{ru} in cell array
    if opt.numisdf == 1
        K_ruk = cell(rank,1);
        for ir = 1:rank
            K_ruk{ir} = zeros(n123,nk_phi2);
        end
        for ir = 1:rank
            Fzeta_tmp = Fzeta(:,ir);
            for ik = 1:nk_phi2
                vc = Fzeta_tmp.*exxgkk(:,ik);
                vc = reshape(vc,n1,n2,n3);
                vc = ifftn(vc);
                K_ruk{ir}(:,ik) = vc(:);
            end
        end
    else
        K_ruk = {cell(rank(1),1);cell(rank(2),1)};
        for is = 1:2
            for ir = 1:rank(is)
                K_ruk{is}{ir} = zeros(n123,nk_phi2);
            end
            for ir = 1:rank(is)
                Fzeta_tmp = Fzeta{is}(:,ir);
                for ik = 1:nk_phi2
                    vc = Fzeta_tmp.*exxgkk(:,ik);
                    vc = reshape(vc,n1,n2,n3);
                    vc = ifftn(vc);
                    K_ruk{is}{ir}(:,ik) = vc(:);
                end
            end
        end
    end
    fprintf('Finish calculation of K_ru\n');
end

if options.store_tensors
    assert(psi_phi_same,'Phi and psi must be same when use Fourier convolution!');
    fprintf('Calculate K_ru and K_uu\n');                
    nkxyz = mol.nkxyz;
    nkx = nkxyz(1);nky = nkxyz(2);nkz = nkxyz(3);
    nk_phi2 = prod(2*nkxyz-1);
    % Calculate and store K_{k-q}_{ru} and K_{k-q}_{uu} in 3D array
    if opt.numisdf == 1
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
                vc = ifftn(vc);
                K_ruk{ik}(:,ir) = vc(:);
            end
        end
        fprintf('Finish calculation of K_ru\n');
        for ik = 1:nk_phi2
            K_vuk{ik} = zeta_mu'*K_ruk{ik};
        end
        fprintf('Finish calculation of K_vu\n');
    else
        K_ruk = {cell(nk_phi2,1);cell(nk_phi2,1)};
        K_vuk = {cell(nk_phi2,1);cell(nk_phi2,1)};
        K_vuk_updw = cell(nk_phi2,1);
        for is = 1:2
            for ik = 1:nk_phi2
                K_ruk{is}{ik} = zeros(n123,rank(is));
            end
            for ir = 1:rank(is)
                Fzeta_tmp = Fzeta{is}(:,ir);
                for ik = 1:nk_phi2
                    vc = Fzeta_tmp.*exxgkk(:,ik);
                    vc = reshape(vc,n1,n2,n3);
                    vc = ifftn(vc);
                    K_ruk{is}{ik}(:,ir) = vc(:);
                end
            end
        end
        fprintf('Finish calculation of K_ru\n');
        for is = 1:2
            for ik = 1:nk_phi2
                K_vuk{is}{ik} = zeta_mu{is}'*K_ruk{is}{ik};
            end
        end
        for ik = 1:nk_phi2
            K_vuk_updw{ik} = zeta_mu{1}'*K_ruk{2}{ik};
        end
        fprintf('Finish calculation of K_vu\n');
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
    fprintf('Finish calculation of K_ru and K_vu\n');
end

xi = cell(nk_psi,1);
% Calculate M_ijk and W_jk according to the tensors obtained by ISDF
% M_updw = \dagger{M_dwup} can be used to reduce computational time
M = cell(nk_psi,1);
W = {cell(nk_psi,1) cell(nk_psi,1)};
for ik = 1:nk_psi
    M{ik} = zeros(nb_psi,nb_psi);
    W{1}{ik} = zeros(ng,nb_psi);
    W{2}{ik} = zeros(ng,nb_psi);
end 
% transform scalar rank to array to avoid too many judgement 
% about opt.numisdf == 1
if opt.numisdf == 1
    rank = [rank rank];
    ind_mu = {ind_mu;ind_mu};
end
for is1 = 1:2
    for is2 = 1:2
        calM = ~(is1 == 2 && is2 == 1);
        % contraction of band index
        P_ru  = zeros(n123,rank(is2),nk_phi);
        for ik=1:nk_phi
            P_ru(:,:,ik)=Phi3{is1}{ik}*(Phi3{is2}{ik}(ind_mu{is2},:)'.*(occ{ik}*wks(ik)));
        end
        if calM
            P_vu  = zeros(rank(is1),rank(is2),nk_phi);
            for ik=1:nk_phi
                P_vu(:,:,ik)=Phi3{is1}{ik}(ind_mu{is1},:)*(Phi3{is2}{ik}(ind_mu{is2},:)'.*(occ{ik}*wks(ik)));
            end
        end

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
            P2_ru = zeros(n123,2*nkx-1,2*nky-1,2*nkz-1);
            VW = zeros(ng,rank(is2),nk_phi);
            if calM
                P2_vu = zeros(rank(is1),2*nkx-1,2*nky-1,2*nkz-1);
                VM = zeros(rank(is1),rank(is2),nk_phi);
            end

            for ir = 1:rank(is2)
                % The Poisson equations have been solved beforehand
                if opt.numisdf == 1
                    K2_ru = K_ruk{ir};
                    if calM
                        K2_vu = zeta_mu'*K2_ru;
                    end
                else
                    K2_ru = K_ruk{is2}{ir};
                    if calM
                        K2_vu = zeta_mu{is1}'*K2_ru;
                    end
                end
                % Put K_k and P_k in kgrid for each mu index
                K2_ru=reshape(K2_ru,n123,2*nkx-1,2*nky-1,2*nkz-1);
                P2_ru(:) = 0;
                P2_ru(:,exxidxnz)=P_ru(:,ir,:);
                if calM
                    K2_vu=reshape(K2_vu,rank(is1),2*nkx-1,2*nky-1,2*nkz-1);
                    P2_vu(:) = 0;
                    P2_vu(:,exxidxnz)=P_vu(:,ir,:);
                end
                % Fourier convolution for calculating M_k
                if calM
                    K2_vu=ifft(K2_vu,[],2); K2_vu=ifft(K2_vu,[],3); K2_vu=ifft(K2_vu,[],4);
                    P2_vu=ifft(P2_vu,[],2); P2_vu=ifft(P2_vu,[],3); P2_vu=ifft(P2_vu,[],4);
                    V_vu=K2_vu.*P2_vu;
                    V_vu=fft(V_vu,[],2); V_vu=fft(V_vu,[],3); V_vu=fft(V_vu,[],4);
                    VM(:,ir,:) = V_vu(:,exxidxnz);
                end
                % Fourier convolution for calculating W_k
                K2_ru=ifft(K2_ru,[],2); K2_ru=ifft(K2_ru,[],3); K2_ru=ifft(K2_ru,[],4);
                P2_ru=ifft(P2_ru,[],2); P2_ru=ifft(P2_ru,[],3); P2_ru=ifft(P2_ru,[],4);
                V_ru=reshape(K2_ru.*P2_ru,n123,nk_phi2);
                V_gu=F*V_ru;
                V_gu=reshape(V_gu,ng,2*nkx-1,2*nky-1,2*nkz-1);
                V_gu=fft(V_gu,[],2); V_gu=fft(V_gu,[],3); V_gu=fft(V_gu,[],4);
                VW(:,ir,:) = V_gu(:,exxidxnz)*sqrt(nk_phi2);
            end
            for ik = 1:nk_psi
                if calM
                    M_tmp = Psi3{is1}(ind_mu{is1},:,ik)'*VM(:,:,ik)*Psi3{is2}(ind_mu{is2},:,ik);
                    if is1 == 1 && is2 == 2 
                        M{ik} = M{ik} + M_tmp + M_tmp';
                    else
                        M{ik} = M{ik} + M_tmp;
                    end
                end
                W{is1}{ik} = W{is1}{ik} + VW(:,:,ik)*Psi3{is2}(ind_mu{is2},:,ik);
            end
        else
            % Construct ACE operator by standard convolution.
            % The code is simpler, but the speed is slower than Fourier convolution version.
            % psi_in and phi_in can be totally different.
            V_ru = zeros(n123,rank(is2));
            K_ru = zeros(n123,rank(is2));
            if calM
                V_vu = zeros(rank(is1),rank(is2));
            end
            for ik1=1:nk_psi
                V_ru(:)=0;
                if calM
                    V_vu(:)=0;
                end
                for ik2=1:nk_phi
                    if ~options.store_tensors
                        %facb equivalent to facb in exx.f90 of QE
                        facb = squeeze(exxgkk(:,ik1,ik2));
                        for ir = 1:rank(is2)
                            vc = Fzeta{is2}(:,ir).*facb;
                            vc = reshape(vc,n1,n2,n3);
                            vc = ifftn(vc);
                            K_ru(:,ir) = vc(:);
                        end
                        if calM
                            K_vu=zeta_mu{is1}'*K_ru;
                        end
                    else
                        idk = mapk(ik1,ik2);
                        if opt.numisdf == 1
                            K_ru = K_ruk{idk};
                        else
                            K_ru = K_ruk{is2}{idk};
                        end
                        if calM
                            if opt.numisdf == 1
                                K_vu = K_vuk{idk};
                            else
                                if is1 == 1 && is2 == 2
                                    K_vu = K_vuk_updw{idk};
                                else
                                    K_vu = K_vuk{is1}{idk};
                                end
                            end
                        end
                    end
                    V_vu=V_vu + P_vu(:,:,ik2).*K_vu;
                    V_ru=V_ru + P_ru(:,:,ik2).*K_ru;
                end

                if calM
                    M_tmp = Psi3{is1}(ind_mu{is1},:,ik1)'*V_vu*Psi3{is2}(ind_mu{is2},:,ik1);
                    if is1 == 1 && is2 == 2
                        M{ik1} = M{ik1} + M_tmp + M_tmp';
                    else
                        M{ik1} = M{ik1} + M_tmp;
                    end
                end
                W{is1}{ik1} = W{is1}{ik1} + (F*V_ru)*Psi3{is2}(ind_mu{is2},:,ik1);
                %W{is1}{ik1} = W{is1}{ik1} + F*(V_ru*Psi3{is2}(ind_mu{is2},:,ik1);
            end
        end
    end
end

for ik = 1:nk_psi
    M{ik} = (M{ik}+M{ik}')/2;
    R = chol(M{ik});
    xi{ik} = [W{1}{ik};W{2}{ik}]/R*(sqrt(n123)/mol.vol);;
end
fprintf('Calculation time of W_k and M_k is %20.4e\n',toc(acet));
end




            





























