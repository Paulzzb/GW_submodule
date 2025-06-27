function [xi, options] = ACEISDF(psi_in,phi_in,mol,options)
% ACEISDF Adaptively Compressed Exchange Operator with Interpolative
% Separable Density Fitting for hybrid functional.
%    [xi, options] = ACEISDF(psi_in,phi_in,mol,options) calculates ACE 
%    operator with the help of ISDF.
%
%    psi_in: Any orbitals to do V_x*psi_in, project V_x into the space 
%            spanned by psi_in.
%    phi_in: Occupied orbitals to construct V_x[{phi_in}]. phi_in = [] 
%            means it is same to psi_in.
%
%    This function is used for Molecule with \Gamma point sampling.

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
ispin = psi_in.ispin;

F = KSFFT(mol);
exxgkk = options.exxgkk;

if ~mol.noncolin
    nocc_max = find(phi_in.occ>eps,1,'last');
    occ = phi_in.occ(1:nocc_max);
        
    Psi = F'*psi_in.psi;     
    if psi_phi_same          
        Phi = Psi(:,1:nocc_max);
    else
        Phi = F'*phi_in(:,1:nocc_max);
    end
        
    if ~options.fixip
        [zeta_mu,ind_mu,options] = isdf(mol,conj(Phi),Psi,[],options,ispin);
    else
        [zeta_mu,ind_mu,options] = isdf(mol,conj(Phi),Psi,options.ind_mu,options,ispin);
    end

    % Fix interpolation points in remaining SCF iterations if dynamic sampling is not used.    
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

    kernel = Poisson_solver(mol,zeta_mu,exxgkk);
    
    % Equal to (F*zeta_mu)'*(F*kernel)       
    K_vu = zeta_mu'*kernel*mol.vol/n123;
   
    P_vu = Phi(ind_mu,:)*(Phi(ind_mu,:)'.*occ);
    M = Psi(ind_mu,:)'*(K_vu.*P_vu)*Psi(ind_mu,:);
    M = (M+M')/2;
    P_ru = Phi*(Phi(ind_mu,:)'.*occ);

    % Here we perform matrix multiplication first and then FFT
    % to save the time for FFT
    FK = (kernel.*P_ru)*Psi(ind_mu,:);
    W = F*FK;
    R = chol(M);        
    xi = W/R;
else
    npw = length(psi_in.idxnz);
    nbnd = mol.nbnd;
    nocc_max = find(phi_in.occ>eps,1,'last');
    occ = phi_in.occ(1:nocc_max);

    Psi = cell(2,1);
    Psi{1} = F'*psi_in.psi(1:npw,:);
    Psi{2} = F'*psi_in.psi(npw+1:end,:);
    Phi = cell(2,1);
    if psi_phi_same
        Phi{1} = Psi{1}(:,1:nocc_max);
        Phi{2} = Psi{2}(:,1:nocc_max);
    else
        Phi{1} = F'*phi_in(1:npw,1:nocc_max);
        Phi{2} = F'*phi_in(npw+1:end,1:nocc_max);
    end

    if options.numisdf == 1
        if ~options.fixip
            [zeta_mu,ind_mu,options] = isdf(mol,Phi,Psi,[],options,ispin);
        else
            [zeta_mu,ind_mu,options] = isdf(mol,Phi,Psi,options.ind_mu,options,ispin);
        end
        kernel = Poisson_solver(mol,zeta_mu,exxgkk);
        K_uv = zeta_mu'*kernel*mol.vol/n123;
    else
        zeta_mu = cell(2,1);
        ind_mu = cell(2,1);
        kernel = cell(2,1);
        fPhi = cell(2,1);
        for i = 1:2
            if ~options.fixip
                [zeta_mu{i},ind_mu{i},options] = isdf(mol,conj(Phi{i}),Psi{i}, [], options,i);
            else
                [zeta_mu{i},ind_mu{i},options] = isdf(mol,conj(Phi{i}),Psi{i}, options.ind_mu{i}, options,i);
            end
            kernel{i} = Poisson_solver(mol,zeta_mu{i},exxgkk);
            fPhi{i} = Phi{i}(ind_mu{i},:)'.*occ;
        end
    end

    if ~options.dynamic_samp
       options.fixip = true;
    end
    options.ind_mu = ind_mu;

    M = zeros(nbnd,nbnd);
    if options.numisdf == 1
        for i = 1:2
            P_uv = Phi{i}(ind_mu,:)*(Phi{i}(ind_mu,:)'.*occ);
            M = M + Psi{i}(ind_mu,:)'*(P_uv.*K_uv)*Psi{i}(ind_mu,:);
        end
        P_uv = Phi{1}(ind_mu,:)*(Phi{2}(ind_mu,:)'.*occ);
        M_updw = Psi{1}(ind_mu,:)'*(P_uv.*K_uv)*Psi{2}(ind_mu,:);
        M = M + M_updw + M_updw';
        M = (M+M')/2;
    else
        for i = 1:2
            K_uv = zeta_mu{i}'*kernel{i}*mol.vol/n123;
            P_uv = Phi{i}(ind_mu{i},:)*fPhi{i};
            M = M + Psi{i}(ind_mu{i},:)'*(P_uv.*K_uv)*Psi{i}(ind_mu{i},:);
        end
        K_uv = zeta_mu{1}'*kernel{2}*mol.vol/n123;
        clear zeta_mu;
        P_uv = Phi{1}(ind_mu{1},:)*fPhi{2};
        M_updw = Psi{1}(ind_mu{1},:)'*(P_uv.*K_uv)*Psi{2}(ind_mu{2},:);
        M = M + M_updw + M_updw';
        M = (M+M')/2;
    end

    W = cell(2,1);
    if options.numisdf == 1
        for is1 = 1:2
            W{is1} = zeros(n123,nbnd);
            for is2 = 1:2
                P_ru = Phi{is1}*(Phi{is2}(ind_mu,:)'.*occ);
                FK=zeros(size(kernel));
                FK = (kernel.*P_ru);
                W{is1} = W{is1} + FK*Psi{is2}(ind_mu,:);
            end
            W{is1} = F*W{is1};
        end

    else
        for is1 = 1:2
            W{is1} = zeros(n123,nbnd);
            for is2 = 1:2
                P_ru = Phi{is1}*fPhi{is2};
                FK = zeros(size(kernel{is2}));
                FK = (kernel{is2}.*P_ru);
                W{is1} = W{is1} + FK*Psi{is2}(ind_mu{is2},:);
            end
            W{is1} = F*W{is1};
        end
    end
    
    R = chol(M);
    xi = [W{1};W{2}]/R; 
end

end
