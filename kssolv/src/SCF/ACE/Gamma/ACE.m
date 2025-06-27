function [xi, options] = ACE(psi_in,phi_in,mol,options)
% ACE Adaptively Compressed Exchange Operator for hybrid functional.
%    [xi, options] = ACE(psi_in,phi_in,mol,options) calculates ACE operator
%    without ISDF, the main calculating time is used in the generation of
%    matrix M_ij = <psi_i|Vexx|psi_j> and w_j = Vexx|psi_j>
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

F = KSFFT(mol);
exxgkk = options.exxgkk;

% neglect the influence of the phi which has very small occupying rate to
% the Fock operator
nocc_max = find(phi_in.occ>eps,1,'last');
occ = phi_in.occ(1:nocc_max);

if ~mol.noncolin
    Psi = F'*psi_in.psi;     
    if psi_phi_same          
        Phi = Psi(:,1:nocc_max);
    else
        Phi = F'*phi_in(:,1:nocc_max);
    end
    
    VexxPsi = zeros(size(Psi));
    for i = 1:size(Psi,2)        
        prod_pair = conj(Phi).*Psi(:,i);
        kernel = Poisson_solver(mol,prod_pair,exxgkk);         
        VexxPsi(:,i) = (Phi.*kernel)*occ;
    end
    VexxPsig = F*VexxPsi;
    %M = psi_in.psi'*VexxPsig;
   
    xi = zeros(size(VexxPsig));
    ndivid = 1;
    np = floor(mol.nbnd/ndivid);
    for i = 1:ndivid-1
        idx = (1:np) + (i-1)*np;
        M = psi_in.psi(:,idx)'*VexxPsig(:,idx);
        M = (M+M')/2;
        R = chol(M); 
        xi(:,idx) = VexxPsig(:,idx)/R;
    end 
    
    idx = (np*(ndivid-1)+1):mol.nbnd;
    M = psi_in.psi(:,idx)'*VexxPsig(:,idx);
    M = (M+M')/2;
    R = chol(M); 
    xi(:,idx) = VexxPsig(:,idx)/R; 
else
    npw = length(psi_in.idxnz);
    Psi_up = F'*psi_in.psi(1:npw,:); 
    Psi_dw = F'*psi_in.psi(npw+1:end,:);
    if psi_phi_same
        Phi_up = Psi_up(:,1:nocc_max);
        Phi_dw = Psi_dw(:,1:nocc_max);
    else
        Phi_up = F'*phi_in(1:npw,1:nocc_max);
        Phi_dw = F'*phi_in(npw+1:end,1:nocc_max);
    end
    VexxPsi = zeros(n123*2,size(Psi_up,2));
    
    for i = 1 : size(Psi_up,2)
        prod_pair = conj(Phi_up).*Psi_up(:,i) + conj(Phi_dw).*Psi_dw(:,i);
        kernel = Poisson_solver(mol,prod_pair,exxgkk);
        VexxPsi(1:n123,i) = (Phi_up.*kernel)*occ;
        VexxPsi(n123+1:end,i) = (Phi_dw.*kernel)*occ;
    end
    VexxPsig = zeros(size(psi_in.psi));
    VexxPsig(1:npw,:) = F*VexxPsi(1:n123,:);
    VexxPsig(npw+1:end,:) = F*VexxPsi(n123+1:end,:);
    M = psi_in.psi'*VexxPsig;
    M = (M+M')/2;
    R = chol(M); 
    xi = VexxPsig/R;
end

% M = (M+M')/2;
% R = chol(M); 
% xi = VexxPsig/R;

end
