function [xi, options] = ACE_k(psi_in,phi_in,mol,options)
% ACE without ISDF, the main calculating time is used in the generation of
% matrix M_ijk = <psi_ik|Vexx|psi_jk> and w_jk = Vexx|psi_jk>
% For Crystal with k-point sampling

psi_phi_same = false;
if isempty(phi_in)
	psi_phi_same = true;
	phi_in = psi_in;
end

n1=mol.n1;n2=mol.n2;n3=mol.n3;
n123  = n1*n2*n3;
nspin = mol.nspin;

F = KSFFT(mol);
F2 = KSFFT(mol,options.exxcut);
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

if ~mol.noncolin
    Psi3 = zeros(n123,nb_psi,nk_psi);
    Phi3 = cell(1,nk_phi);
    for ik=1:nk_psi
        Psi3(:,:,ik)=F'*psi_in{ik}.psi;
        if psi_phi_same
            Phi3{ik} = Psi3(:,1:nocc_max(ik),ik);
        end
    end
    if ~psi_phi_same
        for ik=1:nk_phi
            Phi3{ik} = F'*phi_in{ik}.psi(:,1:nocc_max(ik));
        end
    end        
    % calculate chi_k by a for loop of index k
    xi = cell(nk_psi,1);
    for k = 1 : nk_psi
        fprintf('Construct ACE operator for %d-th k point\n',k);
        VexxPsi = zeros(n123,nb_psi);
        for i = 1:nb_psi
            for q = 1:nk_phi
                facb = squeeze(exxgkk(:,k,q));
                prod_pair = conj(Phi3{q}).*Psi3(:,i,k);
                kernel = Poisson_solver(mol,prod_pair,facb);
                VexxPsi(:,i) = VexxPsi(:,i) + (Phi3{q}.*kernel)*occ{q}*wks(q);
            end
        end
        
        VexxPsig = F*VexxPsi;
        M = Psi3(:,:,k)'*VexxPsi;
        M=(M+M')/2;
        R = chol(M); 
        xi{k} = VexxPsig/R;
    end
else
    Psi3_up = zeros(n123,nb_psi,nk_psi);
    Psi3_dw = zeros(n123,nb_psi,nk_psi);
    Phi3_up = cell(1,nk_phi);
    Phi3_dw = cell(1,nk_phi);
    npw = length(psi_in{1}.idxnz);
    for ik=1:nk_psi
        Psi3_up(:,:,ik)=F'*psi_in{ik}.psi(1:npw,:);
        Psi3_dw(:,:,ik)=F'*psi_in{ik}.psi(npw+1:end,:);
        if psi_phi_same
            Phi3_up{ik} = Psi3_up(:,1:nocc_max(ik),ik);
            Phi3_dw{ik} = Psi3_dw(:,1:nocc_max(ik),ik);
        end
    end
    if ~psi_phi_same
        for ik=1:nk_phi
            Phi3_up{ik} = F'*phi_in{ik}.psi(1:npw,1:nocc_max(ik));
            Phi3_dw{ik} = F'*phi_in{ik}.psi(npw+1:end,1:nocc_max(ik));
        end
    end
    xi = cell(nk_psi,1);
    for k = 1 : nk_psi
        fprintf('Construct ACE operator for %d-th k point\n',k);
        VexxPsi = zeros(n123*2,nb_psi);
        for i = 1: nb_psi    
            for q = 1 : nk_phi   
                facb = squeeze(exxgkk(:,k,q));
                prod_pair = conj(Phi3_up{q}).*Psi3_up(:,i,k)+conj(Phi3_dw{q}).*Psi3_dw(:,i,k);       
                kernel = Poisson_solver(mol,prod_pair,facb);
                VexxPsi(1:n123,i) = VexxPsi(1:n123,i) + (Phi3_up{q}.*kernel)*occ{q}*wks(q);
                VexxPsi(n123+1:end,i) = VexxPsi(n123+1:end,i) + (Phi3_dw{q}.*kernel)*occ{q}*wks(q);
            end
        end
        VexxPsig = zeros(npw*2,nb_psi);
        VexxPsig(1:npw,:) = F*VexxPsi(1:n123,:);
        VexxPsig(npw+1:end,:) = F*VexxPsi(n123+1:end,:);
        M=[Psi3_up(:,:,k);Psi3_dw(:,:,k)]'*VexxPsi;
        M=(M+M')/2;
        R = chol(M); 
        xi{k} = VexxPsig/R;
    end
end   

for ik = 1:nk_psi
    xi{ik} = xi{ik}*(sqrt(n123)/mol.vol);
end

end
