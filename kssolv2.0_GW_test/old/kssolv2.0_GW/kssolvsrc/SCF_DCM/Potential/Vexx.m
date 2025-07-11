function VexxPsi = Vexx(Y, Phi, F, exxgkk, mol)
% Exact Exchange Operator
    
    n123 = mol.n1 * mol.n2 * mol.n3;
    if isa(Y, 'Wavefun')
        Psi = F' * Y.psi;
    else
        Psi = Y;
    end
    
    VexxPsi = zeros(size(Psi));
    F2 = KSFFT(mol,1);
    for i = 1 : size(Psi,2)
        %VexxPsi(:,i) = sum(Phi.*(F2'*bsxfun(@times,F2*bsxfun(@times,conj(Phi),Psi(:,i)),exxgkk)),2);
        VexxPsi(:,i) = sum(Phi.*(F'*bsxfun(@times,F*bsxfun(@times,conj(Phi),Psi(:,i)),exxgkk)),2);
    end
    VexxPsi = -0.25 * VexxPsi * mol.vol;
    
    if isa(Y, 'Wavefun')
        Y.psi = F * VexxPsi;
        VexxPsi = Y;
    end
end
