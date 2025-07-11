function Exx = getExx(X, Vexxf, mol)
    n123 = mol.n1 * mol.n2 * mol.n3;
    if ~isempty(X.occ), X.psi = X.psi(:,logical(X.occ));   end
    VexxPsi = Vexxf(X);
    Exx = real(sum(sum(conj(X).*VexxPsi)));
end
