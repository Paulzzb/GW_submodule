function nel = getnel(mol,rho)
%Calculate the charge by denisty
if ~iscell(rho)
    nel = sum(rho,'all')*mol.vol/mol.n1/mol.n2/mol.n3;
else
    nel = (sum(rho{1},'all')+sum(rho{2},'all'))*mol.vol/mol.n1/mol.n2/mol.n3;
end

