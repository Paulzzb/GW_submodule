function rhout = rho_mix(mol,A,rho1,rho2)
% rhout = rho1 + A*rho2
if mol.nspin == 1 
    rhout = rho1 + A*rho2;
else
    rhout = cell(mol.nspin,1);
    for i = 1:mol.nspin
        rhout{i} = rho1{i} + A*rho2{i};
    end
end
        
