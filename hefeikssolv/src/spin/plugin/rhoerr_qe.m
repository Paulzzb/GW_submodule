function rhoerr = rhoerr_qe(mol,rho,rhoin)
% Convergence criterion for QE's scf iteration. 
% rho should be transformed first.
if mol.nspin == 1
    rho = transform_rho(mol,rho,'real2rec');
    rhoin = transform_rho(mol,rhoin,'real2rec');
    rho = rho - rhoin;
    rhoerr = rho_ddot(mol,rho,rho);
elseif mol.nspin == 2
    rho = transform_rho(mol,rho,'updw2sum');
    rhoin = transform_rho(mol,rhoin,'updw2sum');
    rho = transform_rho(mol,rho,'real2rec');
    rhoin = transform_rho(mol,rhoin,'real2rec');
    rho{1} = rho{1} - rhoin{1};
    rho{2} = rho{2} - rhoin{2};
    rhoerr = rho_ddot(mol,rho,rho);
elseif mol.nspin == 4
    rho = transform_rho(mol,rho,'real2rec');
    rhoin = transform_rho(mol,rhoin,'real2rec');
    for i = 1:4
        rho{i} = rho{i} - rhoin{i};
    end
    rhoerr = rho_ddot(mol,rho,rho);
end
