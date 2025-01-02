function mol = set_nelup_neldw(mol)
% Set the number of up and down electrons  
tot_mag = mol.tot_mag;
nel = mol.nel;

if tot_mag < 0
 mol.nelup = ceil(mol.nel/2);
 mol.neldw = mol.nel - mol.nelup;
else
 mol.nelup = (nel + tot_mag)/2;
 mol.neldw = (nel - tot_mag)/2;
end
 
end

