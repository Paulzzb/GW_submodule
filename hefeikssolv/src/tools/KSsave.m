function filename = KSsave(mol,H,X,info)

filename = [mol.name '.mat'];
save(filename,'mol','H','X','info');

end