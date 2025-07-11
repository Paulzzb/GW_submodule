function savefile = kssave(mol,H,X,info)
% KSSAVE save KSSOLV related data into file.
%
%   savefile = KSSAVE(mol,H,X,info) saves molecule, hamiltonian,
%   wavefunction and scf iteration information into a file named as
%   'mol.name'.mat and returns the name of saved file.
%
%   See also ksload.

savefile = [mol.name '.mat'];
save(savefile,'mol','H','X','info');

end
