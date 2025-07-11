function [mol,H,X,info] = ksload(savefile)
% KSLOAD load KSSOLV related data into workspace.
%
%   [mol,H,X,info] = KSLOAD(savefile) loads molecule, hamiltonian,
%   wavefunction and scf iteration information from a file named as
%   savefile.
%
%   See also kssave.

ksdat = load(savefile);
mol = ksdat.mol;
H = ksdat.H;
X = ksdat.X;
info = ksdat.info;

end
