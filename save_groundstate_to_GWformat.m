function save_groundstate_to_GWformat(mol, H, X0, info, dir)
% 1. Extract groundstate information from some KSSOLV class variables
%    to standard format variables.
% 2. Concrete them as a struct, save it into fileName  

sys = initsys(mol);
rhor = H.rho;
vxc = getVhxc(mol, rhor);
ev = info.Eigvals;
occupation = X0.occ;
psig = X0.psi;

% Calculate Vxc(n)
nb = length(ev); nr = length(vxc(:));
F = KSFFT(mol);
psir = F' * X0.psi;
Vxc = zeros(nb, 1);
for it=1:nb
  Vxc(it) = sumel(vxc(:) .* ((abs(psir(:,it))).^2)) * (sys.vol)^2 / nr;
end
groundstate.rhor = rhor;
groundstate.Vxc = Vxc;
groundstate.eigval = ev;
groundstate.psi = psig;
groundstate.sys = sys;
groundstate.occupation = occupation;

fileName = fullfile(dir, 'groundstate.mat');
save(fileName, 'groundstate');

end


function sys = initsys(mol)
  F = KSFFT(mol);
  [sys.ng, sys.nr] = size(F);
  clear F
  sys.vol = mol.vol;
  sys.xyzlist = mol.xyzlist;
  sys.n1 = mol.n1;
  sys.n2 = mol.n2;
  sys.n3 = mol.n3;
  sys.ne = mol.nel;
  sys.supercell = mol.supercell;
end % EOF