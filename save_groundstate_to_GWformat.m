function save_groundstate_to_GWformat(mol, H, X0, info, dir)
% 1. Extract groundstate information from some KSSOLV class variables
%    to standard format variables.
% 2. Concrete them as a struct, save it into fileName  
ha2ry = 2.0;

sys = initsys(mol);
rhor = H.rho;
vxc = getVhxc(mol, rhor);
ev = info.Eigvals;
occupation = X0.occ;
psig = X0.psi;
Vxc = [];
reciprocal_grid_info = [];


% Calculate Vxc(n)
nb = length(ev); nr = length(vxc(:));
F = KSFFT(mol);
psir = F' * X0.psi;
Vxc = zeros(nb, 1);
for it=1:nb
  Vxc(it) = sumel(vxc(:) .* ((abs(psir(:,it))).^2)) * (sys.vol)^2 / nr;
end

% Prepare reciprocal grid information
ggrid = Ggrid(mol);
C = mol.supercell;
xyz = [ggrid.gkx, ggrid.gky, ggrid.gkz] * C' * 2 * pi;
% KSSOLV use Hartree unit, and our code use Rydberg unit, so we need to convert it.
reciprocal_grid_info = struct('xyz', xyz, 'idxnz', ggrid.idxnz, 'wfncut', ha2ry*ggrid.ecut);
  
% Transform
% Prepare groundstate struct and save
groundstate.rhor = rhor;
groundstate.Vxc = Vxc;
groundstate.ev = ev;
groundstate.psig = psig;
groundstate.sys = sys;
groundstate.occupation = occupation;
groundstate.reciprocal_grid_info = reciprocal_grid_info;

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
  sys.qk = [0, 0, 0];
end % EOF
