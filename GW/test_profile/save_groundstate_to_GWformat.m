function save_groundstate_to_GWformat(mol, H, X0, info, dir)
% 1. Extract groundstate information from some KSSOLV class variables
%    to standard format variables.
% 2. Concrete them as a struct, save it into fileName  
ha2ry = 2.0;

nkpts = 1; nspin = 1; nspinor = 1;
%% Change the following part
kpts = [0,0,0];
groundstate.nkpts = nkpts;
groundstate.kpts = kpts;
groundstate.nspin = nspin;
groundstate.nspinor = nspinor;
fprintf('save_groundstate_to_GWformat: nkpts = %d, nspin = %d\n, nspinor = %d.\n', nkpts, nspin, nspinor);
fprintf('Current set nkpts, nspin, and nspinor always 1\n');
%%

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
xyz = cell(nkpts, 1);
for ik = 1:nkpts
  xyz{ik} = round([ggrid.gkx, ggrid.gky, ggrid.gkz] * C' / 2 / pi);
end
% KSSOLV use Hartree unit, and our code use Rydberg unit, so we need to convert it.
idxnz = cell(nkpts, 1);
for ik = 1:nkpts
  idxnz{ik} = ggrid.idxnz;
end

reciprocal_grid_info = struct('wfncut', ha2ry*ggrid.ecut, ...
'fftgrid', [mol.n1, mol.n2, mol.n3], 'vol', mol.vol);
reciprocal_grid_info.xyz = xyz;
reciprocal_grid_info.idxnz = idxnz;  

% Transform
% Prepare groundstate struct and save
groundstate.rhor = rhor;
groundstate.Vxc = Vxc;
groundstate.ev = ev;
psig_cell = cell(nkpts, nspin);
for ispin = 1:nspin
  for ikpt = 1:nkpts
    psig_cell{ikpt, ispin} = psig;
  end
end
groundstate.psig = psig_cell;
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
