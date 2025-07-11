function GWinfor = construct_GWinfo(data, config)
% Convert from data to @GWinfor
%         data contains fields: 'rhor', 'Vxc', 'ev', 'psig', 'sys', 'occupation',
%                               'reciprocal_grid_info'.

% Step 1: Extract data
nkpts = data.nkpts;
nspin = data.nspin;
nspinor = data.nspinor;


rhor = data.rhor;
Vxc = data.Vxc;
ev = data.ev;
psig = data.psig;
sys = data.sys;
occupation = data.occupation;
reciprocal_grid_info = data.reciprocal_grid_info;
kpts = data.kpts;

% Step 2: Convert to @GWinfo
% @GWinfo contains fields: coulG, coulG0, supercell, bdot, qk, ntot(number of total band calculated),
% vol, ne, nv, gvec, gvec2(double Ecut), gvecrho(for density), Vxc, rho(density in reciprocal space),
% ev, Z
% abandoned: aqs, idxnz
ha2ry = 2.0;
ry2ev = 13.60569253;
GWinfor = GWinfo();



% Compute other values
bmatrix = 2*pi*inv(sys.supercell');
bdot = bmatrix * bmatrix';

% Compute @gvec to GWinfor.gvec
gvecinput = [];
gvecinput.n1 = sys.n1;
gvecinput.n2 = sys.n2;
gvecinput.n3 = sys.n3;
gvecinput.ecut = config.CUTOFFS.coulomb_cutoff;
gvecinput.supercell = sys.supercell;
GWinfor.gvec = gvec(gvecinput);

GWinfor.gvec_list = cell(nkpts, 1);
for ik = 1:nkpts
  gvecinput = [];
  gvecinput.n1 = sys.n1;
  gvecinput.n2 = sys.n2;
  gvecinput.n3 = sys.n3;
  gvecinput.ecut = config.CUTOFFS.coulomb_cutoff;
  gvecinput.supercell = sys.supercell;
  gvecinput.qpoint = kpts(ik, :);
  GWinfor.gvec_list{ik} = gvec(gvecinput);
end
% GWinfor.idxnz = idxnz; 
 
% Compute Coulomb interaction
coulG_list = construct_vcoul_k(data, config, GWinfor.gvec_list);
coulG0 = construct_coulG0(data, config);
coulG = construct_vcoul(data, config);

% If config.FREQUENCY.frequency_dependence == 1 --> GPP approximation
% Prepare rho(G) using rho(R) and config.CUTOFFS.density_cutoff
if config.FREQUENCY.frequency_dependence == 1
  [rhoG, gvecrho] = construct_rhoG(data, config);
  GWinfor.gvecrho = gvecrho;
  GWinfor.rho = rhoG;
end

GWinfor.coulG = coulG;
GWinfor.coulG_list = coulG_list;
GWinfor.coulG0 = coulG0;

GWinfor.supercell = sys.supercell;
GWinfor.bdot = bdot;
GWinfor.vol = sys.vol;
GWinfor.ne = sys.ne;

GWinfor.Vxc = data.Vxc * ha2ry;
GWinfor.ev = data.ev * ha2ry;
GWinfor.psig = psig;
GWinfor.Ggrid4psig = data.reciprocal_grid_info;
GWinfor.occupation = data.occupation;


end
