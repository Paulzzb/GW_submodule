function GWinfor = construct_GWinfo(data, config)
% Convert from data to @GWinfor
%         data contains fields: 'rhor', 'Vxc', 'ev', 'psig', 'sys', 'occupation',
%                               'reciprocal_grid_info'.

% Step 1: Extract data
nkibz = data.nkibz;
nspin = data.nspin;
nspinor = data.nspinor;



rhor = data.rhor;
Vxc = data.Vxc;
ev = data.ev;
psig = data.psig;
sys = data.sys;
occupation = data.occupation;
reciprocal_grid_info = data.reciprocal_grid_info;
kibz = data.kibz;

% Step 2: Convert to @GWinfo
% @GWinfo contains fields: coulG, coulG0, supercell, bdot, qk, ntot(number of total band calculated),
% vol, ne, nv, gvec, gvec2(double Ecut), gvecrho(for density), Vxc, rho(density in reciprocal space),
% ev, Z
% abandoned: aqs, idxnz
ha2ry = 2.0;
ry2ev = 13.60569253;
GWinfor = GWinfo();


GWinfor.nspin = nspin;
GWinfor.nspinor = nspinor;

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

GWinfor.gvec_list = cell(nkibz, 1);
for ik = 1:nkibz
  gvecinput = [];
  gvecinput.n1 = sys.n1;
  gvecinput.n2 = sys.n2;
  gvecinput.n3 = sys.n3;
  gvecinput.ecut = config.CUTOFFS.coulomb_cutoff;
  gvecinput.supercell = sys.supercell;
  gvecinput.qpoint = kibz(ik, :);
  GWinfor.gvec_list{ik} = gvec(gvecinput);
end
% GWinfor.idxnz = idxnz; 
 
% Compute Coulomb interaction
coulG_list = construct_vcoul_k(data, config, GWinfor.gvec_list);
coulG0 = construct_coulG0(data, config);

% set coulG to fit single-kpoints code
if (config.CONTROL.enable_k_points > 0)
  coulG = coulG_list{1};
  fprintf("Formly, fill coulG with coulG_list{1}")
else
  coulG = coulG_list{1};
end

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
bvec = 2*pi*inv(sys.supercell);
GWinfor.bvec = bvec;
GWinfor.bdot = GWinfor.bvec * GWinfor.bvec.';

% Symmetric information if needed

if config.CONTROL.enable_k_points
  syms = data.syms;
  % Get symmetry information
  GWinfor.symminfo = symminfo(syms.ntran, syms.ntranq, syms.mtrx, syms.nrot, syms.indsub, syms.kgzero);
  % Use irreducible k-points to generate full-kpoints set
  [GWinfor.bz_samp] = bz_sampling(data.nkibz, data.kibz, bvec, data.kweight);
  GWinfor.bz_samp = fullbz(GWinfor.bz_samp, GWinfor.symminfo, GWinfor.gvec);
  % Mapping G-Go with fft_grid
  GWinfor.mapping = setmapping(GWinfor);
  % GWinfor.gvec = setGomap(GWinfor.gvec, GWinfor.bz_samp.nGo);
  % Mapping rotation
  % GWinfor.gvec = setsymmGmap(GWinfor);
end

% Construct rhoG if needed
% if config.CONTROL.enable_density == 1
% end


end
