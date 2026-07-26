function data = load_formal_groundstate(dirin, config)
%LOAD_FORMAL_GROUNDSTATE  Synthetic supercell groundstate for complexity benchmarks.
%
% Reads unit-cell QE xml only (no wfc/charge hdf5), scales lattice/FFT grid by
% SUPERCELL [k1,k2,k3], forces nsym=1 (no crystal symmetry), Gamma-only k mesh,
% random G-space wavefunctions and monotonic Eo.

  if nargin < 2 || ~isstruct(config)
    error('load_formal_groundstate:config', 'config struct is required.');
  end
  if ~isfield(config, 'SUPERCELL')
    error('load_formal_groundstate:SUPERCELL', 'FORMAL mode requires &SUPERCELL.');
  end
  sc = config.SUPERCELL;
  ratio = int32([sc.k1, sc.k2, sc.k3]);
  if any(ratio < 1)
    error('load_formal_groundstate:ratio', 'SUPERCELL k1/k2/k3 must be >= 1.');
  end

  formal = local_formal_opts(config);
  meta = local_read_qe_xml_metadata(dirin);

  n1 = meta.n1 * double(ratio(1));
  n2 = meta.n2 * double(ratio(2));
  n3 = meta.n3 * double(ratio(3));
  supercell = meta.supercell;
  supercell(1, :) = supercell(1, :) * double(ratio(1));
  supercell(2, :) = supercell(2, :) * double(ratio(2));
  supercell(3, :) = supercell(3, :) * double(ratio(3));
  vol = det(supercell);
  nelec_sc = meta.nelec * double(prod(ratio));

  nb = local_resolve_nb(config, meta.nbnd, nelec_sc);
  nv = max(1, round(nelec_sc / 2));

  gvecinput = struct();
  gvecinput.n1 = n1;
  gvecinput.n2 = n2;
  gvecinput.n3 = n3;
  gvecinput.ecut = meta.ecutwfc;
  gvecinput.supercell = supercell;
  gvecinput.qpoint = [0, 0, 0];
  gv = gvec(gvecinput);
  ng = gv.ng;
  idxnz = {gv.idxnz};
  mill = gv.components;

  rng(formal.wf_seed);
  psig = cell(1, meta.nspin);
  psig{1, 1} = (randn(ng, nb) + 1i * randn(ng, nb)) / sqrt(2);
  if meta.nspin == 2
    psig{1, 2} = (randn(ng, nb) + 1i * randn(ng, nb)) / sqrt(2);
  end

  ev = zeros(nb, 1, meta.nspin);
  for ib = 1:nb
    ev(ib, 1, :) = formal.eo_e0 + (ib - 1) * formal.eo_delta;
  end

  occ = zeros(nb, 1, meta.nspin);
  occ(1:min(nv, nb), 1, :) = 2.0;

  sys = struct();
  sys.ng = ng;
  sys.nr = n1 * n2 * n3;
  sys.ne = nelec_sc;
  sys.n1 = n1;
  sys.n2 = n2;
  sys.n3 = n3;
  sys.supercell = supercell;
  sys.qk = [0, 0, 0];
  sys.vol = vol;

  reciprocal_grid_info = struct();
  reciprocal_grid_info.fftgrid = [n1, n2, n3];
  reciprocal_grid_info.vol = vol;
  reciprocal_grid_info.idxnz = idxnz;
  reciprocal_grid_info.wfncut = meta.ecutwfc;
  reciprocal_grid_info.xyz = mill;

  data = struct();
  data.rhor = ones(n1, n2, n3) * (nelec_sc / vol);
  data.Vxc = zeros(nb, 1, meta.nspin);
  data.ev = ev;
  data.psig = psig;
  data.sys = sys;
  data.occupation = occ;
  data.reciprocal_grid_info = reciprocal_grid_info;
  data.nkibz = 1;
  data.kibz = [0, 0, 0];
  data.kweight = 1;
  data.nspin = meta.nspin;
  data.nspinor = meta.nspinor;
  data.xyz = meta.xyz;
  data.atom_symbol = meta.atom_symbol;
  data.syms = local_trivial_syms();

  fprintf(['[FORMAL] unit cell %s -> supercell [%d %d %d], fft [%d %d %d], ', ...
    'nb=%d, ng=%d, ne=%g, nsym=1 (no symmetry)\n'], ...
    dirin, ratio(1), ratio(2), ratio(3), n1, n2, n3, nb, ng, nelec_sc);
end

function formal = local_formal_opts(config)
  defaults = default_param_values().FORMAL;
  formal = defaults;
  if isfield(config, 'FORMAL') && isstruct(config.FORMAL)
    keys = fieldnames(config.FORMAL);
    for i = 1:numel(keys)
      formal.(keys{i}) = config.FORMAL.(keys{i});
    end
  end
end

function nb = local_resolve_nb(config, xml_nbnd, nelec_sc)
  nb = [];
  if isfield(config, 'SYSTEM') && isfield(config.SYSTEM, 'energy_band_index_max') ...
      && config.SYSTEM.energy_band_index_max > 0
    nb = round(double(config.SYSTEM.energy_band_index_max));
  end
  if isempty(nb)
    nv = max(1, round(nelec_sc / 2));
    nb = max(double(xml_nbnd), 2 * nv);
  end
  if nb < 1
    error('load_formal_groundstate:nb', 'Resolved nb < 1.');
  end
end

function meta = local_read_qe_xml_metadata(qepath)
  xmlname = fullfile(qepath, 'data-file-schema.xml');
  if ~isfile(xmlname)
    error('load_formal_groundstate:xml', 'Missing %s', xmlname);
  end
  doc = xmlread(xmlname);
  out = doc.getElementsByTagName('output').item(0);

  fft_grid = out.getElementsByTagName('fft_grid').item(0);
  meta.n1 = str2double(fft_grid.getAttribute('nr1'));
  meta.n2 = str2double(fft_grid.getAttribute('nr2'));
  meta.n3 = str2double(fft_grid.getAttribute('nr3'));

  cellobj = out.getElementsByTagName('cell').item(0);
  a1 = str2num(cellobj.getElementsByTagName('a1').item(0).getTextContent); %#ok<ST2NM>
  a2 = str2num(cellobj.getElementsByTagName('a2').item(0).getTextContent); %#ok<ST2NM>
  a3 = str2num(cellobj.getElementsByTagName('a3').item(0).getTextContent); %#ok<ST2NM>
  meta.supercell = [a1; a2; a3];
  meta.ecutwfc = str2double(out.getElementsByTagName('ecutwfc').item(0).getTextContent);

  lsda = out.getElementsByTagName('magnetization').item(0) ...
    .getElementsByTagName('lsda').item(0).getTextContent;
  noncolin = out.getElementsByTagName('magnetization').item(0) ...
    .getElementsByTagName('noncolin').item(0).getTextContent;
  if strcmpi(lsda, 'true')
    meta.nspin = 2;
    if strcmpi(noncolin, 'true')
      meta.nspinor = 2;
    else
      meta.nspinor = 1;
    end
  else
    meta.nspin = 1;
    if strcmpi(noncolin, 'true')
      meta.nspinor = 2;
    else
      meta.nspinor = 1;
    end
  end

  meta.nelec = str2double(out.getElementsByTagName('nelec').item(0).getTextContent);
  if meta.nspin == 2
    nbnd_up = str2double(out.getElementsByTagName('nbnd_up').item(0).getTextContent);
    nbnd_dw = str2double(out.getElementsByTagName('nbnd_dw').item(0).getTextContent);
    meta.nbnd = (nbnd_up + nbnd_dw) / 2;
  else
    meta.nbnd = str2double(out.getElementsByTagName('nbnd').item(0).getTextContent);
  end

  structure = out.getElementsByTagName('atomic_positions').item(0) ...
    .getElementsByTagName('atom');
  xyzlist = [];
  meta.atom_symbol = cell(structure.getLength, 1);
  for i = 0:structure.getLength - 1
    xyzlist = [xyzlist; str2num(structure.item(i).getTextContent)]; %#ok<AGROW,ST2NM>
    meta.atom_symbol{i + 1} = char(structure.item(i).getAttribute('name'));
  end
  meta.xyz = xyzlist;
end

function syms = local_trivial_syms()
  syms = struct();
  syms.nsym = 1;
  syms.is_t_rev = int32(0);
  syms.nrot = 1;
  syms.ntranq = 0;
  syms.mtrx = {eye(3)};
  syms.indsub = 0;
  syms.kgzero = zeros(1, 3);
end
