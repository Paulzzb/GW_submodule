% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/06/12

function report = testfunc(config, id)
%TESTFUNC  Experimental Gamma tildeVq builder (config-driven, mu-batched).
%
%   report = testfunc(config)
%   report = testfunc(config, id)   % build one ISDF pool slot (driver hook)
%
% Requires service modules initialized (wave_functions, isdf, coulomb, ...).
% Optional &TESTFUNC in input:
%   mchq_batch_mib   - real-space MCHq batch cap (default 512 MiB)
%   save_to_pool     - write tildeVq/helperqG back to ISDF pool (default true)

  if nargin < 1 || ~isstruct(config)
    error('testfunc:config', 'config struct is required.');
  end
  opts = local_testfunc_opts(config);

  if nargin >= 2 && ~isempty(id)
    t0 = tic;
    one = local_gamma_tildeVq(config, int32(id), opts);
    one.wall_s = toc(t0);
    report = one;
    return
  end

  ids = local_ids_to_build(config);
  if isempty(ids)
    error('testfunc:ids', 'No ISDF slot to build (check &ISDF compute_* flags).');
  end
  report = repmat(local_empty_case(), 0, 1);
  for k = 1:numel(ids)
    t0 = tic;
    report(end + 1, 1) = local_gamma_tildeVq(config, ids(k), opts); %#ok<AGROW>
    report(end).wall_s = toc(t0);
  end
end

function opts = local_testfunc_opts(config)
  opts = struct();
  opts.mchq_batch_bytes = 512 * 1024^2;
  opts.save_to_pool = true;
  if ~isfield(config, 'TESTFUNC') || ~isstruct(config.TESTFUNC)
    return
  end
  tf = config.TESTFUNC;
  if isfield(tf, 'mchq_batch_mib') && ~isempty(tf.mchq_batch_mib)
    opts.mchq_batch_bytes = double(tf.mchq_batch_mib) * 1024^2;
  end
  if isfield(tf, 'save_to_pool') && ~isempty(tf.save_to_pool)
    opts.save_to_pool = logical(tf.save_to_pool);
  end
end

function ids = local_ids_to_build(config)
  ids = int32([]);
  if ~isfield(config, 'ISDF') || ~isstruct(config.ISDF)
    [id_vc, id_vn, id_nn] = isdf.cohsex_resolve_ids(config);
    ids = local_nonempty_ids([id_vc, id_vn, id_nn]);
    return
  end
  isdf_cfg = config.ISDF;
  want = {};
  if isfield(isdf_cfg, 'compute_vc') && logical(isdf_cfg.compute_vc)
    want{end + 1} = 'vc'; %#ok<AGROW>
  end
  if isfield(isdf_cfg, 'compute_vn') && logical(isdf_cfg.compute_vn)
    want{end + 1} = 'vn'; %#ok<AGROW>
  end
  if isfield(isdf_cfg, 'compute_nn') && logical(isdf_cfg.compute_nn)
    want{end + 1} = 'nn'; %#ok<AGROW>
  end
  [id_vc, id_vn, id_nn] = isdf.cohsex_resolve_ids(config);
  map = struct('vc', id_vc, 'vn', id_vn, 'nn', id_nn);
  for k = 1:numel(want)
    label = want{k};
    if isfield(map, label) && ~isempty(map.(label))
      ids(end + 1) = int32(map.(label)); %#ok<AGROW>
    end
  end
  ids = unique(ids, 'stable');
end

function ids = local_nonempty_ids(id_list)
  ids = int32([]);
  for k = 1:numel(id_list)
    v = id_list(k);
    if ~isempty(v)
      ids(end + 1) = int32(v); %#ok<AGROW>
    end
  end
end

function case_out = local_empty_case()
  case_out = struct('id', int32(0), 'desc', "", 'Nmu', 0, 'ng', 0, ...
    'nc', 0, 'n_batch', 0, 'wall_s', 0, 'tildeVq', []);
end

function case_out = local_gamma_tildeVq(config, id, opts)
  wf_data = wave_functions.get();
  k_data = lattice.manager('k', 'get');
  r_lat_data = lattice.manager('r_lat', 'get');
  d_lat_data = lattice.manager('d_lat', 'get');
  coul_data = coulomb.get();
  fft_data = FFT.get();
  fft_sz = fft_data.fftgrid;
  isdf_data = isdf.get(id);

  if isempty(isdf_data.nrange1) || isempty(isdf_data.nrange2)
    isdf.set_nrange(id, config.SYSTEM);
    isdf_data = isdf.get(id);
  end
  if isempty(isdf_data.nrange1) || isempty(isdf_data.nrange2)
    error('testfunc:nrange', ...
      'Missing nrange for ISDF id=%d; run isdf.set_nrange first.', int32(id));
  end

  nrange1 = double(isdf_data.nrange1);
  nrange2 = double(isdf_data.nrange2);
  nb1 = numel(nrange1);
  nb2 = numel(nrange2);

  DL_vol = d_lat_data.DL_vol;
  nc = double(wf_data.nc);
  ng = double(coul_data.coulomb_ng);
  nibz = int32(k_data.nibz);
  nbz = int32(k_data.nbz);
  R_sampling_RLU = isdf_data.bundle_struct.R_grid_bundle(isdf_data.bundle_struct.sampling2bundle, :);
  Nmu = double(isdf_data.bundle_struct.N_sampling);
  n_batch = local_mchq_batch_size(nc, Nmu, opts.mchq_batch_bytes);
  use_parfor = parallel.enabled();
  g_table_col = fft_data.G_table(:, 1);

  fprintf(['[testfunc] id=%d desc=%s | nc=%.3g ng=%d Nmu=%.0f ', ...
    'nb1=%d nb2=%d nbz=%d | mchq_batch=%d (~%.2f MiB/batch)\n'], ...
    int32(id), char(string(isdf_data.desc)), nc, ng, Nmu, nb1, nb2, nbz, ...
    n_batch, n_batch * nc * 8 / 1024^2);

  tildeVq = zeros(Nmu, Nmu, nibz);
  helperqG = zeros(ng, Nmu);
  isdf_data.helperqG = zeros(ng, Nmu, nibz);

  iqibz = 1;
  vcoul_q = coul_data.vcoul(:, iqibz);
  if iqibz == 1
    vcoul_q(1) = coul_data.vcoul0;
  end

  MCHq_G = zeros(ng, Nmu);
  CCHq = zeros(Nmu, Nmu);
  for i_start = 1:n_batch:Nmu
    i_end = min(i_start + n_batch - 1, Nmu);
    idx_batch = i_start:i_end;
    MCHq_batched_R = zeros(nc, numel(idx_batch));
    accumulate_cchq = (i_start == 1);

    for ikbz = 1:nbz
      ikibz = k_data.bz2ibz(ikbz, 1);
      ikrot = k_data.bz2rot(ikbz, 1);
      ikpbz = r_lat_data.qindx_X(iqibz, ikbz, 1);
      iGo = r_lat_data.qindx_X(iqibz, ikbz, 2);
      ikpibz = k_data.bz2ibz(ikpbz, 1);
      ikprot = k_data.bz2rot(ikpbz, 1);

      wf_1 = zeros(nc, nb1);
      wf_2 = zeros(nc, nb2);
      wf_c_1 = zeros(Nmu, nb1);
      wf_c_2 = zeros(Nmu, nb2);
      for ind_ib = 1:nb1
        ib = nrange1(ind_ib);
        isc1 = int32([ib, ikibz, ikrot, 1]);
        wf_1(:, ind_ib) = wave_functions.WF_apply_symm(isc1);
        wf_c_1(:, ind_ib) = isdf.get_u_xalpha(id, isc1);
      end
      for ind_ib = 1:nb2
        ib = nrange2(ind_ib);
        isc2 = int32([ib, ikpibz, ikprot, 1]);
        wf_2(:, ind_ib) = wave_functions.WF_apply_symm(isc2);
        wf_c_2(:, ind_ib) = isdf.get_u_xalpha(id, isc2);
      end
      wf_1 = conj(wf_1);
      wf_c_1 = conj(wf_c_1);

      Go = double(r_lat_data.Ggrid_RLU(iGo, :));
      phase_shift_sampling = double(exp(-2 * pi * 1i * (R_sampling_RLU ./ double(fft_data.fftgrid)) * Go'));
      phase_shift_fftgrid = double(exp(-2 * pi * 1i * (fft_data.Rgrid_RLU ./ double(fft_data.fftgrid)) * Go'));

      if accumulate_cchq
        CCHq_tmp = (wf_c_1 * wf_c_1') .* (wf_c_2 * wf_c_2') ...
          .* phase_shift_sampling' .* phase_shift_sampling;
        CCHq = CCHq + CCHq_tmp;
      end

      wf_c_1_batch = wf_c_1(idx_batch, :);
      wf_c_2_batch = wf_c_2(idx_batch, :);
      ph_s = phase_shift_sampling(idx_batch).';
      MCHq_tmp = (wf_1 * wf_c_1_batch') .* (wf_2 * wf_c_2_batch') ...
        .* ph_s .* phase_shift_fftgrid;
      MCHq_batched_R = MCHq_batched_R + MCHq_tmp;
    end

    batch_G = local_fft_mchq_batch(MCHq_batched_R, fft_sz, g_table_col, DL_vol, use_parfor);
    MCHq_G(:, idx_batch) = batch_G;
  end

  CCHq = 0.5 * (CCHq + CCHq');
  helperqG = MCHq_G / CCHq;
  isdf_data.helperqG(:, :, iqibz) = helperqG;
  tildeVq(:, :, iqibz) = helperqG' * diag(vcoul_q) * helperqG;
  tildeVq(:, :, iqibz) = 0.5 * (tildeVq(:, :, iqibz) + tildeVq(:, :, iqibz)');

  if opts.save_to_pool
    isdf_data.tildeVq = tildeVq;
    isdf.save2mod(isdf_data, id);
  end

  case_out = struct();
  case_out.id = int32(id);
  case_out.desc = string(isdf_data.desc);
  case_out.Nmu = Nmu;
  case_out.ng = ng;
  case_out.nc = nc;
  case_out.n_batch = n_batch;
  case_out.wall_s = 0;
  case_out.tildeVq = tildeVq;
end

function n_batch = local_mchq_batch_size(nc, Nmu, max_bytes)
  n_batch = max(1, min(Nmu, floor(max_bytes / (nc * 8))));
end

function batch_G = local_fft_mchq_batch(MCHq_batched_R, fft_sz, g_table_col, DL_vol, use_parfor)
  n_bat = size(MCHq_batched_R, 2);
  ng = numel(g_table_col);
  batch_G = zeros(ng, n_bat);
  if use_parfor
    parfor ib = 1:n_bat
      fftbox = reshape(MCHq_batched_R(:, ib), fft_sz);
      fftbox = do_FFT(fftbox, fft_sz, 1) * DL_vol;
      batch_G(:, ib) = fftbox(g_table_col);
    end
  else
    for ib = 1:n_bat
      fftbox = reshape(MCHq_batched_R(:, ib), fft_sz);
      fftbox = do_FFT(fftbox, fft_sz, 1) * DL_vol;
      batch_G(:, ib) = fftbox(g_table_col);
    end
  end
end
