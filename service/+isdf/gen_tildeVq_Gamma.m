% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/06/11 

function gen_tildeVq_Gamma(id, config)
%GEN_TILDEVQ_GAMMA  Gamma-only tildeVq builder (mu-batched MCHq).
%
%   tildeVq = isdf.gen_tildeVq_Gamma(id)
%   tildeVq = isdf.gen_tildeVq_Gamma(id, config)
%
% Used when frequency_dependence == -2. Writes tildeVq / helperqG to the pool.
% Optional config.TESTFUNC.mchq_batch_mib (default 512) caps real-space MCHq batch.

  if nargin < 2
    config = [];
  end

  mchq_batch_bytes = 512 * 1024^2;
  if isstruct(config) && isfield(config, 'TESTFUNC') && isstruct(config.TESTFUNC) ...
      && isfield(config.TESTFUNC, 'mchq_batch_mib') && ~isempty(config.TESTFUNC.mchq_batch_mib)
    mchq_batch_bytes = double(config.TESTFUNC.mchq_batch_mib) * 1024^2;
  end

  wf_data = wave_functions.get();
  k_data = lattice.manager('k', 'get');
  r_lat_data = lattice.manager('r_lat', 'get');
  d_lat_data = lattice.manager('d_lat', 'get');
  coul_data = coulomb.get();
  fft_data = FFT.get();
  fft_sz = fft_data.fftgrid;
  isdf_data = isdf.get(id);

  if isempty(isdf_data.nrange1) || isempty(isdf_data.nrange2)
    if isstruct(config) && isfield(config, 'SYSTEM')
      isdf.set_nrange(id, config.SYSTEM);
      isdf_data = isdf.get(id);
    end
  end
  if isempty(isdf_data.nrange1) || isempty(isdf_data.nrange2)
    error('isdf:gen_tildeVq_Gamma:nrange', ...
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
  n_batch = max(1, min(Nmu, floor(mchq_batch_bytes / (nc * 8))));
  use_parfor = parallel.enabled();
  g_table_col = fft_data.G_table(:, 1);

  output.msg('r', ['gen_tildeVq_Gamma id=%d desc=%s | nc=%.3g ng=%d Nmu=%.0f ', ...
    'nb1=%d nb2=%d nbz=%d | mchq_batch=%d (~%.2f MiB/batch)'], ...
    int32(id), char(string(isdf_data.desc)), nc, ng, Nmu, nb1, nb2, nbz, ...
    n_batch, n_batch * nc * 8 / 1024^2);

  tildeVq = zeros(Nmu, Nmu, nibz);
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

    n_bat = size(MCHq_batched_R, 2);
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
    MCHq_G(:, idx_batch) = batch_G;
  end

  CCHq = 0.5 * (CCHq + CCHq');
  helperqG = MCHq_G / CCHq;
  isdf_data.helperqG(:, :, iqibz) = helperqG;
  tildeVq(:, :, iqibz) = helperqG' * diag(vcoul_q) * helperqG;
  tildeVq(:, :, iqibz) = 0.5 * (tildeVq(:, :, iqibz) + tildeVq(:, :, iqibz)');

  isdf_data.tildeVq = tildeVq;
  isdf.save2mod(isdf_data, id);
  output.msg('r', 'gen_tildeVq_Gamma done for id=%d.', int32(id));
end
