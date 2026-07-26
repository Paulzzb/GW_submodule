% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/02

function tildeVq = gen_tildeVq(id, cfg_isdf)
% GEN_TILDEVQ  SVD-based tildeVq with X*V*X' form.
%
% Accumulates MCHq/CCHq over the BZ, forms truncated C^{-1/2}, stores CCHq and CCHq_inv_sqrt
% per q for isdf.get_rho_xalpha (c_t = C^{-1/2} c).
% cfg_isdf (optional): config.ISDF
%   trunc = inv_param, ratio = inv_ratio.

  wf_data = wave_functions.get();
  k_data = lattice.manager('k', 'get');
  r_lat_data = lattice.manager('r_lat', 'get');
  d_lat_data = lattice.manager('d_lat', 'get');
  coul_data = coulomb.get();
  fft_data = FFT.get();
  fft_sz = fft_data.fftgrid;
  isdf_data = isdf.get(id);
  if isempty(isdf_data.nrange1) || isempty(isdf_data.nrange2)
    error('isdf:gen_tildeVq:nrange', ...
      'Missing cached nrange in ISDF id=%d. Build it first via isdf.set_nrange(id, config.SYSTEM).', ...
      int32(id));
  end
  nrange1 = double(isdf_data.nrange1);
  nrange2 = double(isdf_data.nrange2);
  nb1 = length(nrange1); nb2 = length(nrange2);

  DL_vol = d_lat_data.DL_vol;
  nc = wf_data.nc;
  ng = coul_data.coulomb_ng;
  nibz = int32(k_data.nibz);
  nbz = int32(k_data.nbz);
  R_sampling_RLU = isdf_data.bundle_struct.R_grid_bundle(isdf_data.bundle_struct.sampling2bundle, :);
  
  Nmu = isdf_data.bundle_struct.N_sampling;
  tildeVq = zeros(Nmu, Nmu, nibz);
  isdf_data.helperqG = zeros(ng, Nmu, nibz);
  isdf_data.CCHq = zeros(Nmu, Nmu, nibz);
  isdf_data.CCHq_inv_sqrt = zeros(Nmu, Nmu, nibz);
  isdf_data.CCHq_trunc_factors = cell(1, double(nibz));
  if nargin < 2
    cfg_isdf = [];
  end
  s_cut = 0;
  ratio = 0.5;
  if isstruct(cfg_isdf) && isfield(cfg_isdf, 'inv_param') && ~isempty(cfg_isdf.inv_param)
    s_cut = double(cfg_isdf.inv_param);
  end
  if isstruct(cfg_isdf) && isfield(cfg_isdf, 'inv_ratio') && ~isempty(cfg_isdf.inv_ratio)
    ratio = double(cfg_isdf.inv_ratio);
  end
  isdf_data.svd_s_cut = s_cut;
  isdf_data.svd_ratio = ratio;
  isdf.numerical_cond_report('gen_tildeVq_begin', id, isdf_data.desc, s_cut);
  use_parfor = parallel.enabled();
  g_table_col = fft_data.G_table(:, 1);

  for iqibz = 1:nibz
    helperqG = zeros(ng, Nmu);
    vcoul_q = coul_data.vcoul(:, iqibz);
    if iqibz == 1
      vcoul_q(1) = coul_data.vcoul0;
    end

    MCHq = zeros(nc, Nmu);
    CCHq = zeros(Nmu, Nmu);
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
      %
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
      phase_shift_sampling = double( exp(- 2 * pi * 1i * (R_sampling_RLU ./ double(fft_data.fftgrid)) * Go') );
      phase_shift_fftgrid   = double( exp(- 2 * pi * 1i * (fft_data.Rgrid_RLU ./ double(fft_data.fftgrid)) * Go') );

      MCHq_tmp = (wf_1 * wf_c_1') .* (wf_2 * wf_c_2') ...
        .* phase_shift_sampling' .* phase_shift_fftgrid;
      MCHq = MCHq + MCHq_tmp;
      CCHq_tmp = (wf_c_1 * wf_c_1') .* (wf_c_2 * wf_c_2') ...
        .* phase_shift_sampling' .* ( phase_shift_sampling );
      CCHq = CCHq + CCHq_tmp;
    end % ikbz


    l_keep = isdf.prod_C_inv_t('set', CCHq, s_cut, ratio);
    isdf_data.CCHq(:, :, iqibz) = CCHq;
    [V_trunc_iq, Lambda_trunc_iq] = isdf.prod_C_inv_t('get_factors');
    fac = struct();
    fac.V_trunc = V_trunc_iq;
    fac.Lambda_trunc = Lambda_trunc_iq;
    fac.N_keep = int32(l_keep);
    isdf_data.CCHq_trunc_factors{double(iqibz)} = fac;
    isdf.numerical_cond_report('gen_tildeVq_iq', id, iqibz, CCHq, MCHq, l_keep, s_cut);

    helperqR = isdf.prod_C_inv_t('prod', 'r', MCHq);

    if use_parfor
      parfor imu = 1:Nmu
        fftbox = reshape(helperqR(:, imu), fft_sz);
        fftbox = do_FFT(fftbox, fft_sz, 1) * DL_vol;
        helperqG(:, imu) = fftbox(g_table_col);
      end
    else
      for imu = 1:Nmu
        fftbox = reshape(helperqR(:, imu), fft_sz);
        fftbox = do_FFT(fftbox, fft_sz, 1) * DL_vol;
        helperqG(:, imu) = fftbox(g_table_col);
      end
    end

    isdf_data.helperqG(:, :, iqibz) = helperqG;
    tildeVq(:, :, iqibz) = helperqG' * diag(vcoul_q) * helperqG;
    tildeVq(:, :, iqibz) = tildeVq(:, :, iqibz) / 2 + tildeVq(:, :, iqibz)' / 2;
  end % iqibz


  isdf_data.tildeVq = tildeVq;
  isdf.save2mod(isdf_data, id);
  fprintf('[ISDF] gen_tildeVq done for id=%d (SVD C^{-1/2} path, s_cut=%g).\n', int32(id), s_cut);

end
