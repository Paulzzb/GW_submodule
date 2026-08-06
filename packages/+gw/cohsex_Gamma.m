%
% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/06 ZZ
function [Esx_x, Ecoh] = cohsex_Gamma(config)
  %COHSEX_GAMMA  Static COHSEX self-energy diagonal at Gamma.
  %
  %   [Esx_x, Ecoh] = gw.cohsex_Gamma(config)
  %
  % ISDF path uses service/+isdf (gen_Kq_Gamma, gen_tildeWq_Gamma, ...).
  % Dense no-ISDF path is for G-space validation.
  %
  % Prerequisites (ISDF):
  %   - Service layer ready (wave_functions, system, lattice, FFT, coulomb, isdf).
  %   - ISDF slots desc 'vc' / 'vn' / 'nn' as required by exact_ch.
  %   - nqbz == 1; ispin = 1 (multi-spin not supported).

    cleanup = output.push('+gw/cohsex_Gamma.m'); %#ok<NASGU>
    output.msg('v0s', 'Start computing Esx_x / Ecoh (Gamma).');
    tStart = tic;

    default_Constant = constant_map();
    nameConstants = fieldnames(default_Constant);
    for i = 1:numel(nameConstants)
      eval(sprintf('%s = %.16f;', nameConstants{i}, default_Constant.(nameConstants{i})));
    end

    nbmin = config.SYSTEM.energy_band_index_min;
    nbmax = config.SYSTEM.energy_band_index_max;
    nband = nbmax - nbmin + 1;

    system_data = system.get();
    k_data = lattice.manager('k', 'get');
    q_data = lattice.manager('q', 'get');
    r_lat_data = lattice.manager('r_lat', 'get');
    d_lat_data = lattice.manager('d_lat', 'get');
    coulomb_data = coulomb.get();
    fft_data = FFT.get();

    nkibz = k_data.nibz;
    nqbz = q_data.nbz;

    if nqbz ~= 1
      output.err('nqbz must be 1 for the Gamma-point path.');
    end
    fft_sz = double(fft_data.fftgrid(:).');
    DL_vol = double(lattice.manager('d_lat', 'get').DL_vol);
    dv = d_lat_data.DL_vol / double(fft_data.nr);

    ispin = 1;
    exact_CH = logical(config.COHSEX.exact_ch);

    if (exact_CH && ~config.ISDF.isisdf)
      msg = sprintf(['Exact COH without ISDF is not supported yet. ', ...
        'Use BerkeleyGW as a reference instead.\n']);
      output.err('%s', msg);
    end

    if (config.ISDF.isisdf)
      % ---- ISDF path: resolve COHSEX slot ids ----
      if exact_CH
        [id_vc, id_vn, ~] = isdf.resolve_ids(config);
      else
        [id_vc, id_vn, id_nn] = isdf.resolve_ids(config);
      end

      if isempty(id_vn) && isempty(id_vc)
        output.err( ...
          ['Static COHSEX requires both vc and vn ISDF slots. ' ...
           'Enable config.ISDF.compute_vc / compute_vn, or provide assigned ids.']);
      end
      if ~exact_CH && isempty(id_nn)
        output.err( ...
          ['Static COHSEX (approximate CH) requires an nn ISDF slot. ' ...
           'Enable config.ISDF.compute_nn, or provide an assigned id.']);
      end

  
      isdf.set_nrange(id_vc, config.SYSTEM);
      vc_data = isdf.get(id_vc);
      vn_data = isdf.get(id_vn);

      nrangev = double(vc_data.nrange1);
      fac_vc = vc_data.CCHq_trunc_factors{double(1)};

      fac_vn = vn_data.CCHq_trunc_factors{double(1)};
      s_ratio_vn = vn_data.svd_ratio;

      s2b_vn = vn_data.bundle_struct.sampling2bundle;

      if ~exact_CH
        nn_data = isdf.get(id_nn);
        nrangea = double(nn_data.nrange1);
        fac_nn = nn_data.CCHq_trunc_factors{double(1)};
        s_ratio_nn = nn_data.svd_ratio;
        s2b_nn = nn_data.bundle_struct.sampling2bundle;
      end
  
      iqibz = 1;
      Sigma_sex_x = zeros(nband, nkibz);
      Sigma_coh = zeros(nband, nkibz);
      
      vcoul_q = coulomb_data.vcoul(:, iqibz);
      vcoul_q(1) = coulomb_data.vcoul0;
      Kq_ISDF = isdf.gen_Kq_Gamma(id_vc, iqibz);

      tildeWq_vn = isdf.gen_tildeWq_Gamma(id_vc, iqibz, Kq_ISDF, id_vn); % SEX kernel
      helperqG_vc = double(vc_data.helperqG(:, :, iqibz));
      % Recover helperqR from helperqG (inverse of gen_tildeVq R->G).
      Nkeep_vc = fac_vc.N_keep;
      
      if exact_CH
        nr_fine = prod(fft_sz);
        L_Kq = chol(Kq_ISDF, "lower");
        helperqR_vc = zeros(nr_fine, Nkeep_vc);
        for imu = 1:Nkeep_vc
          fftbox = zeros(fft_sz);
          fftbox(fft_data.G_table(:, 1)) = helperqG_vc(:, imu) .* vcoul_q;
          fftbox = do_FFT(fftbox, fft_sz, -1) / DL_vol;
          % fftbox = do_FFT(fftbox, fft_sz, -1) * r_lat_data.RL_vol;
          helperqR_vc(:, imu) = fftbox(:);
        end
        Wqrr = zeros(nr_fine, 1);
        block_size = 64;
        coeff_delta = double(fft_data.nr);
        coeff_delta = coeff_delta * (d_lat_data.DL_vol^2 / double(fft_data.nr));
        for iblock = 1:block_size:nr_fine
          iend = min(iblock + block_size - 1, nr_fine);
          Hblk = helperqR_vc(iblock:iend, 1:Nkeep_vc);     
          tmp_blk = Hblk / L_Kq';                 
          Wqrr(iblock:iend) = 0.5 * coeff_delta * sum(abs(tmp_blk).^2, 2);
        end
        Wqrr = 0.5 * (Wqrr + conj(Wqrr));
      else
        tildeWq_nn = isdf.gen_tildeWq_Gamma(id_vc, iqibz, Kq_ISDF, id_nn);
      end

      for indib = 1:nband
        ib = nbmin + indib - 1;
        u_ib_xalpha = vn_data.bundle_struct.WF_bundle(s2b_vn, ib, 1, ispin);
        oblist = nrangev;
        u_ob_xalpha = vn_data.bundle_struct.WF_bundle(s2b_vn, oblist, 1, ispin);
        rho_left_vn = conj(u_ib_xalpha) .* u_ob_xalpha;
        rho_left_vn = diag(fac_vn.Lambda_trunc.^(s_ratio_vn-1)) * (fac_vn.V_trunc' * rho_left_vn);
        tmp1 = tildeWq_vn * rho_left_vn;
        sex_t = sum(conj(rho_left_vn) .* tmp1, 1);
        sex_t = real(sex_t);
        tmp = sum(sex_t);

        Sigma_sex_x(indib, 1) = Sigma_sex_x(indib, 1) - tmp;
        % COH contribution
        if exact_CH
          isc = [ib, 1, 1, 1];
          wf_nk = wave_functions.WF_apply_symm(isc);
          % tmp = sum(Wqrr .* abs(wf_nk).^2) * d_lat_data.DL_vol / nr_fine;
          rho_nk = abs(wf_nk).^2;
          tmp = sum(Wqrr .* rho_nk) * dv;
        else
          u_ib_xalpha = nn_data.bundle_struct.WF_bundle(s2b_nn, ib, 1, ispin);
          oblist = nrangea;
          u_ob_xalpha = nn_data.bundle_struct.WF_bundle(s2b_nn, oblist, 1, ispin);
          rho_left_nn = conj(u_ib_xalpha) .* u_ob_xalpha;
          rho_left_nn = diag(fac_nn.Lambda_trunc.^(s_ratio_nn-1)) * (fac_nn.V_trunc' * rho_left_nn);
          tmp1 = tildeWq_nn * rho_left_nn;
          sex_t = sum(conj(rho_left_nn) .* tmp1, 1);
          sex_t = real(sex_t);
          tmp = sum(sex_t);
        end
        Sigma_coh(indib, 1) = Sigma_coh(indib, 1) + tmp;
      end
      Esx_x = - Sigma_sex_x;
      Ecoh = - Sigma_coh;
    else
      % ---- Dense no-ISDF path (validation) ----
      if nkibz ~= 1 || nqbz ~= 1
        msg = sprintf( ...
          ['noISDF path currently ignores k/q loops and uses only (ikibz,iqbz)=(1,1). ', ...
           'Current nkibz=%d, nqbz=%d.\n'], nkibz, nqbz);
        output.warn('%s', msg);
      end
  
      ev_ry = double(system_data.Eo);
      focc = double(system_data.f);
      nb_total = size(ev_ry, 1);
      nv = find(focc(:, 1, 1) > 1 - TOL_SMALL, 1, 'last');
      if isempty(nv)
        output.err('Cannot determine nv from occupations.');
      end
      nsum = min(config.SYSTEM.number_bands_in_summation, nb_total);
      if nsum <= nv
        output.err( ...
          'Invalid summation bands: nv=%d, nsum=%d (need nsum > nv).', nv, nsum);
      end
  
      iqibz_ref = q_data.bz2ibz(1);
      iqrot_ref = q_data.bz2rot(1);
      iGo_ref = r_lat_data.qindx_S(1, 1, 2);
      d_lat_data = lattice.manager('d_lat', 'get');
      vol = double(d_lat_data.DL_vol);
      coul_data = coulomb.get();
      ng = size(coul_data.vcoul, 1);
      vcoul_q = double(coul_data.vcoul(:, iqibz_ref));
      if iqibz_ref == 1
        vcoul_q(1) = double(coul_data.vcoul0);
      end
      Dcoul = spdiags(vcoul_q(:), 0, ng, ng);
  
      startforW = tic;
      chi_acc = zeros(ng, ng);
      scal = 4.0;
      for ind_nv = 1:nv
        Mgvc = zeros(ng, nsum - nv);
        for jc = nv + 1:nsum
          p = struct();
          p.is = [ind_nv, 1, 1, ispin];
          p.os = [jc, 1, 1, ispin];
          p.qs = [iGo_ref, iqibz_ref, iqrot_ref];
          Mgvc(:, jc - nv) = conj(double(SCATTER_Bamp(p)));
        end
        eden = 1 ./ (ev_ry(ind_nv, 1, ispin) - ev_ry(nv + 1:nsum, 1, ispin));
        chi_acc = chi_acc + scal * Mgvc * diag(eden) * Mgvc';
      end
      inveps = eye(ng) - Dcoul * chi_acc;
      msg = sprintf('noISDF time for W = %.4f sec.\n', toc(startforW));
      output.msg('v1s', '%s', msg);
  
      startSigma = tic;
      Esx_x = zeros(nband, nkibz);
      Ecoh = zeros(nband, nkibz);
      for ioper = 1:nv
        Mgvn = zeros(ng, nband);
        for indib = 1:nband
          ib = nbmin + indib - 1;
          p = struct();
          p.is = [ioper, 1, 1, ispin];
          p.os = [ib, 1, 1, ispin];
          p.qs = [iGo_ref, iqibz_ref, iqrot_ref];
          Mgvn(:, indib) = conj(double(SCATTER_Bamp(p)));
        end
        W1Mgvn = Dcoul * Mgvn;
        W1Mgvn = inveps \ W1Mgvn;
        W1Mgvn = W1Mgvn - Dcoul * Mgvn;
        for indib = 1:nband
          diag_term = real(Mgvn(:, indib)' * W1Mgvn(:, indib));
          Esx_x(indib, 1) = Esx_x(indib, 1) - diag_term;
          Ecoh(indib, 1) = Ecoh(indib, 1) + 0.5 * diag_term;
        end
      end
      for ioper = nv + 1:nsum
        Mgcn = zeros(ng, nband);
        for indib = 1:nband
          ib = nbmin + indib - 1;
          p = struct();
          p.is = [ioper, 1, 1, ispin];
          p.os = [ib, 1, 1, ispin];
          p.qs = [iGo_ref, iqibz_ref, iqrot_ref];
          Mgcn(:, indib) = conj(double(SCATTER_Bamp(p)));
        end
        W1Mgcn = Dcoul * Mgcn;
        W1Mgcn = inveps \ W1Mgcn;
        W1Mgcn = W1Mgcn - Dcoul * Mgcn;
        for indib = 1:nband
          diag_term = real(Mgcn(:, indib)' * W1Mgcn(:, indib));
          Ecoh(indib, 1) = Ecoh(indib, 1) + 0.5 * diag_term;
        end
      end
      msg = sprintf('noISDF time for Esx_x / Ecoh = %.4f sec.\n', toc(startSigma));
      output.msg('v1s', '%s', msg);
    end

    output.msg('v0s', 'Finished. Total time: %.2f seconds.', toc(tStart));
end
