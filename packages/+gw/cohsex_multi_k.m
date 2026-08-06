%
% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/27 ZZ
function [Esx_x, Ecoh] = cohsex_multi_k(config)
%GW_COHSEX_MULTI_K  Static COHSEX self-energy diagonal (ISDF or dense G-space).
%
%   [Esx_x, Ecoh] = gw.cohsex_multi_k(GWinfo, config)
%
% ISDF path uses service/+isdf (isdf.get, coeff_seper, gen_tildeVq, cohsex_vcVnn).
% Legacy ISDFDB + isdf_sub is not used here.
%
% ISDF prerequisites:
%   - Relay/service initialized (wave_functions, system, lattice, FFT, coulomb, isdf pool).
%   - Assigned slots desc ''vc'' and ''nn'' with coeff_seper, or config.ISDF.id_vc / id_nn.
%   - Single IBZ k index ikibz = 1 and ispin = 1 for this prototype (multi-k extension later).

  default_Constant = constant_map();
  nameConstants = fieldnames(default_Constant);
  for i = 1:numel(nameConstants)
    eval(sprintf('%s = %.16f;', nameConstants{i}, default_Constant.(nameConstants{i})));
  end

  % from config
  nbmin = config.SYSTEM.energy_band_index_min;
  nbmax = config.SYSTEM.energy_band_index_max;
  nband = nbmax - nbmin + 1;
  
  %
  system_data = system.get();
  k_data = lattice.manager('k', 'get');
  q_data = lattice.manager('q', 'get');
  r_lat_data = lattice.manager('r_lat', 'get');
  d_lat_data = lattice.manager('d_lat', 'get');
  coulomb_data = coulomb.get();
  % 
  nkibz = k_data.nibz;
  nqibz = q_data.nibz;
  nkbz = k_data.nbz;
  nqbz = q_data.nbz;

  ev = system_data.Eo * ry2ev;

  nspin = 1;
  ispin = 1;
  msg = sprintf('Multi-spin is not supported yet.\n');
  output.warn('%s', msg);


  if (config.ISDF.isisdf)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For each q \in IBZ, calculate \Lambda_q, K_q, W_q
    %
    exact_CH = false;
    % exact_ch_mode = "legacy";
    exact_ch_mode = "dv_once_half";
    exact_ch_debug = false;
    if isfield(config, 'COHSEX') && isfield(config.COHSEX, 'exact_ch')
      exact_CH = logical(config.COHSEX.exact_ch);
    end
    if isfield(config, 'COHSEX') && isfield(config.COHSEX, 'exact_ch_mode')
      exact_ch_mode = string(config.COHSEX.exact_ch_mode);
    end
    if isfield(config, 'COHSEX') && isfield(config.COHSEX, 'exact_ch_debug')
      exact_ch_debug = logical(config.COHSEX.exact_ch_debug);
    end
    if ~exact_CH
      [id_vc, id_vn, id_nn] = isdf.cohsex_resolve_ids(config);
    else
      [id_vc, id_vn, ~] = isdf.cohsex_resolve_ids(config);
    end
    
    if isempty(id_vn) && isempty(id_vc)
      output.err( ...
        ['Static COHSEX requires both vc and vn ISDF slots. ' ...
         'Please enable config.ISDF.compute_vc/compute_vn or provide assigned ids.']);
    end
    if ~exact_CH && isempty(id_nn)
      output.err( ...
        ['Non-exact CH branch requires nn ISDF slot. ' ...
         'Enable config.ISDF.compute_nn=true or set COHSEX.exact_ch=true.']);
    end

    isdf.set_nrange(id_vc, config.SYSTEM);
    vc_data = isdf.get(id_vc);
    vn_data = isdf.get(id_vn);
    nrangev = double(vc_data.nrange1);
    nrangec = double(vc_data.nrange2);

    Nisdf_vc = vc_data.nisdf;
    Nisdf_vn = vn_data.nisdf;
    nn_data = [];
    Nisdf_nn = 0;


    if ~exact_CH
      if ~isempty(id_nn)
        nn_data = isdf.get(id_nn);
        Nisdf_nn = nn_data.nisdf;
      end
      Sigma_sex_x = zeros(nband, nkibz);
      Sigma_coh = zeros(nband, nkibz);
      iqibz_old = 0;
      for iqbz = 1:nqbz
        iqibz = q_data.bz2ibz(iqbz);
        iqrot = q_data.bz2rot(iqbz);
        if iqibz ~= iqibz_old
          % tildeW_q = isdf.gen_tildeWq(id_vc, iqibz);
          Kq_ISDF = isdf.gen_Kq(id_vc, iqibz);
          tildeWq_vn = isdf.gen_tildeWq(id_vc, iqibz, Kq_ISDF, id_vn, true); % for SEX
          tildeWq_nn = isdf.gen_tildeWq(id_vc, iqibz, Kq_ISDF, id_nn, true); % for COH
          % verify
          ng = size(vc_data.helperqG, 1);
          if size(vn_data.helperqG, 1) ~= ng || size(nn_data.helperqG, 1) ~= ng
            msg = sprintf('verifyW: inconsistent ng among vc/vn/nn helperqG, skip verification at iqibz=%d.\n', iqibz);
            output.warn('%s', msg);
          else
            d_lat_data = lattice.manager('d_lat', 'get');
            vol = double(d_lat_data.DL_vol);
            ev_ry = double(system_data.Eo);
            spin_id = 1;
            scal = 4.0;

            vcoul_q = double(coulomb.get().vcoul(:, iqibz));
            if iqibz == 1
              vcoul_q(1) = double(coulomb.get().vcoul0);
            end
            Dcoul = spdiags(vcoul_q(:), 0, ng, ng);

            chiq_G = zeros(ng, ng);
            nrangev_row = reshape(nrangev, 1, []);
            nrangec_row = reshape(nrangec, 1, []);
            ncb = numel(nrangec_row);
            for ikbz = 1:nkbz
              ikibz_k = k_data.bz2ibz(ikbz);
              ikrot_k = k_data.bz2rot(ikbz);
              ikq_bz = r_lat_data.qindx_X(iqibz, ikbz, 1);
              iGo_x = r_lat_data.qindx_X(iqibz, ikbz, 2);
              ikq_ibz = k_data.bz2ibz(ikq_bz);
              ikq_rot = k_data.bz2rot(ikq_bz);

              f_c = double(system_data.f(nrangec_row, ikq_ibz, spin_id));
              e_c = ev_ry(nrangec_row, ikq_ibz, spin_id);
              for iv = nrangev_row
                % Build all Mgvc(:, jc) first, then do a double BLAS-3 update.
                Mgvc_blk = zeros(ng, ncb);
                for jc_id = 1:ncb
                  jc = nrangec_row(jc_id);
                  pchk = struct();
                  pchk.is = [iv, ikibz_k, ikrot_k, spin_id];
                  pchk.os = [jc, ikq_ibz, ikq_rot, spin_id];
                  pchk.qs = [iGo_x, iqibz, 1];
                  Mgvc_blk(:, jc_id) = double(SCATTER_Bamp(pchk));
                end

                f_v = double(system_data.f(iv, ikibz_k, spin_id));
                e_v = ev_ry(iv, ikibz_k, spin_id);
                occ = f_v - f_c;
                den = e_v - e_c;
                valid = (abs(occ) >= 1e-8) & (abs(den) >= 1e-12);
                if ~any(valid)
                  continue;
                end

                coeff = occ(valid) ./ den(valid);
                Mgvc_valid = Mgvc_blk(:, valid);
                Mgvc_weighted = Mgvc_valid .* reshape(coeff, 1, []);
                chiq_G = chiq_G + scal * (Mgvc_weighted * Mgvc_valid');
              end
            end

            inveps = eye(ng) - (Dcoul * chiq_G);
            W_dense_mul = full(inveps * Dcoul);
            W_dense_solve = full(inveps \ Dcoul);
            W_v = W_dense_solve - Dcoul;

            helperqG_vc = double(vc_data.helperqG(:, :, iqibz));
            W_isdf =  - diag(vcoul_q) * helperqG_vc * double(inv(Kq_ISDF)) * helperqG_vc' * diag(vcoul_q);

            denom_solve = max(norm(W_v, 'fro'), eps);
            diff_solve = norm(W_v - W_isdf, 'fro');

            msg = sprintf(['[verifyW] iqibz=%d ng=%d Nisdf(vc/vn/nn)=(%d,%d,%d)\n', ...
                     '  solve: ||W-Wisdf||_F=%.6e (rel=%.6e)\n'], ...
                    iqibz, ng, Nisdf_vc, Nisdf_vn, Nisdf_nn, ...
                    diff_solve, diff_solve / denom_solve);
            output.msg('v2s', '%s', msg);
          end
        end
        iqibz_old = iqibz;
        %
        % Calculate Sigma^{SEX_X}_{n\kk} and Sigma^{COH}_{n\kk}, notice summation
        % over q is outermost.
        % <n\kk|Sigma^{SEX_X}|n\kk>
        % = - \sum_{n_v \in \rangev} \sum_{qbz\in\bz} \sum_{\mu, \nu}  f_{n_v\kk-\qq_bz}
        %   \conj{\rho_{n n_v}(\kk, \qq_ibz, S_{qbz}\rr_mu)}
        %   * \tildeW_q^{vn}(\mu, \nu)
        %   * \rho_{n n_v}(\kk, \qq_ibz, S_{qbz}\rr_mu),
        % <n\kk|Sigma^{COH}|n\kk>
        % = 0.5 \sum_{n \in \range} \sum_{qbz\in\bz} \sum_{\mu, \nu}
        %   \conj{\rho_{n n_v}(\kk, \qq_ibz, S_{qbz}\rr_mu)}
        %   * \tildeW_q^{nn}(\mu, \nu)
        %   * \rho_{n n_v}(\kk, \qq_ibz, S_{qbz}\rr_mu),
        % where S_qbz * \qq_ibz = \qq_bz, tildeWq^{vn} and tildeWq^{nn} are obtained from
        % isdf.gen_tildeWq with outer id as 'id_vn' and 'id_nn', respectively.
        % In this non-exact branch, COH uses rho' * tildeWq_nn * rho in ISDF space.
        for ispin = 1:nspin
          for ikibz = 1:nkibz
            ikrot = 1;
            ikqbz = r_lat_data.qindx_S(ikibz,iqbz,1);
            ikqibz = k_data.bz2ibz(ikqbz);
            ikqrot = k_data.bz2rot(ikqbz);
            iGo = r_lat_data.qindx_S(ikibz,iqbz,2);
            param = [];
            param.qs = [iGo, iqibz, iqrot];
            for indib = 1:nband
              ib = nbmin + indib - 1;
              param.is = [ib, ikibz, ikrot, ispin];
              for idob = 1:length(nrangev)
                ob = nrangev(idob);
                param.os = [ob, ikqibz, ikqrot, ispin];
                rho_left_vn = isdf.get_rho_xalpha(id_vn, param);
                tmp1 = rho_left_vn' * tildeWq_vn * rho_left_vn;
                tmp1 = real(tmp1);
                Sigma_sex_x(indib, ikibz) = Sigma_sex_x(indib, ikibz) + tmp1;
                % 
                rho_left_nn = isdf.get_rho_xalpha(id_nn, param);
                tmp2 = rho_left_nn' * tildeWq_nn * rho_left_nn;
                tmp2 = real(tmp2);
                Sigma_coh(indib, ikibz) = Sigma_coh(indib, ikibz) - 0.5 * tmp2;
              end
              for idob = 1:length(nrangec)
                ob = nrangec(idob);
                param.os = [ob, ikqibz, ikqrot, ispin];
                rho_left_nn = isdf.get_rho_xalpha(id_nn, param);
                tmp2 = rho_left_nn' * tildeWq_nn * rho_left_nn;
                tmp2 = real(tmp2);
                Sigma_coh(indib, ikibz) = Sigma_coh(indib, ikibz) - 0.5 * tmp2;
              end 
            end
          end
        end
      end % iqibz
    else % exact_CH
      Sigma_sex_x = zeros(nband, nkibz);
      Sigma_coh = zeros(nband, nkibz);
      exact_ch_dbg_print_count = 0;
      iqibz_old = 0;
      for iqbz = 1:nqbz
        iqibz = q_data.bz2ibz(iqbz);
        iqrot = q_data.bz2rot(iqbz);
        if iqibz ~= iqibz_old
          vcoul_q = coulomb_data.vcoul(:, iqibz);
          if iqibz == 1
            vcoul_q(1) = coulomb_data.vcoul0;
          end
          Kq_ISDF = isdf.gen_Kq(id_vc, iqibz);
          tildeWq_vn = isdf.gen_tildeWq(id_vc, iqibz, Kq_ISDF, id_vn); % for SEX
          helperqG_vc = double(vc_data.helperqG(:, :, iqibz));
          % Recover helperqR from helperqG (inverse of gen_tildeVq's R->G mapping).
          fft_data = FFT.get();
          fft_sz = double(fft_data.fftgrid(:).');
          DL_vol = double(lattice.manager('d_lat', 'get').DL_vol);
          nmu_vc = size(helperqG_vc, 2);
          nr_fine = prod(fft_sz);
          helperqR_vc = zeros(nr_fine, nmu_vc);
          for imu = 1:nmu_vc
            fftbox = zeros(fft_sz);
            fftbox(fft_data.G_table(:, 1)) = helperqG_vc(:, imu) .* vcoul_q;
            fftbox = do_FFT(fftbox, fft_sz, -1) / DL_vol;
            % fftbox = do_FFT(fftbox, fft_sz, -1) * r_lat_data.RL_vol;
            helperqR_vc(:, imu) = fftbox(:);
          end
          % Pick 1-2 imu channels, map helperqR_vn back to G with the same
          % convention as gen_tildeVq, and require it matches inverse target.
          imu_check = unique([1, min(2, nmu_vc)]);
          fft_back_tol = 1e-8;
          for iid = 1:numel(imu_check)
            imu_chk = imu_check(iid);
            fftbox_chk = reshape(helperqR_vc(:, imu_chk), fft_sz);
            fftbox_chk = do_FFT(fftbox_chk, fft_sz, 1) * DL_vol;
            helperqG_back = fftbox_chk(fft_data.G_table(:, 1));
            helperqG_target_tilde = helperqG_vc(:, imu_chk) .* vcoul_q;
            rel_tilde = norm(helperqG_back - helperqG_target_tilde) / max(norm(helperqG_target_tilde), eps);
            if rel_tilde > fft_back_tol
              output.err( ...
                ['helperqR_vc back-to-G check failed at iqibz=%d, imu=%d: ' ...
                 'rel_err=%.3e exceeds tol=%.1e.'], ...
                iqibz, imu_chk, rel_tilde, fft_back_tol);
            end
            if exact_ch_debug
              helperqG_target_raw = helperqG_vc(:, imu_chk);
              rel_raw = norm(helperqG_back - helperqG_target_raw) / max(norm(helperqG_target_raw), eps);
              msg = sprintf(['[exact_CH debug] helperqR<->helperqG check iqibz=%d imu=%d ' ...
                       'rel(back,helperqG*vcoul)=%.6e rel(back,helperqG)=%.6e\n'], ...
                      iqibz, imu_chk, rel_tilde, rel_raw);
              output.msg('v2s', '%s', msg);
            end
          end
          invK = inv(Kq_ISDF);
          Wqrr = zeros(nr_fine, 1);
          block_size = 64;
          % This part is for the delta function factor
          coeff_delta = double(fft_data.nr);
          % This part is for the operator factor.
          % A_R = F^H * A_G * F / nr, while in the previous part, using helper functions
          % we calculate W_R = (F^H / vol) * W_G * (F / vol)
          coeff_delta = coeff_delta * (d_lat_data.DL_vol^2 / double(fft_data.nr));
          for iblock = 1:block_size:nr_fine
            iend = min(iblock + block_size - 1, nr_fine);
            Hblk = helperqR_vc(iblock:iend, :);      % nblk x nmu
            tmp_blk = invK * Hblk.';                 % nmu x nblk
            Wqrr(iblock:iend) = 0.5 * coeff_delta * conj(sum(conj(Hblk.') .* tmp_blk, 1)).';
            % Wqrr(iblock:iend) = 0.5 * conj(sum(conj(Hblk.') .* tmp_blk, 1)).';
          end
          Wqrr = 0.5 * (Wqrr + conj(Wqrr));
        end
        iqibz_old = iqibz;
        %
        % Calculate Sigma^{SEX_X}_{n\kk} and Sigma^{COH}_{n\kk}, where COH uses exact CH.
        %     Sigma^{exact_CH}(r, r') = \delta(r-r') * W(q; r, r') = \delta(r-r') * W(q; r, r).
        % notice that in ISDF
        %     W(q; G, G') = \conj{p^q_\mu(G)} * v(G) * K^{-1}_{\mu, \nu} * p^q_\nu(G') * v(G').
        % We let \tildep^q_\mu(G) = p^q_\mu(G) * v(q; G), then
        %     W(q; G, G') = \conj{\tildep^q_\mu(G)} * K^{-1}_{\mu, \nu} * \tildep^q_\nu(G').
        % and
        %     W(q; r, r) = \conj{\tildep^q_\mu(r)} * K^{-1}_{\mu, \nu} * \tildep^q_\nu(r).
        %
% From BGW
% < n k | \Sigma_{CH} (r, r`; 0) | m k > =
% \frac{1}{2} \sum_{q G G`}
% < n k | e^{i (G - G`) \cdot r} | m k >
% [\eps_{G G`}^{-1} (q; 0) - \delta_{G G`}] v (q + G`)
% 关于 G 的可以不要（卷积核在实空间是直接乘法）
% 但是 0.5 的系数似乎是必须的
        % <n\kk|Sigma^{COH}|n\kk>
        % = \sum_{qbz\in\bz} \sum_{\mu, \nu} |\psi_{n\kk}(\rr)|^2 * W(q; r, r')
        % = \sum_{qbz\in\bz} \sum_{r} |\psi_{n\kk}(\rr)|^2 * W(q; r, r)
        % 
        % where S_qbz * \qq_ibz = \qq_bz, tildeWq^{vn} and tildeWq^{nn} are obtained from
        % isdf.gen_tildeWq with outer id as 'id_vn' and 'id_nn', respectively.
        for ispin = 1:nspin
          for ikibz = 1:nkibz
            ikrot = 1;
            ikqbz = r_lat_data.qindx_S(ikibz,iqbz,1);
            ikqibz = k_data.bz2ibz(ikqbz);
            ikqrot = k_data.bz2rot(ikqbz);
            iGo = r_lat_data.qindx_S(ikibz,iqbz,2);
            param = [];
            param.qs = [iGo, iqibz, iqrot];
            for indib = 1:nband
              ib = nbmin + indib - 1;
              param.is = [ib, ikibz, ikrot, ispin];
              for idob = 1:length(nrangev)
                ob = nrangev(idob);
                param.os = [ob, ikqibz, ikqrot, ispin];
                rho_left_vn = isdf.get_rho_xalpha(id_vn, param);
                tmp1 = rho_left_vn' * tildeWq_vn * rho_left_vn;
                tmp1 = real(tmp1);
                Sigma_sex_x(indib, ikibz) = Sigma_sex_x(indib, ikibz) - tmp1;
              end
              % Calculate exact CH part
              isc = [ib, ikibz, ikrot, ispin];
              wf_nk = wave_functions.WF_apply_symm(isc);
              % tmp = sum(Wqrr .* abs(wf_nk).^2) * d_lat_data.DL_vol / nr_fine;
              rho_nk = abs(wf_nk).^2;
              dv = d_lat_data.DL_vol / double(fft_data.nr);
              coh_legacy = sum((dv * Wqrr) .* rho_nk) * dv;
              coh_dv_once = sum(Wqrr .* rho_nk) * dv;
              coh_dv_once_half = coh_dv_once;
              switch lower(char(exact_ch_mode))
                case 'legacy'
                  tmp = coh_legacy;
                case 'dv_once'
                  tmp = coh_dv_once;
                case 'dv_once_half'
                  tmp = coh_dv_once_half;
                otherwise
                  output.err( ...
                    'Unknown COHSEX.exact_ch_mode = "%s". Use legacy | dv_once | dv_once_half.', ...
                    char(exact_ch_mode));
              end
              if exact_ch_debug && exact_ch_dbg_print_count < 12
                exact_ch_dbg_print_count = exact_ch_dbg_print_count + 1;
                norm_plain = sum(rho_nk);
                norm_dv = sum(rho_nk) * dv;
                msg = sprintf(['[exact_CH debug] mode=%s iqibz=%d ikibz=%d ib=%d ' ...
                         'norm_plain=%.6e norm_dv=%.6e coh_legacy=%.6e ' ...
                         'coh_dv_once=%.6e coh_dv_once_half=%.6e ratio(legacy/half)=%.6e\n'], ...
                        char(exact_ch_mode), iqibz, ikibz, ib, ...
                        norm_plain, norm_dv, coh_legacy, coh_dv_once, ...
                        coh_dv_once_half, coh_legacy / max(abs(coh_dv_once_half), eps));
                output.msg('v2s', '%s', msg);
              end
              Sigma_coh(indib, ikibz) = Sigma_coh(indib, ikibz) + tmp;
            end
          end
        end
      end % iqibz
    end % exact_CH


    Esx_x = Sigma_sex_x;
    Ecoh = Sigma_coh;
  else
    % Dense no-ISDF prototype (double-k/double-q path for validation).
    if nkibz ~= 1 || nqbz ~= 1
      msg = sprintf( ...
        ['noISDF path currently ignores k-q loops and uses only (ikibz,iqbz)=(1,1). ', ...
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
        'Invalid summation bands: nv=%d, nsum=%d (need nsum>nv).', nv, nsum);
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
    msg = sprintf('[noISDF] timeforW = %.4f sec.\n', toc(startforW));
    output.msg('v0s', '%s', msg);

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
    msg = sprintf('[noISDF] timeforEsx_xExEch = %.4f sec.\n', toc(startSigma));
    output.msg('v0s', '%s', msg);
  end

  % Esx_x = real(diag(Esx_x));
  % Ecoh = real(diag(Ecoh));
end
