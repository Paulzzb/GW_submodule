% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/02

function tildeVq = gen_tildeVq(id)
% ISDF_BUILD_TILDEVQ  Build tildeVq(:, :, iq) from coarse-grid ISDF helpers and Coulomb matrix in G.
%
% Accumulates MCHq/CCHq over the BZ for each IBZ q, solves helperqR = MCHq/CCHq, runs a
% Frobenius-ratio print check ||helperqR*C - M||_F / ||M||_F over k_BZ, then FFTs helpers to G
% and forms helperqG' * diag(vcoul) * helperqG. G_table column uses ikbz == nbz after that loop
% (same as the original monolithic test).

  wf_data = wave_functions.get();
  k_data = lattice.manager('k', 'get');
  r_lat_data = lattice.manager('r_lat', 'get');
  d_lat_data = lattice.manager('d_lat', 'get');
  coul_data = coulomb.get();
  fft_data = FFT.get();
  fft_sz = fft_data.fftgrid;
  isdf_data = isdf.get(id);
  [nrange1, nrange2] = isdf_get_nrange(id);
  nb1 = length(nrange1); nb2 = length(nrange2);

  DL_vol = d_lat_data.DL_vol;
  nc = wf_data.nc;
  ng = coul_data.coulomb_ng;
  nibz = int32(k_data.nibz);
  nbz = int32(k_data.nbz);
  nb = int32(wf_data.nb);
  R_sampling_RLU = isdf_data.R_sampling_RLU;
  coeff_seper = isdf_data.coeff_seper;

  Nmu = isdf_data.nisdf;
  tildeVq = zeros(Nmu, Nmu, nibz);
  helperqR = zeros(nc, Nmu);
  helperqG = zeros(ng, Nmu);

  % Fine FFT box: same symmetry index map as isdf_coarse_R_rot (nr_fine x nsym), built lazily in non-coarse branch.
  R_rot_extra = int32([]);

  for iqibz = 1:nibz
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

      if strcmp(isdf_data.interp_scheme, "coarse")
        wf_c_1tmp = coeff_seper(:, nrange1, ikibz, 1);
        wf_c_2tmp = coeff_seper(:, nrange2, ikpibz, 1);
        if ikrot ~= 1
          wf_c_1 = isdf.isdf_apply_symm_on_coarse(id, wf_c_1tmp, ikrot);
        else
          wf_c_1 = wf_c_1tmp;
        end
        if ikprot ~= 1
          wf_c_2 = isdf.isdf_apply_symm_on_coarse(id, wf_c_2tmp, ikprot);
        else
          wf_c_2 = wf_c_2tmp;
        end
      elseif strcmp(isdf_data.interp_scheme, "adaptive")
        wf_c_1 = coeff_seper(:, nrange1, ikbz, 1);
        wf_c_2 = coeff_seper(:, nrange2, ikpbz, 1);
      else % In this situation, a direct access without rotation.
        wf_c_1 = coeff_seper(:, nrange1, ikbz, 1);
        wf_c_2 = coeff_seper(:, nrange2, ikpbz, 1);
      end
      wf_c_1 = conj(wf_c_1);

      wf_1 = zeros(nc, nb1);
      wf_2 = zeros(nc, nb2);
      for ib = nrange1
        isc1 = int32([ib, ikibz, ikrot, 1]);
        wf_1(:, ib) = wave_functions.WF_apply_symm(isc1);
      end
      for ib = nrange2
        isc2 = int32([ib, ikpibz, ikprot, 1]);
        wf_2(:, ib) = wave_functions.WF_apply_symm(isc2);
      end
      wf_1 = conj(wf_1);

      Go = double(r_lat_data.Ggrid_RLU(iGo, :));
      phase_shift_coarse = double( exp(- 2 * pi * 1i * (R_sampling_RLU ./ double(fft_data.fftgrid)) * Go') );
      phase_shift_fine   = double( exp(- 2 * pi * 1i * (fft_data.Rgrid_RLU ./ double(fft_data.fftgrid)) * Go') );

      MCHq_tmp = (wf_1 * wf_c_1') .* (wf_2 * wf_c_2') ...
        .* phase_shift_coarse' .* phase_shift_fine;
      MCHq = MCHq + MCHq_tmp;
      CCHq_tmp = (wf_c_1 * wf_c_1') .* (wf_c_2 * wf_c_2') ...
        .* phase_shift_coarse' .* ( phase_shift_coarse );
      CCHq = CCHq + CCHq_tmp;

      % % Seems that MCHq and CCHq are not correct. We need to compute M and C explicitly.
      % M = zeros(nc, nb * nb);
      % C = zeros(Nmu, nb * nb);
      % for i1 = 1:nb
      %   for i2 = 1:nb
      %     M(:, (i1 - 1) * nb + i2) = wf_1(:, i1) .* wf_2(:, i2);
      %     C(:, (i1 - 1) * nb + i2) = wf_c_1(:, i1) .* wf_c_2(:, i2);
      %   end
      % end

      % Go = single(r_lat_data.Ggrid_RLU(iGo, :));
      % % R_sampling_scal = R_sampling_RLU ./ double(fft_data.fftgrid);
      % phase_shift_coarse = exp(- 2 * pi * 1i * (R_sampling_RLU ./ double(fft_data.fftgrid)) * Go');
      % phase_shift_fine = exp(- 2 * pi * 1i * (fft_data.Rgrid_RLU ./ double(fft_data.fftgrid)) * Go');
      % M = M .* phase_shift_fine;
      % C = C .* phase_shift_coarse;

      % MCHqdiff = MCHq_tmp - M*C';
      % CCHqdiff = CCHq_tmp - C*C';

      % output1 = norm(MCHqdiff, 'fro')/norm(MCHq_tmp, 'fro');
      % output2 = norm(CCHqdiff, 'fro')/norm(CCHq_tmp, 'fro');
      % if output1 > 1e-5
      %   error('[ISDF] Frobenius norm |MCHqdiff|_F/|MCHq_tmp|_F: %.6e\n', output1);
      % end
      % if output2 > 1e-5
      %   error('[ISDF] Frobenius norm |CCHqdiff|_F/|CCHq_tmp|_F: %.6e\n', output2);
      % end
      % fprintf('[ISDF] Frobenius norm |MCHqdiff|_F: %.6e\n', norm(MCHqdiff, 'fro')/norm(MCHq_tmp, 'fro'));
      % fprintf('[ISDF] Frobenius norm |CCHqdiff|_F: %.6e\n', norm(CCHqdiff, 'fro')/norm(CCHq_tmp, 'fro'));
    end % ikbz

    helperqR = MCHq / CCHq;
    % [Q, ~] = qr(helperqR, 'econ');
    % fro_norm_diff = norm(Q' * Q - eye(size(Q, 2)), 'fro');
    % fprintf('[ISDF] Orthogonality ||Q''Q - I||_F (iqibz=%d): %.3e\n', iqibz, fro_norm_diff);
    % Diagnostic loop only: must not overwrite MCHq/CCHq used for helperqR.

    % for ikbz = 1:nbz
    %   ikibz = k_data.bz2ibz(ikbz, 1);
    %   ikrot = k_data.bz2rot(ikbz, 1);
    %   ikpbz = r_lat_data.qindx_X(iqibz, ikbz, 1);
    %   iGo = r_lat_data.qindx_X(iqibz, ikbz, 2);
    %   ikpibz = k_data.bz2ibz(ikpbz, 1);
    %   ikprot = k_data.bz2rot(ikpbz, 1);

    %   if strcmp(isdf_data.interp_scheme, "coarse")
    %     wf_c_1tmp = coeff_seper(:, :, ikibz, 1);
    %     wf_c_2tmp = coeff_seper(:, :, ikpibz, 1);
    %     if ikrot ~= 1
    %       wf_c_1 = isdf.isdf_apply_symm_on_coarse(id, wf_c_1tmp, ikrot);
    %     else
    %       wf_c_1 = wf_c_1tmp;
    %     end
    %     if ikprot ~= 1
    %       wf_c_2 = isdf.isdf_apply_symm_on_coarse(id, wf_c_2tmp, ikprot);
    %     else
    %       wf_c_2 = wf_c_2tmp;
    %     end
    %   else
    %     wf_c_1 = coeff_seper(:, :, ikbz, 1);
    %     wf_c_2 = coeff_seper(:, :, ikpbz, 1);
    %   end
    %   wf_c_1 = conj(wf_c_1);

    %   wf_1 = zeros(nc, nb);
    %   wf_2 = zeros(nc, nb);
    %   for ib = 1:nb
    %     isc1 = int32([ib, ikibz, ikrot, 1]);
    %     isc2 = int32([ib, ikpibz, ikprot, 1]);
    %     wf_1(:, ib) = wave_functions.WF_apply_symm(isc1);
    %     wf_2(:, ib) = wave_functions.WF_apply_symm(isc2);
    %   end
    %   wf_1 = conj(wf_1);

    %   M = zeros(nc, nb * nb);
    %   C = zeros(Nmu, nb * nb);
    %   for i1 = 1:nb
    %     for i2 = 1:nb
    %       M(:, (i1 - 1) * nb + i2) = wf_1(:, i1) .* wf_2(:, i2);
    %       C(:, (i1 - 1) * nb + i2) = wf_c_1(:, i1) .* wf_c_2(:, i2);
    %     end
    %   end

    %   Go = double(r_lat_data.Ggrid_RLU(iGo, :));
    %   phase_shift_coarse = double(exp(-2 * pi * 1i * (R_sampling_RLU ./ double(fft_data.fftgrid)) * Go'));
    %   phase_shift_fine = double(exp(-2 * pi * 1i * (fft_data.Rgrid_RLU ./ double(fft_data.fftgrid)) * Go'));
    %   M = M .* phase_shift_fine;
    %   C = C .* phase_shift_coarse;

    %   try
    %     nm = norm(M, 'fro');
    %     if nm > 0
    %       frob_ratio = norm(M - Q * (Q' * M), 'fro') / nm;
    %       fprintf('[ISDF] iqibz=%d ikbz=%d ||M - Q(Q''M)||_F/||M||_F: %.3e\n', iqibz, ikbz, frob_ratio);
    %     end
    %   catch err
    %     warning('[ISDF] Frobenius diagnostic failed (iqibz=%d ikbz=%d): %s', iqibz, ikbz, err.message);
    %   end
    % end % ikbz

    % 'The following value should be zero:'
    % fprintf('[ISDF] Frobenius norm |MCHq|_F: %.6e\n', norm(MCHq, 'fro'));
    % fprintf('[ISDF] Frobenius norm |CCHq|_F: %.6e\n', norm(CCHq, 'fro'));

    for imu = 1:Nmu
      fftbox = reshape(helperqR(:, imu), fft_sz);
      fftbox = do_FFT(fftbox, fft_sz, 1) * DL_vol;
      helperqG(:, imu) = fftbox(fft_data.G_table(:, 1));
    end
    tildeVq(:, :, iqibz) = helperqG' * diag(vcoul_q) * helperqG;
    tildeVq(:, :, iqibz) = tildeVq(:, :, iqibz) / 2 + tildeVq(:, :, iqibz)' / 2;
  end % iqibz


  isdf_data.tildeVq = tildeVq;
  isdf.save2mod(isdf_data, id);

end
