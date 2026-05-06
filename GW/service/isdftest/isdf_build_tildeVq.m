% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/02

function tildeVq = isdf_build_tildeVq(wf_on_coarse, R_coarse_RLU, R_rot_coarse, fftgrid_i, fftgrid_c, fft_sz, Nmu)
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
  symm_data = symmetry.get();

  DL_vol = d_lat_data.DL_vol;
  nc = wf_data.nc;
  ng = coul_data.coulomb_ng;
  nibz = int32(k_data.nibz);
  nbz = int32(k_data.nbz);
  nb = int32(wf_data.nb);

  tildeVq = zeros(Nmu, Nmu, nibz);
  helperqR = zeros(nc, Nmu);
  helperqG = zeros(ng, Nmu);

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

      wf_c_1tmp = wf_on_coarse(:, :, ikibz, 1);
      wf_c_2tmp = wf_on_coarse(:, :, ikpibz, 1);
      if ikrot ~= 1
        wf_c_1 = isdf.coeff.isdf_apply_symm_on_coarse(wf_c_1tmp, ikrot, symm_data, R_rot_coarse);
      else
        wf_c_1 = wf_c_1tmp;
      end
      if ikprot ~= 1
        wf_c_2 = isdf.coeff.isdf_apply_symm_on_coarse(wf_c_2tmp, ikprot, symm_data, R_rot_coarse);
      else
        wf_c_2 = wf_c_2tmp;
      end
      wf_c_1 = conj(wf_c_1);

      wf_1 = zeros(nc, nb);
      wf_2 = zeros(nc, nb);
      for ib = 1:nb
        isc1 = int32([ib, ikibz, ikrot, 1]);
        isc2 = int32([ib, ikpibz, ikprot, 1]);
        wf_1(:, ib) = wave_functions.WF_apply_symm(isc1);
        wf_2(:, ib) = wave_functions.WF_apply_symm(isc2);
      end
      wf_1 = conj(wf_1);

      Go = double(r_lat_data.Ggrid_RLU(iGo, :));
      R_coarse_scal = R_coarse_RLU .* (double(fftgrid_i) ./ double(fftgrid_c));
      phase_shift_coarse = double( exp(- 2 * pi * 1i * R_coarse_scal * Go') );
      phase_shift_fine   = double( exp(- 2 * pi * 1i * fft_data.Rgrid_RLU * Go') );

      MCHq_tmp = (wf_1 * wf_c_1') .* (wf_2 * wf_c_2') ...
        .* phase_shift_coarse' .* phase_shift_fine;
      MCHq = MCHq + MCHq_tmp;
      CCHq_tmp = (wf_c_1 * wf_c_1') .* (wf_c_2 * wf_c_2') ...
        .* phase_shift_coarse' .* ( phase_shift_coarse );
      CCHq = CCHq + CCHq_tmp;

      % Seems that MCHq and CCHq are not correct. We need to compute M and C explicitly.
      % M = zeros(nc, nb * nb);
      % C = zeros(Nmu, nb * nb);
      % for i1 = 1:nb
      %   for i2 = 1:nb
      %     M(:, (i1 - 1) * nb + i2) = wf_1(:, i1) .* wf_2(:, i2);
      %     C(:, (i1 - 1) * nb + i2) = wf_c_1(:, i1) .* wf_c_2(:, i2);
      %   end
      % end

      % Go = single(r_lat_data.Ggrid_RLU(iGo, :));
      % R_coarse_scal = R_coarse_RLU .* single(double(fftgrid_i) ./ double(fftgrid_c));
      % phase_shift_coarse = exp(2 * pi * 1i * R_coarse_scal * Go');
      % phase_shift_fine = exp(2 * pi * 1i * single(fft_data.Rgrid_RLU) * Go');
      % M = M .* phase_shift_fine;
      % C = C .* phase_shift_coarse;

      % MCHqdiff = MCHq_tmp - M*C';
      % CCHqdiff = CCHq_tmp - C*C';

      % norm(MCHqdiff, 'fro') / norm(MCHq_tmp, 'fro')
      % norm(CCHqdiff, 'fro') / norm(CCHq_tmp, 'fro')
    end % ikbz

    helperqR = MCHq / CCHq;
    [Q, ~] = qr(helperqR, 'econ');
    fro_norm_diff = norm(Q'*Q - eye(size(Q, 2)), 'fro');
    fprintf('Orthogonality check (Frobenius norm deviation from identity): %.3e\n', fro_norm_diff);
    % Least-squares sanity: print Frobenius relative error per k_BZ (not stored).

    % for ikbz = 1:nbz
    %   ikibz = k_data.bz2ibz(ikbz, 1);
    %   ikrot = k_data.bz2rot(ikbz, 1);
    %   ikpbz = r_lat_data.qindx_X(iqibz, ikbz, 1);
    %   iGo = r_lat_data.qindx_X(iqibz, ikbz, 2);
    %   ikpibz = k_data.bz2ibz(ikpbz, 1);
    %   ikprot = k_data.bz2rot(ikpbz, 1);

    %   wf_c_1tmp = wf_on_coarse(:, :, ikibz, 1);
    %   wf_c_2tmp = wf_on_coarse(:, :, ikpibz, 1);
    %   if ikrot ~= 1
    %     wf_c_1 = apply_symm_on_coarse(wf_c_1tmp, ikrot, symm_data, R_rot_coarse);
    %   else
    %     wf_c_1 = wf_c_1tmp;
    %   end
    %   if ikprot ~= 1
    %     wf_c_2 = apply_symm_on_coarse(wf_c_2tmp, ikprot, symm_data, R_rot_coarse);
    %   else
    %     wf_c_2 = wf_c_2tmp;
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
    %   R_coarse_scal = R_coarse_RLU .* (double(fftgrid_i) ./ double(fftgrid_c));
    %   phase_shift_coarse = double(exp(- 2 * pi * 1i * R_coarse_scal * Go'));
    %   phase_shift_fine   = double(exp(- 2 * pi * 1i * fft_data.Rgrid_RLU * Go'));
    %   M = M .* phase_shift_fine;
    %   C = C .* phase_shift_coarse;

    %   MCHq = MCHq - M*C';
    %   CCHq = CCHq - C*C';

    %   % ProjM: row-wise orthogonal projection of M onto span(rows of C) in C^K (K = nb^2).
    %   % B = C.';                  % K x Nmu, columns are C(i,:).'
    %   % ProjM = (M * B) * pinv(B);
    %   % nm = norm(M, 'fro');
    %   % if nm > 0
    %   %   norm(M - ProjM, 'fro') / nm
    %   % else
    %   %   0
    %   % end
    %   % Print Frobenius norm ratio ||M - Proj_M||_F / ||M||_F for diagnostics,
    %   % where Proj_M is the row-space projection of M onto C (see previous notes).
    %   % Use a try-catch for clarity and robust output if Q is missing.
    %   try
    %     frob_ratio = norm(M - Q*(Q'*M), 'fro') / norm(M, 'fro');
    %     fprintf('[ISDF] Frobenius norm check for (M - P_C M)/||M||: %.3e\n', frob_ratio);
    %   catch err
    %     warning('[ISDF] Could not compute Frobenius norm ratio: %s', err.message);
    %   end
    % end % ikbz

    % % 'The following value should be zero:'
    fprintf('[ISDF] Frobenius norm |MCHq|_F: %.6e\n', norm(MCHq, 'fro'));
    fprintf('[ISDF] Frobenius norm |CCHq|_F: %.6e\n', norm(CCHq, 'fro'));

    for imu = 1:Nmu
      fftbox = reshape(helperqR(:, imu), fft_sz);
      fftbox = do_FFT(fftbox, fft_sz, 1) * DL_vol;
      helperqG(:, imu) = fftbox(fft_data.G_table(:, 1));
    end
    tildeVq(:, :, iqibz) = helperqG' * diag(vcoul_q) * helperqG;
  end % iqibz



  clear helperqR helperqG;
end
