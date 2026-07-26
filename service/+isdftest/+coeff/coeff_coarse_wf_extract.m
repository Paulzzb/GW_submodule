% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function wf_on_coarse = coeff_coarse_wf_extract(fftgrid_c, R_coarse_RLU, id)
% Wavefunctions on coarse RLU sites: fine real-space WF -> FFT -> G_table -> dense map to coarse RLU.
%
%   Optional coincident-site check: isdftest.debug.init_from_config (see ISDF_debug.md),
%   tag coeff/coeff_coarse_wf_extract.

  fft_data = FFT.get();
  fftgrid_i = int32(fft_data.fftgrid(:).');
  fft_sz = double(fftgrid_i);
  wf_data = wave_functions.get();
  k_data = lattice.manager('k', 'get');
  r_lat_data = lattice.manager('r_lat', 'get');
  d_lat_data = lattice.manager('d_lat', 'get');

  DL_vol = d_lat_data.DL_vol;
  Nmu = int32(prod(double(fftgrid_c(:)')));

  scale_vec = double(fftgrid_i) ./ double(fftgrid_c(:)');
  R_scaled = double(R_coarse_RLU) .* scale_vec ./ double(fftgrid_i);
  phase_arg = 2 * pi * (R_scaled * double(r_lat_data.Ggrid_RLU)');
  scal = 1 ./ ( prod(fft_sz) * double(DL_vol) );
  % scal = prod(fft_sz) ./ double(DL_vol);
  % iFFTmat = exp(1i * phase_arg) * scal / prod(fft_sz);
  iFFTmat = exp(1i * phase_arg) * scal;

  nibz = int32(k_data.nibz);
  nb = int32(wf_data.nb);
  nspin = int32(wf_data.nspin);

  % Fine-RLU coincidence: R_fine = R_coarse_RLU .* (fftgrid_i ./ fftgrid_c) (same as gen_coeff_coarse scaling).
  % When R_fine is integer mod fine torus, coarse iFFT path must match direct wf_b(lin_fine).
  tol_coincide = 1e-4;
  wf_coincide_tol = 1e-3;
  Gi = double(fftgrid_i(:)).';
  Gc = double(fftgrid_c(:)).';
  Rf_all = double(R_coarse_RLU) .* (Gi ./ Gc);
  on_fine = max(abs(Rf_all - round(Rf_all)), [], 2) < tol_coincide;
  lin_on_fine = zeros(double(Nmu), 1, 'int32');

  if any(on_fine)
    for imu = 1:double(Nmu)
      if ~on_fine(imu)
        continue;
      end
      Rfm = mod(round(Rf_all(imu, :)) + Gi, Gi);
      li = 1 + Rfm(1) + Gi(1) * Rfm(2) + Gi(1) * Gi(2) * Rfm(3);
      lin_on_fine(imu) = int32(round(li));
      if lin_on_fine(imu) < 1 || lin_on_fine(imu) > fft_data.nr
        error('isdftest.coeff:coeff_coarse_wf_extract:LinFine', ...
          'Coincident coarse site imu=%d maps to invalid fine lin=%d (nr=%d).', imu, lin_on_fine(imu), fft_data.nr);
      end
    end
  end

  wf_on_coarse = zeros(Nmu, nb, nibz, nspin);

  for ik_ibz = 1:nibz
    ik_rot = 1;
    for ispin = 1:nspin
      for ib = 1:nb
        isc = int32([ib, ik_ibz, ik_rot, ispin]);
        wf_b = wave_functions.WF_apply_symm(isc);
        fftbox = double( reshape(wf_b, fft_sz) );
        fftbox = do_FFT(fftbox, fft_sz, -1) * double(DL_vol);
        psig = fftbox(fft_data.G_table(:, 1));
        wf_b_coarse = (iFFTmat * psig).';
        wf_col = wf_b_coarse(:);
        wf_on_coarse(:, ib, ik_ibz, ispin) = wf_col;
        if isdftest.debug.on('coeff/coeff_coarse_wf_extract')
          if any(on_fine)
            ix = find(on_fine);
            wbf = wf_b(:, 1);
            dmx = max(abs(wf_col(ix) - wbf(double(lin_on_fine(ix)))));
            msg = sprintf(['Coarse iFFT extract differs from fine WF at coincident sites: ', ...
              'max|diff|=%g (ib=%d ik_ibz=%d ispin=%d).'], dmx, ib, ik_ibz, ispin);
            isdftest.debug.react(dmx > wf_coincide_tol, msg, 'coeff_coarse_wf_extract_coincide');
          end
        end
      end
    end
  end

  if nargin >= 3 && ~isempty(id)
    local_warn_cchq0_condition(wf_on_coarse, id);
  end

end

function local_warn_cchq0_condition(wf_on_coarse, id)
% Mirror adaptive Schur-init CCH construction on q=0 slice and warn if ill-conditioned.
  isdf_data = isdftest.get(id);
  if isempty(isdf_data.nrange1) || isempty(isdf_data.nrange2)
    return;
  end

  nrange1 = double(isdf_data.nrange1(:).');
  nrange2 = double(isdf_data.nrange2(:).');
  if isempty(nrange1) || isempty(nrange2)
    return;
  end

  if size(wf_on_coarse, 3) < 1 || size(wf_on_coarse, 4) < 1
    return;
  end

  Psixga = wf_on_coarse(:, nrange1, 1, 1);
  Phixga = wf_on_coarse(:, nrange2, 1, 1);
  CCHq0 = isdftest.prod(Psixga, Psixga, Phixga, Phixga);
  cnd = condest(CCHq0);
  if ~isfinite(cnd) || cnd > 1e6
    warning('isdf:coeff_coarse_wf_extract:CCHq0IllConditioned', ...
      'CCHq(q=0) is ill-conditioned before adaptive update for id=%d: condest=%.6e (>1e6).', ...
      int32(id), cnd);
  end
end
