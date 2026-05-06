% License-Identifier: GPL
%
% Copyright (C) 2026 ZZ
%
% Last modified: 2026/04/19

function gen_bundle(id)

  if nargin < 1 || isempty(id)
    error('isdf.rsymm:adaptive_isdf_closed_r_bundle:Id', 'ISDF id is required.');
  end

  isdf_data = isdf.get(id);
  fft_data = FFT.get();
  wf_data = wave_functions.get();
  k_data = lattice.manager('k', 'get');
  symm_data = symmetry.get();

  nr = double(fft_data.nr);
  nsym = symm_data.nsym;
  fftgrid = single(fft_data.fftgrid);
  if nsym < 1
    error('isdf.rsymm:adaptive_isdf_closed_r_bundle:R_rot', ...
    'FFT.R_rot has no symmetry columns.');
  end

  Nisdf = double(isdf_data.nisdf);
  if Nisdf < 1
    error('isdf.rsymm:adaptive_isdf_closed_r_bundle:Empty',...
         'Empty ISDF sampling (nisdf<1).');
  end

  if isdf_data.interp_scheme ~= "adaptive"
    N_old = 0;
    N_new = Nisdf;
  else
    N_old = isdf_data.N_coarse;
    N_new = isdf_data.N_extra;
  end
  Nrange = N_old+1:N_new+N_old;

  g = double(fftgrid(:)).';
  R_rot = zeros(N_new, nsym);

  for isym = 1:nsym
    rot_mtrx = symm_data.rot_mtrx_RLU_R(:, :, isym);
    SR = round( isdf_data.R_sampling_RLU(N_old+1:N_new+N_old, :) * rot_mtrx );
    ind_SR = mod(SR+fftgrid, fftgrid);
    iv = double(ind_SR);
    R_rot(:, isym) = 1 + iv(:,1) + iv(:,2)*g(1) + iv(:,3)*g(1)*g(2);
  end
  % Then select unique point list in R_rot
  ind_bundle_points = unique(R_rot(:));
  N_bundle = length(ind_bundle_points);
  R_grid_bundle = zeros(N_bundle+N_old, 3);
  R_grid_bundle(N_old+1:N_old+N_bundle, :) = fft_data.Rgrid_RLU(ind_bundle_points, :);
  R_grid_bundle(1:N_old, :) = isdf_data.R_sampling_RLU(1:N_old, :);
  

  % finegrid2bundle is a map, from indices of R_grid_RLU to indices of R_bundle_RLU
  % if is not an element of R_bundle_RLU, be -1
  finegrid2bundle = zeros(nr, 1) - 1;
  bundle2finegrid = zeros(N_bundle, 1);

  for irb = 1:N_bundle
    r_bundle = R_grid_bundle(irb+N_old, :);
    irbinfft_i = 1 + r_bundle(1) + fftgrid(1) * r_bundle(2) ...
          + fftgrid(1) * fftgrid(2) * r_bundle(3);
    finegrid2bundle(irbinfft_i) = irb+N_old;
    bundle2finegrid(irb) = irbinfft_i;
  end % for irb

  % sampling2bundle is a map, from ISDF sampling rows to indices of R_grid_bundle
  sampling2bundle = zeros(N_old+N_new, 1);
  sampling2bundle(1:N_old) = 1:N_old;
  sampling2bundle(N_old+1:N_old+N_new) = finegrid2bundle(R_rot(:, 1));

  % R_rot_in_bundle:
  % R_grid_bundle(R_rot_in_bundle(indR, isym), :) = R_grid_bundle(indR, :) * S(isym)
  % R_grid_RLU(R_bundle2findgrid(indR, isym), :) = R_grid_bundle(indR, :) * S(isym)
  R_rot_in_bundle = zeros(N_bundle+N_old, nsym);
  R_bundle2findgrid = zeros(N_bundle, nsym);
  R_rot_in_bundle(1:N_old, :)= isdf_data.R_rot_coarse;
  for isym = 1:nsym
    rot_mtrx = symm_data.rot_mtrx_RLU_R(:, :, isym);
    M2_r_RLU = R_grid_bundle(N_old+1:N_old+N_bundle, :) * rot_mtrx;
    M2_r_RLU = round(M2_r_RLU);
    iv_mod = double(mod(M2_r_RLU + fftgrid, fftgrid));
    indR = 1 + iv_mod(:,1) + iv_mod(:,2)*fftgrid(1) ...
             + iv_mod(:,3)*fftgrid(1)*fftgrid(2);
    R_bundle2findgrid(:, isym) = indR;
    R_rot_in_bundle(N_old+1:N_old+N_bundle, isym) = finegrid2bundle(indR);
  end


  % Get wavefunction on the bundle
  WF_bundle = zeros(N_bundle + N_old, wf_data.nb, wf_data.nk, wf_data.nspin);
  % WF_bundle(1:N_old, :, :, :) = isdf_data.coeff_seper(:, :, :, :);
  WF_bundle(1:N_old, :, :, :) = isdf_data.tmp(:, :, :, :);
  for ispin = 1:wf_data.nspin
    for ik = 1:wf_data.nk
      for ib = 1:wf_data.nb
        wf_t = wf_data.c(:, ib, ik, ispin);
        WF_bundle(N_old+1:N_old+N_bundle, ib, ik, ispin) = wf_t(bundle2finegrid);
      end
    end
  end

  isdf_data.bundle_struct.N_bundle = N_bundle;
  isdf_data.bundle_struct.R_grid_bundle = R_grid_bundle;
  isdf_data.bundle_struct.R_rot_in_bundle = R_rot_in_bundle;
  isdf_data.bundle_struct.WF_bundle = WF_bundle;
  isdf_data.bundle_struct.sampling2bundle = sampling2bundle;
  isdf.save2mod(isdf_data, id);

  N_coarse = double(isdf_data.N_coarse);
  N_extra = double(isdf_data.N_extra);

  % Do a test, test if the wavefunction on the bundle is correct
  is_t_rev = symm_data.is_t_rev;
  for ispin = 1:wf_data.nspin
    for ikbz = 1:k_data.nbz
      for ib = 3:5
        ikibz = k_data.bz2ibz(ikbz);
        ikrot = double(k_data.bz2rot(ikbz));
        wf_ibz_in_bundle_a = WF_bundle(:, ib, ikibz, ispin);
        % wf_ibz = wf_data.c(:, ib, ikibz, ispin);
        if ikrot > nsym / (is_t_rev + 1)
          ikrot_wf = ikrot -  nsym / (is_t_rev + 1);
          ikrot_wf = symm_data.inv_rot_index(ikrot_wf);
          conjflag = true;
        else
          ikrot_wf = ikrot;
          ikrot_wf = symm_data.inv_rot_index(ikrot_wf);
          conjflag = false;
        end
        indices = R_rot_in_bundle(:, ikrot_wf);
        wf_bz_in_bundle = wf_ibz_in_bundle_a(indices);
        if conjflag
          wf_bz_in_bundle = conj(wf_bz_in_bundle);
        end
        wf_tmp = wf_bz_in_bundle(N_old+1:N_old+N_bundle);
        % direct way, wf_apply_symm + extraction
        isc = [ib, ikibz, ikrot, ispin];
        wf_dir = wave_functions.WF_apply_symm(isc);
        wf_bz_in_bundle_dir = wf_dir(bundle2finegrid);
        if norm(wf_tmp - wf_bz_in_bundle_dir) > 8e-5
          fprintf('The wavefunction on the bundle is not correct\n');
          fprintf('isc = [%d, %d, %d, %d]\n', ib, ikibz, ikrot, ispin);
          fprintf('norm(wf_bz_in_bundle - wf_bz_in_bundle_dir) = %f\n', norm(wf_bz_in_bundle - wf_bz_in_bundle_dir));
          error('The wavefunction on the bundle is not correct');
        end
        % Apply S_q to wf_bundle
        for iq = 1:k_data.nbz
          iqrot = double(k_data.bz2rot(iq));
          % ind = isdf_data.R_rot_extra(N_old+1:N_new+N_old, iqrot);
          % wf_Sq_dir = wf_dir(ind);
          %
          ind = R_rot_in_bundle(:, iqrot);
          wf_Sq_bundle = wf_bz_in_bundle(ind);
          wf_Sq_xalpha = wf_Sq_bundle(sampling2bundle(N_old+1:N_old+N_new));
          % wf_tmp = wf_Sq_xalpha(N_old+1:N_old+N_bundle);
          % Direct way.
          ind = isdf_data.R_rot_extra(N_coarse + 1:N_coarse + N_extra, iqrot);
          wf_Sq_xalpha_dir = wf_dir(ind);
          if norm(wf_Sq_xalpha - wf_Sq_xalpha_dir) > 8e-5
            error('The wavefunction on the bundle is not correct');
          end
        end
      end

    end
  end

end
