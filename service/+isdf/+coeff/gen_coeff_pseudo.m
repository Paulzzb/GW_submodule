% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/23

function gen_coeff_pseudo(id)
% ISDF pseudo index initializer:
%   select one seed by maximal adaptive weight on the full fine grid, then
%   take the whole symmetry orbit of that seed as the pseudo initial set.
%   Orbits with CCH condition number >= 1e+10 are rejected; the next-largest
%   seed outside rejected orbits is tried until a well-conditioned set is found.

  fft_data = FFT.get();
  wf_data = wave_functions.get();
  symm_data = symmetry.get();
  isdf_data = isdf.get(id);

  if isempty(isdf_data.nrange1) || isempty(isdf_data.nrange2)
    error('isdf:gen_coeff_pseudo:nrange', ...
      'Missing cached nrange in ISDF id=%d. Build it first via isdf.set_nrange(id, config.SYSTEM).', int32(id));
  end

  nrange1 = double(isdf_data.nrange1(:).');
  nrange2 = double(isdf_data.nrange2(:).');
  if isempty(nrange1) || isempty(nrange2)
    error('isdf:gen_coeff_pseudo:nrange', 'nrange1/nrange2 must be non-empty.');
  end

  k_data = lattice.manager('k', 'get');
  Nw = double(wf_data.nc);
  Psi = zeros(Nw, numel(nrange1) * double(k_data.nbz), 'double');
  Phi = zeros(Nw, numel(nrange2) * double(k_data.nbz), 'double');
  count1 = 1;
  count2 = 1;
  ispin = 1;
  for ikbz = 1:k_data.nbz
    ikibz = k_data.bz2ibz(ikbz, 1);
    ikrot = k_data.bz2rot(ikbz, 1);
    for ib1 = 1:numel(nrange1)
      param1 = [nrange1(ib1), ikibz, ikrot, ispin];
      Psi(:, count1) = double(wave_functions.WF_apply_symm(param1));
      count1 = count1 + 1;
    end
    for ib2 = 1:numel(nrange2)
      param2 = [nrange2(ib2), ikibz, ikrot, ispin];
      Phi(:, count2) = double(wave_functions.WF_apply_symm(param2));
      count2 = count2 + 1;
    end
  end

  w = zeros(Nw, 1, 'double');
  for ic = 1:Nw
    psi_r = Psi(ic, :);
    phi_r = Phi(ic, :);
    w(ic) = real(isdf.prod(psi_r, psi_r, phi_r, phi_r));
  end
  gc_s = double(fft_data.fftgrid(:).');
  max_cch_cond = 1e+10;
  excluded_mask = false(double(fft_data.nr), 1);
  orbit_idx = int32([]);
  seed_idx = 0;

  while true
    w_pick = w;
    w_pick(excluded_mask(1:Nw)) = -inf;
    [~, seed_idx] = max(w_pick);
    if ~isfinite(w_pick(seed_idx)) || w_pick(seed_idx) < 0
      error('isdf:gen_coeff_pseudo:NoValidOrbit', ...
        'No pseudo orbit with CCH condition number below %.1e.', max_cch_cond);
    end

    orbit_mask = local_build_symmetry_orbit(seed_idx, fft_data, symm_data, gc_s);
    orbit_idx = int32(find(orbit_mask));
    if local_orbit_cch_ok(orbit_idx, Psi, Phi, max_cch_cond)
      break;
    end

    excluded_mask(orbit_mask) = true;
    fprintf(2, ...
      'gen_coeff_pseudo: reject seed %d (orbit size %d, CCH cond >= %.1e).\n', ...
      seed_idx, numel(orbit_idx), max_cch_cond);
  end

  Nmu = int32(numel(orbit_idx));
  if Nmu < 1
    error('isdf:gen_coeff_pseudo:EmptyOrbit', 'Pseudo orbit is empty for seed index %d.', seed_idx);
  end

  R_coarse_RLU = double(fft_data.Rgrid_RLU(orbit_idx, :));
  wf_on_coarse = wf_data.c(orbit_idx, :, :, :);

  lin2loc = zeros(double(fft_data.nr), 1, 'int32');
  lin2loc(double(orbit_idx)) = int32(1:double(Nmu));
  R_rot_coarse = zeros(double(Nmu), symm_data.nsym, 'int32');
  for isym = 1:symm_data.nsym
    M2 = symm_data.rot_mtrx_RLU_R(:, :, isym);
    M2_r_RLU = (R_coarse_RLU * M2);
    if norm(M2_r_RLU - round(M2_r_RLU)) > double(1e-4)
      error('isdf:gen_coeff_pseudo:NonIntegerMap', ...
        'Non-integer mapping under rotation; check rot_mtrx_RLU_R and fftgrid.');
    end
    M2_r_RLU = round(M2_r_RLU);
    iv_mod = int32(mod(M2_r_RLU + gc_s, gc_s));
    dm = double(iv_mod);
    g1 = gc_s(1);
    g12 = gc_s(1) * gc_s(2);
    rot_lin = int32(round(1 + dm(:, 1) + dm(:, 2) * g1 + dm(:, 3) * g12));
    rot_loc = lin2loc(double(rot_lin));
    if any(rot_loc <= 0)
      error('isdf:gen_coeff_pseudo:OrbitNotClosed', ...
        'Pseudo orbit is not closed under symmetry isym=%d.', isym);
    end
    R_rot_coarse(:, isym) = rot_loc;
  end

  isdf_data.nisdf = Nmu;
  isdf_data.N_coarse = Nmu;
  isdf_data.N_extra = int32(0);
  isdf_data.coeff_seper = wf_on_coarse;
  isdf_data.fftgrid_c = int32(fft_data.fftgrid(:).');
  isdf_data.R_sampling_RLU = R_coarse_RLU;
  isdf_data.interp_scheme = "pseudo";
  isdf_data.R_rot_coarse = R_rot_coarse;

  bundle_struct = struct();
  bundle_struct.N_bundle = isdf_data.N_coarse;
  bundle_struct.N_coarse = isdf_data.N_coarse;
  bundle_struct.N_sampling = isdf_data.N_coarse;
  bundle_struct.R_grid_bundle = R_coarse_RLU;
  bundle_struct.R_rot_in_bundle = R_rot_coarse;
  bundle_struct.WF_bundle = wf_on_coarse;
  bundle_struct.sampling2bundle = int32((1:double(isdf_data.N_coarse)).');
  bundle_struct.fine_grid_lin = orbit_idx;
  isdf_data.bundle_struct = bundle_struct;

  isdf.save2mod(isdf_data, id);
end

function orbit_mask = local_build_symmetry_orbit(seed_idx, fft_data, symm_data, gc_s)
  orbit_mask = false(double(fft_data.nr), 1);
  orbit_mask(seed_idx) = true;
  frontier = int32(seed_idx);
  while ~isempty(frontier)
    nxt = int32([]);
    for isym = 1:symm_data.nsym
      M2 = symm_data.rot_mtrx_RLU_R(:, :, isym);
      R_frontier = double(fft_data.Rgrid_RLU(frontier, :));
      M2_r_RLU = (R_frontier * M2);
      if norm(M2_r_RLU - round(M2_r_RLU)) > double(1e-4)
        error('isdf:gen_coeff_pseudo:NonIntegerMap', ...
          'Non-integer mapping under rotation; check rot_mtrx_RLU_R and fftgrid.');
      end
      M2_r_RLU = round(M2_r_RLU);
      iv_mod = int32(mod(M2_r_RLU + gc_s, gc_s));
      dm = double(iv_mod);
      g1 = gc_s(1);
      g12 = gc_s(1) * gc_s(2);
      rot_lin = int32(round(1 + dm(:, 1) + dm(:, 2) * g1 + dm(:, 3) * g12));
      new_lin = rot_lin(~orbit_mask(double(rot_lin)));
      if ~isempty(new_lin)
        orbit_mask(double(new_lin)) = true;
        nxt = [nxt; new_lin(:)]; %#ok<AGROW>
      end
    end
    frontier = unique(nxt, 'stable');
  end
end

function ok = local_orbit_cch_ok(orbit_idx, Psi, Phi, max_cch_cond)
  idx = double(orbit_idx);
  if isempty(idx)
    ok = false;
    return;
  end
  Psi_orbit = Psi(idx, :);
  Phi_orbit = Phi(idx, :);
  CCH = isdf.prod(Psi_orbit, Psi_orbit, Phi_orbit, Phi_orbit);
  CCH = 0.5*CCH + 0.5*CCH';
  if condest(CCH) >= max_cch_cond
    ok = false;
    return;
  end
  try
    chol(CCH, 'lower');
    ok = true;
  catch
    ok = false;
  end
end
