% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/20 ZZ

function idnew = adaptiveisdf(id, cfg_isdf)
  % Wall-clock time for the main adaptive phase through the orbit summary (before Verification).
  t_phase1 = tic;

  if nargin < 2
    cfg_isdf = struct();
  end

  % 1. Starting grid: from coarse fft grid
  system_data = system.get();
  nspin = system_data.nspin;
  k_data = lattice.manager('k', 'get');
  if nspin > 1
    error('adaptiveisdf:nspin', 'nspin > 1 is not supported.');
  end
  wf_data = wave_functions.get();

  [ck_loaded, idnew] = isdf.adaptive_double.adaptive_checkpoint_try_load(id);
  if ck_loaded
    return;
  end

  import_root = local_get_import_isdf_root(cfg_isdf);
  if ~isempty(import_root)
    [ck_loaded, idnew] = isdf.adaptive_double.adaptive_checkpoint_import_from_isdf(id, import_root);
    if ck_loaded
      return;
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Remove redundant points
  % isdf.adaptive_double.isdf_exclude_point(id);
  isdf_data = isdf.get(id);
  Nisdf = (isdf_data.nisdf);
  params = adaptiveisdf_resolve_params(isdf_data.desc, cfg_isdf);
  threshold = params.threshold;
  num_add = params.num_add;
  ratio = params.candidate_ratio;
  adaptive_backend = 'adaptive_double';
  adaptive_arithmetic = 'double';
  adaptiveisdf_print_run_config(id, isdf_data, params, threshold, num_add, ratio, ...
    adaptive_backend, adaptive_arithmetic);
  isdf.adaptive_double.adaptive_weight('set_batch_size', params.weight_batch_size);
  
  if isempty(isdf_data.nrange1) || isempty(isdf_data.nrange2)
    error('adaptiveisdf:nrange', ...
      'Missing cached nrange in ISDF id=%d. Build it first via isdf.set_nrange(id, config.SYSTEM).', ...
      int32(id));
  end
  nrange1 = double(isdf_data.nrange1);
  nrange2 = double(isdf_data.nrange2);
  nmu_cap = params.isdf_ratio * sqrt(k_data.nbz) * params.max_add_frac * sqrt(double(numel(nrange1)) * double(numel(nrange2)));
  Naddmax = int32(max(0, ceil(nmu_cap - double(Nisdf))));
  R_sampling_indices = isdf_data.bundle_struct.sampling2bundle;
  if isfield(isdf_data.bundle_struct, 'fine_grid_lin') ...
      && ~isempty(isdf_data.bundle_struct.fine_grid_lin)
    R_sampling_indices = isdf_data.bundle_struct.fine_grid_lin;
  end

  % Build full-grid Psi/Phi (all BZ k via WF_apply_symm), then sample centroids.
  isdf.adaptive_double.adaptive_weight('init', id);
  [Psixga, Phixga] = isdf.adaptive_double.adaptive_weight('get_wf_xga', R_sampling_indices);

  isdf.adaptive_double.isdf_schur_update('init', Nisdf, Psixga, Phixga, ...
    params.max_cond, params.use_cond_guard);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 
  % 3. Calculate the weight based on new set of isdf points,
  %    and calculate the loss function 
  w = isdf.adaptive_double.adaptive_weight('get');
  loss0 = sum(w);
  fprintf('Initial loss: %f\n', loss0);
  isdf.adaptive_double.adaptive_weight('init_update');
  w = isdf.adaptive_double.adaptive_weight('get');
  loss = sum(w);
  fprintf('loss after initial update: %f\n', loss);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Now, do the update
  % 

  num_add_i = double(num_add);
  isdf_new_indices = zeros(double(Naddmax), 1);
  global_index_remain = 1:length(w);
  N_candidate = int32(max(1, ceil(double(num_add) * double(ratio))));
  candidate_indices = zeros(double(N_candidate), 1);
  Nisdf_new = 0;
  Nisdf_new_max = Naddmax;
  selected_indices = zeros(num_add_i, 1);
  w_candidate = zeros(double(N_candidate), 1);
  fft_data = FFT.get();
  symm_data = symmetry.get();
  n_iter = 0;
  loss_history = zeros(double(Naddmax) + 1, 1);
  rel_loss_history = zeros(double(Naddmax) + 1, 1);
  loss_history(1) = loss;
  if loss0 > 0
    rel_loss_history(1) = loss / loss0;
  else
    rel_loss_history(1) = NaN;
  end


  breakflag = false;
  schur_converged = false;
  while Nisdf_new < Nisdf_new_max && loss0 > 0 && (loss / loss0 > threshold) && ~isempty(global_index_remain)
    n_iter = n_iter + 1;
    % 1. Select candidate set
    N_candidate_eff = min(double(N_candidate), numel(global_index_remain));
    [~, idx_sorted] = sort(w(global_index_remain), 'descend');
    % candidate set
    candidate_indices(1:N_candidate_eff) = global_index_remain(idx_sorted(1:N_candidate_eff));
    % 2. Select num_candidate from candidate set, also update the schur matrices
    % selected_set = isdf_select_from_candidate();
    Nremain = int32(N_candidate_eff);
    selected_indices(1:num_add_i) = 0;
    candidate_indices_saved = candidate_indices(1:N_candidate_eff);
    n_selected_iter = 0;
    n_schur_fail_this_iter = 0;
    for iadd = 1:num_add_i
      if Nremain <= 0
        break;
      end
      schur_ok = false;
      while Nremain > 0 && ~schur_ok
        % 2.1 select the point with the highest weight among remaining candidates
        w_candidate(1:Nremain) = isdf.adaptive_double.adaptive_weight('get', candidate_indices(1:Nremain));
        if all(w_candidate(1:Nremain) < 0)
          break;
        end
        [~, idx_local] = sort(w_candidate(1:Nremain), 'ascend');
        candidate_indices(1:Nremain) = candidate_indices(idx_local);
        selected_indices(iadd) = candidate_indices(Nremain);
        % 2.2 update the Gram matrix (all BZ k via WF_apply_symm in adaptive_weight)
        [Psi_on_new_grid, Phi_on_new_grid] = ...
          isdf.adaptive_double.adaptive_weight('get_wf_xga', selected_indices(iadd));
        schur_ok = isdf.adaptive_double.isdf_schur_update('update', 1, Psi_on_new_grid, Phi_on_new_grid, selected_indices(iadd));
        if ~schur_ok
          if loss0 > 0
            rel_now = loss / loss0;
          else
            rel_now = NaN;
          end
          fprintf('adaptiveisdf: skip index %d (Schur guard failed), relative loss = %.6e\n', ...
            selected_indices(iadd), rel_now);
          n_schur_fail_this_iter = n_schur_fail_this_iter + 1;
          Nremain = Nremain - 1;
          selected_indices(iadd) = 0;
          continue;
        end
        % 2.3 update the weight
        isdf.adaptive_double.adaptive_weight('update', candidate_indices(1:Nremain), 1);
        Nremain = Nremain - 1;
        Nisdf_new = Nisdf_new + 1;
        isdf_new_indices(Nisdf_new) = selected_indices(iadd);
        n_selected_iter = iadd;
        if Nisdf_new == Nisdf_new_max
          breakflag = true;
          break;
        end
      end
      if ~schur_ok
        break;
      end
      if breakflag
        break;
      end
    end
    if n_selected_iter == 0 && n_schur_fail_this_iter > 0
      fprintf('adaptiveisdf: all candidate Schur updates failed this iteration; treating as converged.\n');
      schur_converged = true;
      break;
    end
    % 3. Update the global weight
    %    In case early deflation, only iadd new indices are chosen.
    selected_now = selected_indices(1:n_selected_iter);
    selected_now = selected_now(selected_now > 0);
    global_index_remain = setdiff(global_index_remain, selected_now);
    tocal = setdiff(global_index_remain, candidate_indices_saved);
    isdf.adaptive_double.adaptive_weight('update', tocal, n_selected_iter);
    % 4. Update the loss
    w(global_index_remain) = isdf.adaptive_double.adaptive_weight('get', global_index_remain);
    indt = w(global_index_remain) < 0;
    zero_weight_indices = global_index_remain(indt);
    % zero_weight_indices = find(w < 0); % Because of truncated error
    global_index_remain = setdiff(global_index_remain, zero_weight_indices);
    loss = sum(w(global_index_remain));
    loss_history(n_iter + 1) = loss;
    rel_loss_history(n_iter + 1) = loss / loss0;
    fprintf('Loss: %f\n', loss);
    fprintf('Relative loss: %f\n', loss / loss0);
    if breakflag || schur_converged
      break;
    end
  end


  % Generate a concise report for adaptive updating.
  if loss0 > 0
    rel_loss = loss / loss0;
  else
    rel_loss = NaN;
  end

  if ~(loss0 > 0)
    stop_reason = 'invalid_initial_loss';
  elseif Nisdf_new_max <= 0
    stop_reason = 'naddmax_zero';
  elseif rel_loss <= threshold
    stop_reason = 'reach_threshold';
  elseif isempty(global_index_remain)
    stop_reason = 'no_remaining_points';
  elseif schur_converged
    stop_reason = 'schur_converged';
  else
    stop_reason = 'loop_guard';
  end

  n_report = min(Nisdf_new, 10);
  report = struct();
  report.initial_nisdf = Nisdf;
  report.added_nisdf = Nisdf_new;
  report.final_nisdf = Nisdf + Nisdf_new;
  report.iterations = n_iter;
  report.loss_initial = loss0;
  report.loss_final = loss;
  report.relative_loss = rel_loss;
  report.threshold = threshold;
  report.num_add = num_add;
  report.candidate_ratio = ratio;
  report.max_add_frac = params.max_add_frac;
  report.nmu_cap = nmu_cap;
  report.naddmax = Naddmax;
  report.isdf_ratio = params.isdf_ratio;
  report.max_cond = params.max_cond;
  report.use_cond_guard = logical(params.use_cond_guard);
  report.param_source = params.source;
  report.stop_reason = stop_reason;
  report.adaptive_backend = adaptive_backend;
  report.adaptive_arithmetic = adaptive_arithmetic;
  report.isdf_desc = char(string(isdf_data.desc));
  report.new_indices = isdf_new_indices(1:Nisdf_new);
  report.loss_history = loss_history(1:n_iter+1);
  report.relative_loss_history = rel_loss_history(1:n_iter+1);

  fprintf('\n=== Adaptive ISDF Update Report ===\n');
  fprintf('Backend            : %s\n', adaptive_backend);
  fprintf('Arithmetic         : %s\n', adaptive_arithmetic);
  fprintf('ISDF desc          : %s\n', char(string(isdf_data.desc)));
  fprintf('Coarse ISDF id     : %d\n', int32(id));
  fprintf('Initial Nisdf      : %d\n', Nisdf);
  fprintf('Added points       : %d\n', Nisdf_new);
  fprintf('Final Nisdf        : %d\n', report.final_nisdf);
  fprintf('Iterations         : %d\n', n_iter);
  fprintf('Initial loss       : %.8e\n', loss0);
  fprintf('Final loss         : %.8e\n', loss);
  fprintf('Relative loss      : %.8e\n', rel_loss);
  fprintf('Threshold          : %.8e\n', threshold);
  fprintf('num_add            : %d\n', num_add);
  fprintf('candidate ratio    : %.4f\n', ratio);
  fprintf('max add frac       : %.4f\n', params.max_add_frac);
  fprintf('ISDF ratio         : %.4f\n', params.isdf_ratio);
  fprintf('max cond           : %.4e\n', params.max_cond);
  fprintf('use cond guard     : %d\n', logical(params.use_cond_guard));
  fprintf('Nisdf cap          : %.8e\n', nmu_cap);
  fprintf('Naddmax            : %d\n', Naddmax);
  fprintf('param source       : %s\n', params.source);
  fprintf('Stop reason        : %s\n', stop_reason);
  if n_report > 0
    fprintf('New indices (first %d): ', n_report);
    fprintf('%d ', isdf_new_indices(1:n_report));
    fprintf('\n');
  else
    fprintf('New indices        : <none>\n');
  end
  fprintf('Stepwise loss/loss0 after each +Nadd update:\n');
  if n_iter == 0
    fprintf('  Step 0: loss/loss0 = %.8e\n', rel_loss_history(1));
  else
    for istep = 1:n_iter
      fprintf('  Step %d: loss = %.8e, loss/loss0 = %.8e\n', ...
              istep, loss_history(istep+1), rel_loss_history(istep+1));
    end
  end
  fprintf('===================================\n\n');

  % --- Symmetry orbits on the fine FFT grid for newly added centroids ---
  % (1) How many distinct orbits among the new seeds.
  % (2) Total fine-grid points in the union of those orbits.
  % (3) How many of those orbit points lie in the current centroid set (initial + new).
  Nadded = double(Nisdf_new);
  if Nadded > 0
    seeds = unique(isdf_new_indices(1:Nadded));
    fft_data_orb = FFT.get();
    R_rot = double(fft_data_orb.R_rot);
    nr_orb = double(fft_data_orb.nr);
    U_orb = adaptiveisdf_orbit_union_bfs(seeds, R_rot, nr_orb);
    n_pts_orbit_union = nnz(U_orb);

    canon = zeros(numel(seeds), 1);
    for ks = 1:numel(seeds)
      Om = adaptiveisdf_orbit_mask_bfs(seeds(ks), R_rot, nr_orb);
      canon(ks) = min(find(Om));
    end
    n_orbits_among_new = numel(unique(canon));

    isdf_data_cur = isdf.get(id);
    Nold_mu = double(isdf_data_cur.nisdf);
    fine_old = adaptiveisdf_r_sampling_rows_to_lin(isdf_data_cur, fft_data_orb, Nold_mu);
    current_mask = false(nr_orb, 1);
    current_mask(fine_old) = true;
    current_mask(seeds) = true;
    n_in_current = nnz(U_orb & current_mask);

    current_old_only = false(nr_orb, 1);
    current_old_only(fine_old) = true;
    n_in_initial_only = nnz(U_orb & current_old_only);

    fprintf('--- Orbit report (fine FFT grid, spatial symmetries via FFT.R_rot) ---\n');
    fprintf('  New centroid count (with possible repeats): %d\n', Nadded);
    fprintf('  Distinct new seeds (fine linear indices): %d\n', numel(seeds));
    fprintf('  (1) Distinct symmetry orbits among new seeds: %d\n', n_orbits_among_new);
    fprintf('  (2) Total fine-grid points in union of those orbits: %d\n', n_pts_orbit_union);
    fprintf('  (3a) Of those orbit points, in initial centroid set only (first Nisdf=%d): %d\n', ...
      Nisdf, n_in_initial_only);
    fprintf('  (3b) Of those orbit points, in full current set (initial + new seeds): %d\n', ...
      n_in_current);
    fprintf('--- end orbit report ---\n\n');
    report.orbit_has_data = true;
    report.orbit_Nadded = Nadded;
    report.orbit_seeds_count = numel(seeds);
    report.orbit_n_orbits = n_orbits_among_new;
    report.orbit_n_pts_union = n_pts_orbit_union;
    report.orbit_n_in_initial_only = n_in_initial_only;
    report.orbit_n_in_current = n_in_current;
  else
    fprintf('--- Orbit report: no new centroids (Nadded=0), skip. ---\n\n');
    report.orbit_has_data = false;
  end

  elapsed_phase1 = toc(t_phase1);
  report.coarse_isdf_id = id;
  report.elapsed_phase1_seconds = elapsed_phase1;
  fpath_r = adaptiveisdf_write_phase1_report(report);
  fprintf(1, 'adaptiveisdf: wrote phase-1 report to %s (elapsed %.6f s, through orbit summary)\n', ...
    fpath_r, elapsed_phase1);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Verification
  [Nisdf_new, ~, ~, ~, Psi_on_grid, Phi_on_grid] = ...
         isdf.adaptive_double.isdf_schur_update('get');
  
  Nextra = Nisdf_new - Nisdf;
  % Check if the new-added points correcly indiced.
  % for i = 1:Nextra
  %   ind_new = isdf_new_indices(i);
  %   [Psi_on_new_grid, Phi_on_new_grid] = isdf.adaptive_double.adaptive_weight('get_wf_xga', ind_new);
  %   Psi_imported = Psi_on_grid(Nisdf+i, :);
  %   Phi_imported = Phi_on_grid(Nisdf+i, :);
  %   if norm(Psi_on_new_grid - Psi_imported) > 1e-6
  %     error('adaptiveisdf:new_indices', 'The new-added points are not correctly indexed.');
  %   end
  %   if norm(Phi_on_new_grid - Phi_imported) > 1e-6
  %     error('adaptiveisdf:new_indices', 'The new-added points are not correctly indexed.');
  %   end
  % end
  %
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CCH = isdf.prod(Psi_on_grid(1:Nisdf_new, :), Psi_on_grid(1:Nisdf_new, :), ...
  %                 Phi_on_grid(1:Nisdf_new, :), Phi_on_grid(1:Nisdf_new, :));
  % L_CCH_dir = chol(CCH(1:Nisdf_new, 1:Nisdf_new), "lower");
  % output = norm(L_CCH_dir - L_CCH(1:Nisdf_new, 1:Nisdf_new), 'fro') / norm(L_CCH(1:Nisdf_new, 1:Nisdf_new), 'fro');
  % fprintf('Difference between direct Cholesky and adaptive Cholesky: %f\n', output);

  % output = norm(...
  %          invL_CCH(1:Nisdf_new, 1:Nisdf_new) * L_CCH(1:Nisdf_new, 1:Nisdf_new)...
  %           - eye(Nisdf_new), 'fro') / norm(eye(Nisdf_new), 'fro');
  % fprintf('Difference between invL_CCH * L_CCH and eye: %f\n', output);

  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  idnew = isdf.isdf_add(isdf_data.desc);
  % Build a new ISDF object at idnew, keep original id unchanged.
  isdf_data_new = isdf.get(idnew);
  isdf_data_new.nrange1 = isdf_data.nrange1;
  isdf_data_new.nrange2 = isdf_data.nrange2;
  isdf_data_new.Nnrange1 = isdf_data.Nnrange1;
  isdf_data_new.Nnrange2 = isdf_data.Nnrange2;
  isdf_data_new = isdf_data;
  isdf_data_new.id = idnew;
  % isdf_data_new.desc = sprintf('adaptiveisdf from id=%d', int32(id));
  isdf_data_new.desc = isdf_data.desc;
  % isdf_data_new.desc = sprintf('');
  isdf_data_new.nisdf = int32(Nisdf_new);
  isdf_data_new.interp_scheme = 'adaptive';
  isdf_data_new.N_coarse = int32(Nisdf);
  isdf_data_new.N_extra = int32(Nextra);

  % Append newly selected R points.
  if Nextra > 0
    new_lin = int32(isdf_new_indices(1:Nextra));
    R_new = double(fft_data.Rgrid_RLU(new_lin, :));
    isdf_data_new.R_sampling_RLU = double([isdf_data.R_sampling_RLU; R_new]);
    isdf_data_new.R_rot_extra = zeros(Nisdf_new, symm_data.nsym, 'int32');
    for is = 1:int32(symm_data.nsym)
      % if is > symm_data.nsym / (1 + symm_data.is_t_rev)
      %   is_t = is - symm_data.nsym / (1 + symm_data.is_t_rev);
      % else
      %   is_t = is;
      % end
      is_t = is;
      isdf_data_new.R_rot_extra(1:Nisdf, is) = isdf_data.R_rot_coarse(1:Nisdf, is_t);
      % Select appropriate rotation matrix
      mtrx_RLU_R = symm_data.rot_mtrx_RLU_R(:, :, is);
      M2 = double( mtrx_RLU_R );
      M2_r_RLU = (fft_data.Rgrid_RLU * M2);  % nr x 3
      if norm(M2_r_RLU - round(M2_r_RLU)) > 1e-3
        error('Non-integer mapping found in rotation. Check the rotation matrices and FFT grid.');
      end
      M2_r_RLU = int32(round(M2_r_RLU));
      % Update R_rot_extra for extra grid points
      iv_mod = int32(mod(M2_r_RLU + fft_data.fftgrid, fft_data.fftgrid));  % nr x 3
      i4 = 1 + iv_mod(:,1) + iv_mod(:,2)*fft_data.fftgrid(1) + iv_mod(:,3)*fft_data.fftgrid(1)*fft_data.fftgrid(2);  % nr x 1
      isdf_data_new.R_rot_extra(Nisdf+1:Nisdf_new, is) = int32(i4(isdf_new_indices(1:Nextra)));
      % % Update R_rot_extra for coarse grid points
      % R_coarse_scal = double(isdf_data.fftgrid_c) ./ double(fft_data.fftgrid);
      % R_coarse_scal_RLU = isdf_data.R_sampling_RLU(1:Nisdf, :) .* R_coarse_scal;
      % if norm(R_coarse_scal_RLU - round(R_coarse_scal_RLU)) > 1e-4
      %   error('Non-integer mapping found in rotation. Check the rotation matrices and FFT grid.');
      % end
      % R_coarse_scal_RLU = round(R_coarse_scal_RLU);
      % M2_r_RLU_c = round(R_coarse_scal_RLU * M2);
      % if norm(M2_r_RLU_c - round(M2_r_RLU_c)) > 1e-4
      %   error('Non-integer mapping found in rotation. Check the rotation matrices and FFT grid.');
      % end
      % M2_r_RLU_c = round(M2_r_RLU_c);
      % iv_mod_c = int32(mod(M2_r_RLU_c + fftgrid, fftgrid));  % nr x 3
      % i4_c = 1 + iv_mod_c(:,1) + iv_mod_c(:,2)*isdf_data.fftgrid_c(1) ...
      % + iv_mod_c(:,3)*isdf_data.fftgrid_c(1)*isdf_data.fftgrid_c(2);
      % isdf_data_new.R_rot_extra(1:Nisdf, is) = ...
      % int32(i4_c);
    end 


  end

  % Rebuild coeff_seper with enlarged first dimension.
  coeff_old = isdf_data.coeff_seper;
  nb = size(coeff_old, 2);
  % nk = size(coeff_old, 3);
  nbz = k_data.nbz; 
  nspin_coeff = size(coeff_old, 4);
  coeff_new = zeros(double(Nisdf_new), nb, nbz, nspin_coeff);
  % coeff_new(1:double(Nisdf), :, :, :) = coeff_old(1:double(Nisdf), :, :, :);
  for ikbz = 1:nbz
    ikibz = k_data.bz2ibz(ikbz, 1);
    ikrot = k_data.bz2rot(ikbz, 1);
    wf_c_1= coeff_old(:, :, ikibz, 1);
    if ikrot ~= 1
      wf_c_1 = isdf.coeff.isdf_apply_symm_on_coarse(id, wf_c_1, ikrot);
    end
    coeff_new(1:Nisdf, :, ikbz, 1) = wf_c_1;
  end

  ispin = 1;
  for ikbz = 1:nbz
    ikibz = k_data.bz2ibz(ikbz, 1);
    ikrot = k_data.bz2rot(ikbz, 1);
    for ib = 1:nb
      isc = [ib, ikibz, ikrot, ispin];
      wf_1 = wave_functions.WF_apply_symm(isc);
      coeff_new(Nisdf+1:Nisdf_new, ib, ikbz, 1) = wf_1(isdf_new_indices(1:Nextra));
    end
  end
  % isdf_data_new.coeff_seper = isdf_data.coeff_seper;
  isdf_data_new.coeff_seper = coeff_new;
  isdf_data_new.tmp = isdf_data.coeff_seper;
  %
  % R_new = double(fft_data.Rgrid_RLU(isdf_new_indices(1:Nextra), :));
  % isdf_data_new.R_sampling_RLU = [isdf_data.R_sampling_RLU; R_new];
  isdf_data_new.R_rot_coarse = isdf_data.R_rot_coarse;

  % if Nextra > 0
  %   for isym = 1:symm_data.nsym
  %     isdf_data_new.R_rot_extra(:, isym) = ...
  %     fft_data.R_rot(isdf_new_indices(1:Nextra), isym);
  %   end
  % end
  
  isdf_data_new.bundle_struct = isdf_data.bundle_struct;
  isdf_data_new.assigned = true;
  isdf.save2mod(isdf_data_new, idnew);
  % idnew = isdf.isdf_add(isdf_data.desc);
  isdf.rsymm.bundle_refresh(id, Nextra, isdf_new_indices, idnew);
  % isdf.rsymm.gen_bundle(idnew);

  % Mark heavy q-dependent tensors dirty to avoid using stale dimensions.
  % isdf_data_new.tildeVq = zeros(0, 0, 0, 0);
  % isdf_data_new.helperqG = zeros(0, 0, 0, 0);
  

  % validation using changing
  % isdf.changing(idnew);
  % isdf.gen_tildeVq(idnew);
  % isdf.validation.isdf_validation(idnew);

  isdf.adaptive_double.adaptive_checkpoint_save(id, idnew);
  isdf.numerical_cond_report('adaptive', id, idnew, isdf_data.desc, ...
    loss_history(1), loss, report.final_nisdf);
  fprintf('Saved adaptive ISDF to new id = %d\n', idnew);

end

function fpath = adaptiveisdf_write_phase1_report(r)
% Write adaptive loop + orbit summary to adaptiveisdf_id<coarse_isdf_id>.txt in pwd.

  cid = double(r.coarse_isdf_id);
  fname = sprintf('adaptiveisdf_id%d.txt', cid);
  fpath = fullfile(pwd, fname);
  fid = fopen(fpath, 'w');
  if fid < 0
    error('adaptiveisdf:reportOpen', 'Cannot open for write: %s', fpath);
  end
  oc = onCleanup(@() fclose(fid)); %#ok<NASGU>

  fprintf(fid, '=== adaptiveisdf phase-1 report (coarse ISDF id = %d) ===\n', cid);
  fprintf(fid, 'Generated: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
  fprintf(fid, ['Elapsed wall time (tic/toc): from function start through orbit summary ', ...
    '(excludes Verification and later): %.9f s\n\n'], double(r.elapsed_phase1_seconds));

  fprintf(fid, '\n=== Adaptive ISDF Update Report ===\n');
  if isfield(r, 'adaptive_backend')
    fprintf(fid, 'Backend            : %s\n', char(string(r.adaptive_backend)));
  end
  if isfield(r, 'adaptive_arithmetic')
    fprintf(fid, 'Arithmetic         : %s\n', char(string(r.adaptive_arithmetic)));
  end
  if isfield(r, 'isdf_desc')
    fprintf(fid, 'ISDF desc          : %s\n', char(string(r.isdf_desc)));
  end
  fprintf(fid, 'Coarse ISDF id     : %d\n', cid);
  fprintf(fid, 'Initial Nisdf      : %d\n', int32(r.initial_nisdf));
  fprintf(fid, 'Added points       : %d\n', int32(r.added_nisdf));
  fprintf(fid, 'Final Nisdf        : %d\n', int32(r.final_nisdf));
  fprintf(fid, 'Iterations         : %d\n', int32(r.iterations));
  fprintf(fid, 'Initial loss       : %.8e\n', r.loss_initial);
  fprintf(fid, 'Final loss         : %.8e\n', r.loss_final);
  fprintf(fid, 'Relative loss      : %.8e\n', r.relative_loss);
  fprintf(fid, 'Threshold          : %.8e\n', r.threshold);
  fprintf(fid, 'num_add            : %d\n', int32(r.num_add));
  fprintf(fid, 'candidate ratio    : %.8e\n', r.candidate_ratio);
  fprintf(fid, 'max add frac       : %.8e\n', r.max_add_frac);
  fprintf(fid, 'ISDF ratio         : %.8e\n', r.isdf_ratio);
  fprintf(fid, 'max cond           : %.8e\n', r.max_cond);
  if isfield(r, 'use_cond_guard')
    fprintf(fid, 'use cond guard     : %d\n', logical(r.use_cond_guard));
  end
  fprintf(fid, 'Nisdf cap          : %.8e\n', r.nmu_cap);
  fprintf(fid, 'Naddmax            : %d\n', int32(r.naddmax));
  fprintf(fid, 'param source       : %s\n', char(string(r.param_source)));
  fprintf(fid, 'Stop reason        : %s\n', char(string(r.stop_reason)));
  n_rep = min(int32(r.added_nisdf), int32(10));
  if n_rep > 0
    fprintf(fid, 'New indices (first %d): ', double(n_rep));
    ni = r.new_indices(1:double(n_rep));
    fprintf(fid, '%d ', ni);
    fprintf(fid, '\n');
  else
    fprintf(fid, 'New indices        : <none>\n');
  end
  fprintf(fid, 'Stepwise loss/loss0 after each +Nadd update:\n');
  n_it = int32(r.iterations);
  lh = r.loss_history;
  rlh = r.relative_loss_history;
  if n_it == 0
    fprintf(fid, '  Step 0: loss/loss0 = %.8e\n', rlh(1));
  else
    for istep = 1:double(n_it)
      fprintf(fid, '  Step %d: loss = %.8e, loss/loss0 = %.8e\n', istep, lh(istep + 1), rlh(istep + 1));
    end
  end
  fprintf(fid, '===================================\n\n');

  if isfield(r, 'orbit_has_data') && r.orbit_has_data
    fprintf(fid, '--- Orbit report (fine FFT grid, spatial symmetries via FFT.R_rot) ---\n');
    fprintf(fid, '  New centroid count (with possible repeats): %d\n', int32(r.orbit_Nadded));
    fprintf(fid, '  Distinct new seeds (fine linear indices): %d\n', int32(r.orbit_seeds_count));
    fprintf(fid, '  (1) Distinct symmetry orbits among new seeds: %d\n', int32(r.orbit_n_orbits));
    fprintf(fid, '  (2) Total fine-grid points in union of those orbits: %d\n', int32(r.orbit_n_pts_union));
    fprintf(fid, '  (3a) Of those orbit points, in initial centroid set only (first Nisdf=%d): %d\n', ...
      int32(r.initial_nisdf), int32(r.orbit_n_in_initial_only));
    fprintf(fid, '  (3b) Of those orbit points, in full current set (initial + new seeds): %d\n', ...
      int32(r.orbit_n_in_current));
    fprintf(fid, '--- end orbit report ---\n\n');
  else
    fprintf(fid, '--- Orbit report: no new centroids (Nadded=0), skip. ---\n\n');
  end

  fprintf(fid, '=== end adaptiveisdf phase-1 report ===\n');
end

function Om = adaptiveisdf_orbit_mask_bfs(seed, R_rot, nr)
  Om = false(nr, 1);
  dq = seed;
  Om(seed) = true;
  head = 1;
  while head <= numel(dq)
    i = dq(head);
    head = head + 1;
    for is = 1:size(R_rot, 2)
      j = R_rot(i, is);
      if j >= 1 && j <= nr && ~Om(j)
        Om(j) = true;
        dq(end + 1) = j; %#ok<AGROW>
      end
    end
  end
end

function U = adaptiveisdf_orbit_union_bfs(seeds, R_rot, nr)
  U = false(nr, 1);
  for k = 1:numel(seeds)
    U = U | adaptiveisdf_orbit_mask_bfs(seeds(k), R_rot, nr);
  end
end

function lin = adaptiveisdf_r_sampling_rows_to_lin(isdf_data, fft_data, Nmu)
  ni = double(fft_data.fftgrid(:)).';
  Rgrid = double(fft_data.Rgrid_RLU);
  lin = zeros(Nmu, 1);
  Rs = double(isdf_data.R_sampling_RLU(1:Nmu, :));
  for i = 1:Nmu
    v = round(Rs(i, :));
    v = mod(v, ni);
    [is_hit, k] = ismember(v, Rgrid, 'rows');
    if ~is_hit
      dd = zeros(size(Rgrid, 1), 3);
      for ddim = 1:3
        t = abs(Rgrid(:, ddim) - v(ddim));
        dd(:, ddim) = min(t, min(abs(t - ni(ddim)), abs(t + ni(ddim))));
      end
      [~, k] = min(sum(dd, 2));
    end
    lin(i) = k;
  end
end 

function adaptiveisdf_print_run_config(id, isdf_data, params, threshold, num_add, ratio, ...
    adaptive_backend, adaptive_arithmetic)
  fprintf('\n=== Adaptive ISDF start ===\n');
  fprintf('Backend            : %s\n', adaptive_backend);
  fprintf('Arithmetic         : %s\n', adaptive_arithmetic);
  fprintf('ISDF desc          : %s\n', char(string(isdf_data.desc)));
  fprintf('Coarse ISDF id     : %d\n', int32(id));
  fprintf('Initial Nisdf      : %d\n', int32(isdf_data.nisdf));
  fprintf('Threshold          : %.8e\n', threshold);
  fprintf('num_add            : %d\n', int32(num_add));
  fprintf('candidate ratio    : %.4f\n', ratio);
  fprintf('ISDF ratio         : %.4f\n', params.isdf_ratio);
  fprintf('max cond           : %.4e\n', params.max_cond);
  fprintf('use cond guard     : %d\n', logical(params.use_cond_guard));
  fprintf('weight batch size  : %d\n', int32(params.weight_batch_size));
  fprintf('param source       : %s\n', params.source);
  fprintf('===========================\n\n');
end

function params = adaptiveisdf_resolve_params(desc, cfg_isdf)
  desc = lower(strtrim(char(string(desc))));
  params = struct();
  params.threshold = 2e-4;
  params.num_add = int32(16);
  params.candidate_ratio = 2.0;
  params.max_add_frac = 1.0;
  params.isdf_ratio = 8.0;
  params.max_cond = 1e6;
  params.use_cond_guard = true;
  params.weight_batch_size = 256;
  params.source = "legacy";

  if isempty(cfg_isdf) || ~isstruct(cfg_isdf)
    return;
  end

  if strcmp(desc, 'vc')
    suffix = 'type1';
  elseif strcmp(desc, 'vn')
    suffix = 'type2';
  elseif strcmp(desc, 'nn')
    suffix = 'type3';
  else
    return;
  end

  params.threshold = adaptiveisdf_get_cfg_positive(cfg_isdf, ['adaptive_threshold_' suffix], params.threshold);
  num_add_d = adaptiveisdf_get_cfg_integer(cfg_isdf, ['adaptive_num_add_' suffix], double(params.num_add));
  params.num_add = int32(num_add_d);
  params.candidate_ratio = adaptiveisdf_get_cfg_positive(cfg_isdf, ['adaptive_candidate_ratio_' suffix], params.candidate_ratio);
  params.max_add_frac = adaptiveisdf_get_cfg_positive(cfg_isdf, ['adaptive_max_add_frac_' suffix], params.max_add_frac);
  params.isdf_ratio = adaptiveisdf_get_cfg_positive(cfg_isdf, ['isdf_ratio_' suffix], params.isdf_ratio);
  params.max_cond = adaptiveisdf_get_cfg_positive(cfg_isdf, ['adaptive_max_cond_' suffix], params.max_cond);
  params.use_cond_guard = adaptiveisdf_get_cfg_logical(cfg_isdf, 'adaptive_use_cond_guard', params.use_cond_guard);
  params.weight_batch_size = adaptiveisdf_get_cfg_integer(cfg_isdf, 'adaptive_weight_batch_size', params.weight_batch_size);
  params.weight_batch_size = adaptiveisdf_get_cfg_integer(cfg_isdf, 'adaptive_batch_size', params.weight_batch_size);
  params.source = "config";
end

function val = adaptiveisdf_get_cfg_positive(cfg, field_name, fallback)
  val = fallback;
  if ~isfield(cfg, field_name)
    return;
  end
  v = double(cfg.(field_name));
  if isfinite(v) && v > 0
    val = v;
  end
end

function val = adaptiveisdf_get_cfg_integer(cfg, field_name, fallback)
  val = fallback;
  if ~isfield(cfg, field_name)
    return;
  end
  v = double(cfg.(field_name));
  if isfinite(v) && v >= 1
    val = max(1, round(v));
  end
end

function val = adaptiveisdf_get_cfg_logical(cfg, field_name, fallback)
  val = fallback;
  if ~isfield(cfg, field_name)
    return;
  end
  v = cfg.(field_name);
  if islogical(v)
    val = logical(v);
    return;
  end
  if isnumeric(v) && isfinite(v) && isscalar(v)
    val = logical(v ~= 0);
    return;
  end
  if ischar(v) || isstring(v)
    s = lower(strtrim(char(string(v))));
    if any(strcmp(s, {'true', '.true.', '1', 'yes', 'on'}))
      val = true;
    elseif any(strcmp(s, {'false', '.false.', '0', 'no', 'off'}))
      val = false;
    end
  end
end

function root = local_get_import_isdf_root(cfg_isdf)
  root = '';
  if ~isstruct(cfg_isdf) || ~isfield(cfg_isdf, 'import_isdf_adaptive_root')
    return;
  end
  v = cfg_isdf.import_isdf_adaptive_root;
  if isempty(v)
    return;
  end
  if isstring(v) && isscalar(v)
    v = char(v);
  end
  if ~(ischar(v) && isvector(v) && ~isempty(strtrim(v)))
    return;
  end
  root = char(strtrim(v));
  if exist(root, 'file') ~= 2 && exist(root, 'dir') ~= 7
    fprintf(1, 'adaptiveisdf: import_isdf_adaptive_root not found: %s\n', root);
    root = '';
  end
end
