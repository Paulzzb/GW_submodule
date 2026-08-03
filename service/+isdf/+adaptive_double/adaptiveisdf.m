% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/20 ZZ

function idnew = adaptiveisdf(id, cfg_isdf)
  % Wall-clock time through rebuild + orbit summary + phase-1 report.
  t_phase1 = tic;

  if nargin < 2 || isempty(cfg_isdf) || ~isstruct(cfg_isdf)
    error('adaptiveisdf:cfg', ...
      ['cfg_isdf (config.ISDF) is required. ', ...
       'Initialize once via default_param_values / set_default_param_value.']);
  end

  % 1. Starting grid: from coarse fft grid
  system_data = system.get();
  nspin = system_data.nspin;
  k_data = lattice.manager('k', 'get');
  if nspin > 1
    error('adaptiveisdf:nspin', 'nspin > 1 is not supported.');
  end

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
  params = isdf.adaptive.adaptive_param(isdf_data.desc, cfg_isdf);
  threshold = params.threshold;
  num_add = params.num_add;
  ratio = params.candidate_ratio;
  adaptive_backend = 'adaptive_double';
  adaptive_arithmetic = 'double';
  isdf.adaptive.adaptiveisdf_print_run_config(id, isdf_data, params, threshold, num_add, ratio, ...
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
  output.msg('v2l', 'adaptiveisdf: initial loss = %.8e', loss0);
  isdf.adaptive_double.adaptive_weight('init_update');
  w = isdf.adaptive_double.adaptive_weight('get');
  loss = sum(w);
  output.msg('v2l', 'adaptiveisdf: loss after initial update = %.8e', loss);
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
  n_iter = 0;
  n_schur_skip = 0;
  schur_skip_warned = false;
  loss_history = zeros(double(Naddmax) + 1, 1);
  rel_loss_history = zeros(double(Naddmax) + 1, 1);
  loss_history(1) = loss;
  if loss0 > 0
    rel_loss_history(1) = loss / loss0;
  else
    rel_loss_history(1) = NaN;
  end

  live_on = false;
  if double(Nisdf_new_max) > 0
    try
      timing.get();
    catch %#ok<CTCH>
      timing.driver();
    end
    tm_live = timing.get();
    tm_live.live.nhash = int32(20);
    tm_live.live.live_report_min_seconds = 0;
    tm_live.live.show_expected = false;  % non-linear cost: no continuous (X)
    timing.save2mod(tm_live);
    if loss0 > 0
      live_label0 = sprintf('Adaptive ISDF rel=%.3e', loss / loss0);
    else
      live_label0 = 'Adaptive ISDF rel=NaN';
    end
    timing.LIVE_timing(live_label0, double(Nisdf_new_max));
    live_on = true;
    cleanup_live = onCleanup(@() timing.LIVE_timing()); %#ok<NASGU>
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
          n_schur_skip = n_schur_skip + 1;
          output.msg('v2l', ...
            'adaptiveisdf: skip index %d (Schur guard failed), relative loss = %.6e', ...
            selected_indices(iadd), rel_now);
          if ~schur_skip_warned && double(Nisdf_new_max) > 0 ...
              && n_schur_skip > 0.5 * double(Nisdf_new_max)
            output.warn([ ...
              'adaptiveisdf: schur_skips=%d > 0.5*Naddmax=%d - method may be stuck. ', ...
              'Try: (1) decrease adaptive_batch_size  (2) set exxmethod=''pseudo''  then recompute.'], ...
              n_schur_skip, double(Nisdf_new_max));
            schur_skip_warned = true;
          end
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
      output.msg('rs', ...
        'adaptiveisdf: all candidate Schur updates failed this iteration; treating as converged.');
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
    output.msg('v2l', 'adaptiveisdf: iter=%d  loss=%.8e  rel=%.8e  added=%d/%d', ...
      n_iter, loss, rel_loss_history(n_iter + 1), Nisdf_new, Nisdf_new_max);
    if live_on
      label = sprintf('Adaptive ISDF rel=%.3e', rel_loss_history(n_iter + 1));
      timing.LIVE_timing(double(numel(selected_now)), label);
    end
    if breakflag || schur_converged
      break;
    end
  end

  if live_on
    timing.LIVE_timing();
    live_on = false;
    clear cleanup_live
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Verification
  [Nisdf_new, ~, ~, ~, Psi_on_grid, Phi_on_grid] = ...
         isdf.adaptive_double.isdf_schur_update('get');
  Nextra = Nisdf_new - Nisdf;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Build a new ISDF object at idnew, keep original id unchanged.
  idnew = isdf.isdf_add(isdf_data.desc);

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

  % Append newly selected R points. (Symmetry for extras lives in bundle_struct after refresh.)
  isdf_data_new.R_rot_extra = int32(zeros(0, 0));
  if Nextra > 0
    new_lin = int32(isdf_new_indices(1:Nextra));
    R_new = double(fft_data.Rgrid_RLU(new_lin, :));
    isdf_data_new.R_sampling_RLU = double([isdf_data.R_sampling_RLU; R_new]);
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
  isdf_data_new.R_rot_coarse = isdf_data.R_rot_coarse;
  isdf_data_new.bundle_struct = isdf_data.bundle_struct;
  isdf_data_new.assigned = true;
  isdf.save2mod(isdf_data_new, idnew);
  isdf.rsymm.bundle_refresh(id, Nextra, isdf_new_indices, idnew);
  isdf.adaptive_double.adaptive_checkpoint_save(id, idnew);





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

  report = struct();
  report.iterations = n_iter;
  report.loss_initial = loss0;
  report.loss_final = loss;
  report.relative_loss = rel_loss;
  report.stop_reason = stop_reason;
  report.adaptive_backend = adaptive_backend;
  report.adaptive_arithmetic = adaptive_arithmetic;
  report.loss_history = loss_history(1:n_iter+1);
  report.relative_loss_history = rel_loss_history(1:n_iter+1);
  report.schur_skips = n_schur_skip;
  report = isdf.adaptive.adaptive_fill_report(report, id, idnew, params);

  elapsed_phase1 = toc(t_phase1);
  report.elapsed_phase1_seconds = elapsed_phase1;
  fpath_r = isdf.adaptive.adaptiveisdf_write_phase1_report(report);

  output.msg('nrs', '[Adaptive ISDF] desc=%s  id=%d  backend=%s', ...
    char(string(isdf_data.desc)), int32(id), adaptive_backend);
  output.msg('rs', '  Nisdf %d -> %d  (+%d)  iters=%d  stop=%s', ...
    Nisdf, report.final_nisdf, Nextra, n_iter, stop_reason);
  output.msg('rs', '  loss0=%.6e  loss=%.6e  rel=%.6e  (thr=%.6e)', ...
    loss0, loss, rel_loss, threshold);
  output.msg('rs', '  schur_skips=%d  (details at verbose>=2 / log)', n_schur_skip);
  output.msg('rs', '  wrote %s  (%.1f s)', fpath_r, elapsed_phase1);

  isdf.numerical_cond_report('adaptive', id, idnew, isdf_data.desc, ...
    loss_history(1), loss, report.final_nisdf);
  output.msg('rs', '  Saved adaptive ISDF -> id=%d', idnew);

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
    output.msg('rs', 'adaptiveisdf: import_isdf_adaptive_root not found: %s', root);
    root = '';
  end
end
