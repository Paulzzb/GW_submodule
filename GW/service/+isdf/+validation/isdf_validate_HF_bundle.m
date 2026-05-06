% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function [Esum2, EsumISDF2, DiffEsum2, report] = isdf_validate_HF(id, outDir)
% ISDF_COARSE_VALIDATE_ENERGIES  Compare direct Coulomb exchange-style energy vs ISDF tildeVq contraction.
%
% Text report (band summary, preview, stats, global sums) is written by isdf.report.gen_report to
%   <outDir>/isdf_validate_HF_id<id>.txt  (outDir defaults to pwd; see report.report_file when nargout>3).
% A one-line message with the absolute path is printed to the command window after a successful write.
%
% Mirrors the accumulation loops in gen_indices_coarse_test (SCATTER_Bamp vs c_rho' * tildeVq * c_rho).
% Optional trailing args:
%   One matrix R_sampling_RLU (Nwf x 3): same phase as gen_tildeVq, row-aligned with wf_on_coarse.
%   Or {R_coarse_RLU, fftgrid_i, fftgrid_c}: integer coarse box scaled by fftgrid_i./fftgrid_c (isdftest path).
%
% Progress and elapsed / estimated total time use timing.LIVE_timing (CPU time via timing.timing_string).
% tildeVq, wf_on_coarse, R_rot_coarse, varargin

  if nargin < 2
    outDir = '';
  end

  wf_data = wave_functions.get();
  k_data = lattice.manager('k', 'get');
  r_lat_data = lattice.manager('r_lat', 'get');
  coul_data = coulomb.get();
  system_data = system.get();
  fft_data = FFT.get();
  symm_data = symmetry.get();
  nsym = double(symm_data.nsym);

  nibz = int32(k_data.nibz);
  nbz = int32(k_data.nbz);
  nb = int32(wf_data.nb);
  nspin = int32(wf_data.nspin);

  


  isdf_data = isdf.get(id);
  R_sampling_RLU = isdf_data.R_sampling_RLU;
  
  Esum2 = 0.0;
  EsumISDF2 = 0.0;
  DiffEsum2 = 0.0;
  Esum = 0.0;
  DiffEsum = 0.0;

  % Online stats (Welford): O(1) memory; no full sample buffer.
  st_Ex_t = isdf_validate_HF_stream_init();
  st_Ex_ISDF = isdf_validate_HF_stream_init();
  st_Diff = isdf_validate_HF_stream_init();
  st_AbsDiff = isdf_validate_HF_stream_init();
  st_RelSigned = isdf_validate_HF_stream_init();
  st_RelAbs = isdf_validate_HF_stream_init();

  % Optional preview only: O(nb^2 * nibz * nspin), not O(nibz * nbz * nb^2).
  preview_budget = double(nb) * double(nb) * double(nibz) * double(nspin);
  preview_rows = min(max(preview_budget, 1), 1);
  prev_Ex_t = nan(preview_rows, 1);
  prev_Ex_ISDF = nan(preview_rows, 1);
  n_ob_sample = 0;

  total_triples = max(1, double(nibz) * double(nspin) * double(nb));

  isdftest_ensure_timing_initialized();
  tm_live = timing.get();
  tm_live.live.nhash = int32(20);
  tm_live.live.live_report_min_seconds = 0;
  timing.save2mod(tm_live);

  timing.LIVE_timing('isdf validate (ik,ispin,ib)', total_triples);
  cleanup_live = onCleanup(@() timing.LIVE_timing());

  E_HF = zeros(nb, nibz, nspin);
  E_HF_ISDF = zeros(nb, nibz, nspin);
  ibnb = 1;
  obnb = 4;
  if strcmp(isdf_data.interp_scheme, "adaptive") || strcmp(isdf_data.interp_scheme, "coarse")
    N_coarse = isdf_data.N_coarse;
    N_extra = isdf_data.N_extra;
  end
  for ik = 1:nibz
    ikibz = ik;
    ikbz = k_data.ibz2bz(ikibz);
    ikrot = 1;
    for ispin = 1:nspin
      for ib = 1:ibnb
        for iqbz = 1:nbz
          iqibz = k_data.bz2ibz(iqbz, 1);
          iqrot = k_data.bz2rot(iqbz, 1);
          inviqrot = symm_data.inv_rot_index(iqrot);

          ikpbz = r_lat_data.qindx_S(ik, iqbz, 1);
          iGo = r_lat_data.qindx_S(ik, iqbz, 2);
          ikpibz = k_data.bz2ibz(ikpbz, 1);
          ikprot = k_data.bz2rot(ikpbz, 1);

          isc = [ib, ik, 1, ispin];


          vcoul_q = coul_data.vcoul(:, iqibz);
          if iqibz == 1
            vcoul_q(1) = coul_data.vcoul0;
          end


          for ob = 1:obnb
            occupation = system_data.f(ob, ikpibz, ispin);
            if occupation < 1e-6
              continue;
            end


            % Can we get Psi_{S_qibz^{-1}*k_1}(x_alpha) instead?
            % We need <S_qibz * x_alpha |ib, ikibz, ob, ikpbz >,
            % where 
            %    -  {x_alpha} is the isdf grid points,
            %    -   S_qibz is the symmetry operator that q_bz = q_ibz * S_qibz
            %    -   |ib, ikibz, ob, ikpbz > is the indices
            % which means
            % 1. Get indices ib, ikibz, ob, ikpbz, ispin
            if strcmp(isdf_data.interp_scheme, "coarse")
              % u_{ib, ikibz}(r)
              wf_1_c = isdf_data.coeff_seper(:, ib, ikibz, ispin); 
              % u_{ib, ikbz}(r) = u_{ib, ikibz}(S_ikrot^{-1} * r)
              if ikrot ~= 1
                wf_1_c = isdf.coeff.isdf_apply_symm_on_coarse(id, wf_1_c, ikrot);
              end
              % apply rotation of iqrot
              % u_{ib, ikbz}(S_iqrot * r)
              if iqrot ~= 1
                % wf_1_c = isdf.coeff.isdf_apply_symm_on_coarse(id, wf_1_c, iqrot);
                ind = isdf_data.R_rot_extra(1:N_coarse, iqrot);
                wf_1_c = wf_1_c(ind);
              end
            elseif strcmp(isdf_data.interp_scheme, "adaptive")
              ikbz = double(k_data.ibz2bz(ikibz));
              wf_1_c = isdf_data.coeff_seper(:, ib, ikbz, ispin);
              % This is u_{ib, ikbz}(r)
              if iqrot ~= 1
                % Update u_{ib, ikbz}(S_iqrot * r) for coarse grid points
                ind = isdf_data.R_rot_extra(1:N_coarse, iqrot);
                wf_1_c(1:N_coarse) = wf_1_c(ind);
                % Update u_{ib, ikbz}(S_iqrot^{-1} * r) for extra grid points
                isc1 = int32([ib, ikibz, ikrot, 1]);
                wf_c_1_tmp = wave_functions.WF_apply_symm(isc1);
                % inv_isym = symm_data.inv_rot_index(iqrot);
                % ind = isdf_data.R_rot_extra(inv_isym, :);
                ind = isdf_data.R_rot_extra(N_coarse + 1:N_coarse + N_extra, iqrot);
                wf_1_c(N_coarse + 1:N_coarse + N_extra) = wf_c_1_tmp(ind);
              end
            end
            % Similar for u_{ob, ikpbz}(S_{iqrot} * r)
            if strcmp(isdf_data.interp_scheme, "coarse")
              wf_2_c = isdf_data.coeff_seper(:, ob, ikpibz, ispin);
              if ikprot ~= 1
                wf_2_c = isdf.coeff.isdf_apply_symm_on_coarse(id, wf_2_c, ikprot);
              end
              if iqrot ~= 1
                ind = isdf_data.R_rot_extra(1:N_coarse, iqrot);
                wf_2_c = wf_2_c(ind);
              end
            elseif strcmp(isdf_data.interp_scheme, "adaptive")
              % This is u_{ob, ikpbz}(r)
              wf_2_c = isdf_data.coeff_seper(:, ob, ikpbz, ispin);
              if iqrot ~= 1
                % Update u_{ob, ikpbz}(S_{iqrot} * r) for coarse grid points
                ind = isdf_data.R_rot_extra(1:N_coarse, iqrot);
                wf_2_c(1:N_coarse) = wf_2_c(ind);
                % Update u_{ob, ikpbz}(S_{iqrot} * r) for extra grid points
                isc2 = int32([ob, ikpibz, ikprot, 1]);
                wf_c_2_tmp = wave_functions.WF_apply_symm(isc2);
                ind = isdf_data.R_rot_extra(N_coarse + 1:N_coarse + N_extra, iqrot);
                wf_2_c(N_coarse + 1:N_coarse + N_extra) = wf_c_2_tmp(ind);
              end
            end
            % c_rho = isdf.coeff.isdf_get_coeff(id, wf_1_c, wf_2_c, ikrot, ikprot);
            c_rho = conj(wf_1_c) .* wf_2_c;
            
            Go = single(r_lat_data.Ggrid_RLU(iGo, :));
            % SqR_scal = single(R_sampling_RLU*symm_data.rot_mtrx_RLU_R(:, :, iqrot)) ./ double(fft_data.fftgrid);
            invSqGo = single(Go * symm_data.rot_mtrx_RLU_G(:, :, inviqrot));
            phase_shift_coarse = exp(-2 * pi * 1i * (R_sampling_RLU ./ double(fft_data.fftgrid)) * invSqGo');
            % phase_shift_coarse = exp(-2 * pi * 1i * SqR_scal * Go');
            c_rho = c_rho .* phase_shift_coarse;

            iscp = [ob, ikpibz, ikprot, ispin];
            param = [];
            param.is = isc;
            param.os = iscp;
            param.qs = [iGo, iqibz, iqrot];
            ngrho_left = SCATTER_Bamp(param);

            Ex_t = sum(vcoul_q .* abs(ngrho_left).^2);
            Ex_ISDF = c_rho' * isdf_data.tildeVq(:, :, iqibz) * c_rho;
            Ex_ISDF = real(Ex_ISDF);
            if abs(Ex_ISDF - Ex_t) > 1e-4
              warning('isdf_validate_HF: Ex_ISDF - Ex_t = %f', Ex_ISDF - Ex_t);
            end
            E_HF(ib, ik, ispin) = E_HF(ib, ik, ispin) + Ex_t;
            E_HF_ISDF(ib, ik, ispin) = E_HF_ISDF(ib, ik, ispin) + Ex_ISDF;
            n_ob_sample = n_ob_sample + 1;
            diff_v = Ex_t - Ex_ISDF;
            abs_diff_v = abs(diff_v);
            thres = 5e-3;
            if abs(Ex_t) >= thres
              st_Ex_t = isdf_validate_HF_stream_push(st_Ex_t, Ex_t);
              st_Ex_ISDF = isdf_validate_HF_stream_push(st_Ex_ISDF, Ex_ISDF);
              st_Diff = isdf_validate_HF_stream_push(st_Diff, diff_v);
              st_AbsDiff = isdf_validate_HF_stream_push(st_AbsDiff, abs_diff_v);
              st_RelSigned = isdf_validate_HF_stream_push(st_RelSigned, diff_v / abs(Ex_t));
              st_RelAbs = isdf_validate_HF_stream_push(st_RelAbs, abs_diff_v / abs(Ex_t));
            end
            if n_ob_sample <= preview_rows
              prev_Ex_t(n_ob_sample) = Ex_t;
              prev_Ex_ISDF(n_ob_sample) = Ex_ISDF;
            end
            Esum2 = Esum2 + Ex_t^2;
            EsumISDF2 = EsumISDF2 + Ex_ISDF^2;
            DiffEsum2 = DiffEsum2 + (Ex_t - Ex_ISDF)^2;
            Esum = Esum + Ex_t + Ex_ISDF;
            DiffEsum = DiffEsum + abs(Ex_t - Ex_ISDF);
          end
        end % iqbz

        timing.LIVE_timing(1);
      end % ib
    end % ispin
  end % ik

  % --- Band-resolved report: E_HF(ib,ik,ispin) vs E_HF_ISDF (sums over ob paths) ---
  sum_E_HF = sum(E_HF(:));
  sum_E_ISDF = sum(E_HF_ISDF(:));
  diff_E = E_HF - E_HF_ISDF;
  sum_abs_diff_E = sum(abs(diff_E(:)));
  max_abs_diff_E = max(abs(diff_E(:)));
  frob_hf = norm(E_HF(:), 2);
  frob_diff = norm(diff_E(:), 2);
  ncells = double(nb) * double(nibz) * double(nspin);
  mean_E_HF = sum_E_HF / ncells;
  mean_E_ISDF = sum_E_ISDF / ncells;
  [max_E_HF, lin_mx_hf] = max(E_HF(:));
  [max_E_ISDF, lin_mx_isdf] = max(E_HF_ISDF(:));
  [ib_mx_hf, ik_mx_hf, is_mx_hf] = ind2sub([double(nb), double(nibz), double(nspin)], lin_mx_hf);
  [ib_mx_isdf, ik_mx_isdf, is_mx_isdf] = ind2sub([double(nb), double(nibz), double(nspin)], lin_mx_isdf);
  [min_E_HF, lin_mn_hf] = min(E_HF(:));
  [min_E_ISDF, lin_mn_isdf] = min(E_HF_ISDF(:));
  [ib_mn_hf, ik_mn_hf, is_mn_hf] = ind2sub([double(nb), double(nibz), double(nspin)], lin_mn_hf);
  [ib_mn_isdf, ik_mn_isdf, is_mn_isdf] = ind2sub([double(nb), double(nibz), double(nspin)], lin_mn_isdf);

  % Band-resolved summary: keep the table rows short so disp(table) aligns in the command window.
  row_label = strings(0, 1);
  row_result = strings(0, 1);
  row_label(end + 1, 1) = "sum(E_HF)";
  row_result(end + 1, 1) = string(sprintf('%.8e', sum_E_HF));
  row_label(end + 1, 1) = "sum(E_HF_ISDF)";
  row_result(end + 1, 1) = string(sprintf('%.8e', sum_E_ISDF));
  row_label(end + 1, 1) = "sum(E_HF) - sum(E_HF_ISDF)";
  row_result(end + 1, 1) = string(sprintf('%.8e', sum_E_HF - sum_E_ISDF));
  row_label(end + 1, 1) = "mean(E_HF) (over nb*nibz*nspin cells)";
  row_result(end+1, 1) = string(sprintf('%.8e', mean_E_HF));
  row_label(end+1, 1) = "mean(E_HF_ISDF) (over nb*nibz*nspin cells)";
  row_result(end+1, 1) = string(sprintf('%.8e', mean_E_ISDF));
  row_label(end+1, 1) = "sum(|E_HF - E_HF_ISDF|)";
  row_result(end+1, 1) = string(sprintf('%.8e', sum_abs_diff_E));
  row_label(end+1, 1) = "max(|E_HF - E_HF_ISDF|)";
  row_result(end+1, 1) = string(sprintf('%.8e', max_abs_diff_E));
  row_label(end+1, 1) = "||vec(E_HF - E_HF_ISDF)||_2";
  row_result(end+1, 1) = string(sprintf('%.8e', frob_diff));
  row_label(end+1, 1) = "||vec(E_HF - E_HF_ISDF)||_2 / ||vec(E_HF)||_2";
  row_result(end+1, 1) = string(sprintf('%.8e', frob_diff / max(frob_hf, eps('double'))));
  row_label(end+1, 1) = "E_HF minimum";
  row_result(end+1, 1) = string(sprintf('%.8e at (ib=%d, ik=%d, ispin=%d)', min_E_HF, ib_mn_hf, ik_mn_hf, is_mn_hf));
  row_label(end+1, 1) = "E_HF maximum";
  row_result(end+1, 1) = string(sprintf('%.8e at (ib=%d, ik=%d, ispin=%d)', max_E_HF, ib_mx_hf, ik_mx_hf, is_mx_hf));
  row_label(end+1, 1) = "E_HF_ISDF minimum";
  row_result(end+1, 1) = string(sprintf('%.8e at (ib=%d, ik=%d, ispin=%d)', min_E_ISDF, ib_mn_isdf, ik_mn_isdf, is_mn_isdf));
  row_label(end+1, 1) = "E_HF_ISDF maximum";
  row_result(end+1, 1) = string(sprintf('%.8e at (ib=%d, ik=%d, ispin=%d)', max_E_ISDF, ib_mx_isdf, ik_mx_isdf, is_mx_isdf));
  if nspin > 1
    for is = 1:double(nspin)
      s_hf = sum(sum(E_HF(:, :, is)));
      s_is = sum(sum(E_HF_ISDF(:, :, is)));
      s_ad = sum(sum(abs(E_HF(:, :, is) - E_HF_ISDF(:, :, is))));
      row_label(end+1, 1) = string(sprintf('ispin=%d: sum(E_HF)', is));
      row_result(end+1, 1) = string(sprintf('%.8e', s_hf));
      row_label(end+1, 1) = string(sprintf('ispin=%d: sum(E_HF_ISDF)', is));
      row_result(end+1, 1) = string(sprintf('%.8e', s_is));
      row_label(end+1, 1) = string(sprintf('ispin=%d: sum(|E_HF - E_HF_ISDF|)', is));
      row_result(end+1, 1) = string(sprintf('%.8e', s_ad));
    end
  end
  band_summary_tbl = table(row_label, row_result, 'VariableNames', {'Item', 'Result'});

  n_rel_samples = st_RelSigned.n;

  vars = {'Ex_t', 'Ex_ISDF', 'Diff', 'AbsDiff', 'Diff_over_abs_Ex_t', 'AbsDiff_over_abs_Ex_t'};
  streams = {st_Ex_t, st_Ex_ISDF, st_Diff, st_AbsDiff, st_RelSigned, st_RelAbs};
  stat_names = {'Mean'; 'Std'; 'Var'; 'Max'};
  M = zeros(4, numel(vars));
  for j = 1:numel(vars)
    [M(1, j), M(2, j), M(3, j), M(4, j)] = isdf_validate_HF_stream_finalize(streams{j});
  end
  stats_tbl = array2table(M, 'VariableNames', vars, 'RowNames', stat_names);

  nshow = min(preview_rows, n_ob_sample);
  if nshow > 0
    diff_p = prev_Ex_t(1:nshow) - prev_Ex_ISDF(1:nshow);
    abs_dp = abs(diff_p);
    tol_zero = max(eps('double'), 1e-30);
    rel_s = nan(nshow, 1);
    rel_a = nan(nshow, 1);
    for ii = 1:nshow
      abs_ex = abs(prev_Ex_t(ii));
      if abs_ex >= tol_zero
        rel_s(ii) = diff_p(ii) / abs_ex;
        rel_a(ii) = abs_dp(ii) / abs_ex;
      end
    end
    preview_tbl = table(prev_Ex_t(1:nshow), prev_Ex_ISDF(1:nshow), diff_p, abs_dp, rel_s, rel_a, ...
      'VariableNames', {'Ex_t', 'Ex_ISDF', 'Diff', 'AbsDiff', 'Diff_over_abs_Ex_t', 'AbsDiff_over_abs_Ex_t'});
  else
    preview_tbl = table();
  end

  report = struct();
  report.isdf_id = id;
  report.desc = isdf_data.desc;
  report.interp_scheme = isdf_data.interp_scheme;
  report.nisdf = isdf_data.nisdf;
  report.band_summary = band_summary_tbl;
  report.preview = preview_tbl;
  report.stats = stats_tbl;
  report.n_samples = n_ob_sample;
  report.n_samples_rel = n_rel_samples;
  report.preview_rows = preview_rows;
  report.Esum2 = Esum2;
  report.EsumISDF2 = EsumISDF2;
  report.DiffEsum2 = DiffEsum2;
  report.E_HF = E_HF;
  report.E_HF_ISDF = E_HF_ISDF;
  report.E_diff = diff_E;
  report.sum_E_HF = sum_E_HF;
  report.sum_E_HF_ISDF = sum_E_ISDF;
  report.sum_abs_diff_E = sum_abs_diff_E;
  report.max_abs_diff_E = max_abs_diff_E;
  report.l2_diff_over_l2_HF = frob_diff / max(frob_hf, eps('double'));
  report.mean_E_HF = mean_E_HF;
  report.mean_E_HF_ISDF = mean_E_ISDF;
  report.max_E_HF = struct('value', max_E_HF, 'ib', ib_mx_hf, 'ik', ik_mx_hf, 'ispin', is_mx_hf);
  report.max_E_HF_ISDF = struct('value', max_E_ISDF, 'ib', ib_mx_isdf, 'ik', ik_mx_isdf, 'ispin', is_mx_isdf);
  report.min_E_HF = struct('value', min_E_HF, 'ib', ib_mn_hf, 'ik', ik_mn_hf, 'ispin', is_mn_hf);
  report.min_E_HF_ISDF = struct('value', min_E_ISDF, 'ib', ib_mn_isdf, 'ik', ik_mn_isdf, 'ispin', is_mn_isdf);

  fpath = isdf.report.gen_report(report, outDir);
  report.report_file = fpath;
  fprintf(1, 'isdf_validate_HF: wrote HF report to:\n  %s\n', fpath);

  if nargout <= 3
    report = []; %#ok<NASGU>
  end
end

function s = isdf_validate_HF_stream_init()
  s = struct('n', 0, 'mean', 0, 'M2', 0, 'mx', NaN, 'have_max', false);
end

function s = isdf_validate_HF_stream_push(s, x)
  if ~isfinite(x)
    return;
  end
  s.n = s.n + 1;
  delta = x - s.mean;
  s.mean = s.mean + delta / s.n;
  delta2 = x - s.mean;
  s.M2 = s.M2 + delta * delta2;
  if ~s.have_max
    s.mx = x;
    s.have_max = true;
  else
    s.mx = max(s.mx, x);
  end
end

function [mn, st, vr, mx] = isdf_validate_HF_stream_finalize(s)
  if s.n < 1
    mn = NaN;
    st = NaN;
    vr = NaN;
    mx = NaN;
  elseif s.n == 1
    mn = s.mean;
    st = 0;
    vr = 0;
    mx = s.mx;
  else
    mn = s.mean;
    vr = s.M2 / (s.n - 1);
    st = sqrt(vr);
    mx = s.mx;
  end
end

function isdftest_ensure_timing_initialized()
  try
    timing.get();
  catch %#ok<CTCH>
    timing.driver();
  end
end
