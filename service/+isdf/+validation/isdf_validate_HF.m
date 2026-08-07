% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08

function [Esum2, EsumISDF2, DiffEsum2, report] = isdf_validate_HF(id)
% ISDF_COARSE_VALIDATE_ENERGIES  Compare direct Coulomb exchange-style energy vs ISDF tildeVq contraction.
%
% Text report is written by isdf.report.hf under filename_map().isdf_report_dir
% (see +report/NAMING.md).
% A one-line message with the absolute path is printed to the command window after a successful write.
%
% Progress and elapsed / estimated total time use timing.LIVE_timing (CPU time via timing.timing_string).

  wf_data = wave_functions.get();
  k_data = lattice.manager('k', 'get');
  r_lat_data = lattice.manager('r_lat', 'get');
  coul_data = coulomb.get();
  system_data = system.get();

  nibz = int32(k_data.nibz);
  nbz = int32(k_data.nbz);
  nb = int32(wf_data.nb);
  nspin = int32(wf_data.nspin);
  isdf_data = isdf.get(id);

  % nv rule:
  % first band where, at the last index of 2nd dimension (k-index),
  % there exists at least one occupation-flag == false.
  tol_occ = 1e-6;
  occ_flag = double(system_data.f) < tol_occ;
  for i = 2:nb
    if all(occ_flag(i, :, :))
      nv = i-1;
      break;
    end
  end

  cfg_ibmax = isdf_data.nrange2(end);
  cfg_ibmin = isdf_data.nrange2(1);

  ob_list = int32(1:min(double(nv), double(nb)));
  if strcmpi(char(isdf_data.desc), 'vc')
    ib_start = double(nv) + 1;
    ib_end = 2 * double(nv);
  else
    ib_start = cfg_ibmin;
    ib_end = cfg_ibmax;
  end
  ib_start = max(1, ib_start);
  ib_end = min(double(nb), ib_end);
  if ib_end < ib_start
    ib_list = int32([]);
  else
    ib_list = int32(ib_start:ib_end);
  end
  
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

  n_ob_sample = 0;
  diff_warn_tol = 1e-4;
  n_mismatch = 0;
  max_abs_mismatch = 0.0;

  total_triples = max(1, double(nibz) * double(nspin) * double(numel(ib_list)));

  isdf_ensure_timing_initialized();
  tm_live = timing.get();
  tm_live.live.nhash = int32(20);
  tm_live.live.live_report_min_seconds = 0;
  timing.save2mod(tm_live);

  timing.LIVE_timing('isdf validate (ik,ispin,ib)', total_triples);
  cleanup_live = onCleanup(@() timing.LIVE_timing());

  E_HF = zeros(nb, nibz, nspin);
  E_HF_ISDF = zeros(nb, nibz, nspin);

  for ik = 1:nibz
    ikibz = ik;
    ikrot = 1;
    for ispin = 1:nspin
      for ib = reshape(ib_list, 1, [])
        for iqbz = 1:nbz
          iqibz = k_data.bz2ibz(iqbz, 1);
          iqrot = k_data.bz2rot(iqbz, 1);
          ikpbz = r_lat_data.qindx_S(ik, iqbz, 1);
          iGo = r_lat_data.qindx_S(ik, iqbz, 2);
          ikpibz = k_data.bz2ibz(ikpbz, 1);
          ikprot = k_data.bz2rot(ikpbz, 1);

          isc = [ib, ik, 1, ispin];


          vcoul_q = coul_data.vcoul(:, iqibz);
          if iqibz == 1
            vcoul_q(1) = coul_data.vcoul0;
          end
          fac_q = isdf_data.CCHq_trunc_factors{iqibz};
          N_keep_q = fac_q.N_keep;

          for ob = reshape(ob_list, 1, [])
            occupation = system_data.f(ob, ikpibz, ispin);
            if occupation < 1e-6
              continue;
            end

            iscp = [ob, ikpibz, ikprot, ispin];
            param = [];
            param.is = isc;
            param.os = iscp;
            param.qs = [iGo, iqibz, iqrot];
            ngrho_left = SCATTER_Bamp(param);
            c_rho = isdf.get_rho_xalpha(id, param);

            Ex_t = sum(vcoul_q .* abs(ngrho_left).^2);
            Ex_ISDF = c_rho(1:N_keep_q)' * isdf_data.tildeVq(1:N_keep_q, 1:N_keep_q, iqibz) * c_rho(1:N_keep_q);
            Ex_ISDF = real(Ex_ISDF);
            d_isdf_t = Ex_ISDF - Ex_t;
            ad_isdf_t = abs(d_isdf_t);
            if ad_isdf_t > diff_warn_tol
              n_mismatch = n_mismatch + 1;
              if ad_isdf_t > max_abs_mismatch
                max_abs_mismatch = ad_isdf_t;
              end
              output.msg('v2l', ...
                ['isdf_validate_HF: |Ex_ISDF-Ex_t|=%.6e (Ex_ISDF-Ex_t=%+.6e) ', ...
                 'ib=%d ik=%d iqibz=%d ob=%d'], ...
                ad_isdf_t, d_isdf_t, ib, ik, iqibz, ob);
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

  output.msg('rs', ...
    'isdf_validate_HF: mismatches(|Ex_ISDF-Ex_t|>%.0e)=%d / %d samples, max|d|=%.6e', ...
    diff_warn_tol, n_mismatch, n_ob_sample, max_abs_mismatch);

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

  % Band summary as plain label/value cell columns (no table).
  band_items = {};
  band_vals = {};
  band_items{end + 1} = 'sum(E_HF)'; %#ok<AGROW>
  band_vals{end + 1} = sprintf('%.8e', sum_E_HF); %#ok<AGROW>
  band_items{end + 1} = 'sum(E_HF_ISDF)';
  band_vals{end + 1} = sprintf('%.8e', sum_E_ISDF);
  band_items{end + 1} = 'sum(E_HF) - sum(E_HF_ISDF)';
  band_vals{end + 1} = sprintf('%.8e', sum_E_HF - sum_E_ISDF);
  band_items{end + 1} = 'mean(E_HF) (over nb*nibz*nspin cells)';
  band_vals{end + 1} = sprintf('%.8e', mean_E_HF);
  band_items{end + 1} = 'mean(E_HF_ISDF) (over nb*nibz*nspin cells)';
  band_vals{end + 1} = sprintf('%.8e', mean_E_ISDF);
  band_items{end + 1} = 'sum(|E_HF - E_HF_ISDF|)';
  band_vals{end + 1} = sprintf('%.8e', sum_abs_diff_E);
  band_items{end + 1} = 'max(|E_HF - E_HF_ISDF|)';
  band_vals{end + 1} = sprintf('%.8e', max_abs_diff_E);
  band_items{end + 1} = '||vec(E_HF - E_HF_ISDF)||_2';
  band_vals{end + 1} = sprintf('%.8e', frob_diff);
  band_items{end + 1} = '||vec(E_HF - E_HF_ISDF)||_2 / ||vec(E_HF)||_2';
  band_vals{end + 1} = sprintf('%.8e', frob_diff / max(frob_hf, eps('double')));
  band_items{end + 1} = 'E_HF minimum';
  band_vals{end + 1} = sprintf('%.8e at (ib=%d, ik=%d, ispin=%d)', ...
    min_E_HF, ib_mn_hf, ik_mn_hf, is_mn_hf);
  band_items{end + 1} = 'E_HF maximum';
  band_vals{end + 1} = sprintf('%.8e at (ib=%d, ik=%d, ispin=%d)', ...
    max_E_HF, ib_mx_hf, ik_mx_hf, is_mx_hf);
  band_items{end + 1} = 'E_HF_ISDF minimum';
  band_vals{end + 1} = sprintf('%.8e at (ib=%d, ik=%d, ispin=%d)', ...
    min_E_ISDF, ib_mn_isdf, ik_mn_isdf, is_mn_isdf);
  band_items{end + 1} = 'E_HF_ISDF maximum';
  band_vals{end + 1} = sprintf('%.8e at (ib=%d, ik=%d, ispin=%d)', ...
    max_E_ISDF, ib_mx_isdf, ik_mx_isdf, is_mx_isdf);
  if nspin > 1
    for is = 1:double(nspin)
      s_hf = sum(sum(E_HF(:, :, is)));
      s_is = sum(sum(E_HF_ISDF(:, :, is)));
      s_ad = sum(sum(abs(E_HF(:, :, is) - E_HF_ISDF(:, :, is))));
      band_items{end + 1} = sprintf('ispin=%d: sum(E_HF)', is); %#ok<AGROW>
      band_vals{end + 1} = sprintf('%.8e', s_hf); %#ok<AGROW>
      band_items{end + 1} = sprintf('ispin=%d: sum(E_HF_ISDF)', is);
      band_vals{end + 1} = sprintf('%.8e', s_is);
      band_items{end + 1} = sprintf('ispin=%d: sum(|E_HF - E_HF_ISDF|)', is);
      band_vals{end + 1} = sprintf('%.8e', s_ad);
    end
  end

  n_rel_samples = st_RelSigned.n;

  stats_vars = {'Ex_t', 'Ex_ISDF', 'Diff', 'AbsDiff', 'Diff_over_abs_Ex_t', 'AbsDiff_over_abs_Ex_t'};
  streams = {st_Ex_t, st_Ex_ISDF, st_Diff, st_AbsDiff, st_RelSigned, st_RelAbs};
  stats_rows = {'Mean', 'Std', 'Var', 'Max'};
  stats_M = zeros(4, numel(stats_vars));
  for j = 1:numel(stats_vars)
    [stats_M(1, j), stats_M(2, j), stats_M(3, j), stats_M(4, j)] = ...
      isdf_validate_HF_stream_finalize(streams{j});
  end

  report = struct();
  report.isdf_id = id;
  report.desc = isdf_data.desc;
  report.interp_scheme = isdf_data.interp_scheme;
  report.nisdf = isdf_data.nisdf;
  report.band_items = band_items;
  report.band_vals = band_vals;
  report.stats_vars = stats_vars;
  report.stats_rows = stats_rows;
  report.stats_M = stats_M;
  report.n_samples = n_ob_sample;
  report.n_samples_rel = n_rel_samples;
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
  report.n_mismatch = int32(n_mismatch);
  report.max_abs_mismatch = double(max_abs_mismatch);

  fpath = isdf.report.hf(report);
  report.report_file = fpath;
  isdf.report.run_summary('hf', report);
  fprintf(1, 'isdf_validate_HF: wrote HF report to:\n  %s\n', fpath);

  if nargout <= 3
    report = [];
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

function isdf_ensure_timing_initialized()
  try
    timing.get();
  catch %#ok<CTCH>
    timing.driver();
  end
end
