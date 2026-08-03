% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/06/11 ZZ

function numerical_cond_report(action, varargin)
%NUMERICAL_COND_REPORT  Append numeric diagnostics to numerical_cond_report.txt (pwd).
%
%   isdf.numerical_cond_report('init')
%   isdf.numerical_cond_report('adaptive', coarse_id, idnew, desc, loss_before_step1, loss_after, final_nisdf)
%   isdf.numerical_cond_report('gen_tildevq_begin', id, desc, s_cut)
%   isdf.numerical_cond_report('gen_tildevq_iq', id, iqibz, CCHq, MCHq, l_keep, s_cut)

  action = lower(strtrim(char(string(action))));

  switch action
    case 'init'
      local_write_init();

    case 'adaptive'
      if numel(varargin) < 6
        error('numerical_cond_report:adaptive', ...
          'Requires coarse_id, idnew, desc, loss_before_step1, loss_after, final_nisdf.');
      end
      local_write_adaptive(varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5}, varargin{6});

    case 'gen_tildevq_begin'
      if numel(varargin) < 3
        error('numerical_cond_report:gen_tildevq_begin', 'Requires id, desc, s_cut.');
      end
      local_write_gen_tildeVq_begin(varargin{1}, varargin{2}, varargin{3});

    case 'gen_tildevq_iq'
      if numel(varargin) < 6
        error('numerical_cond_report:gen_tildevq_iq', ...
          'Requires id, iqibz, CCHq, MCHq, l_keep, s_cut.');
      end
      local_write_gen_tildeVq_iq(varargin{1}, varargin{2}, varargin{3}, ...
        varargin{4}, varargin{5}, varargin{6});

    otherwise
      error('numerical_cond_report:action', 'Unknown action ''%s''.', action);
  end
end

function fpath = local_path()
  fpath = fullfile(pwd, 'numerical_cond_report.txt');
end

function s = local_fe(x)
  s = sprintf('%.4e', double(x));
end

function local_append(lines)
  fpath = local_path();
  fid = fopen(fpath, 'a');
  if fid < 0
    error('numerical_cond_report:io', 'Cannot open %s for append.', fpath);
  end
  oc = onCleanup(@() fclose(fid)); %#ok<NASGU>
  for k = 1:numel(lines)
    fprintf(fid, '%s\n', lines{k});
  end
end

function local_write_init()
  fpath = local_path();
  fid = fopen(fpath, 'w');
  if fid < 0
    error('numerical_cond_report:io', 'Cannot open %s for write.', fpath);
  end
  oc = onCleanup(@() fclose(fid)); %#ok<NASGU>
  fprintf(fid, '=== isdf numerical_cond_report ===\n');
  fprintf(fid, 'Generated: %s\n', datestr(now, 31));
  fprintf(fid, 'pwd: %s\n\n', pwd);
end

function local_write_adaptive(coarse_id, idnew, desc, loss_before_step1, loss_after, final_nisdf)
  lines = {
    '--- adaptiveisdf ---'
    sprintf('coarse_id: %d', round(double(coarse_id)))
    sprintf('idnew: %d', round(double(idnew)))
    sprintf('desc: %s', char(string(desc)))
    sprintf('final_nisdf: %d', round(double(final_nisdf)))
    sprintf('loss_before_step1: %s', local_fe(loss_before_step1))
    sprintf('loss_after_complete: %s', local_fe(loss_after))
    ''
  };
  local_append(lines);
end

function local_write_gen_tildeVq_begin(id, desc, s_cut)
  lines = {
    '--- gen_tildeVq ---'
    sprintf('id: %d', round(double(id)))
    sprintf('desc: %s', char(string(desc)))
    sprintf('s_cut (inv_param): %s', local_fe(s_cut))
    ''
  };
  local_append(lines);
end

function stats = local_cchq_stats(CCHq)
  C_sym = (CCHq + CCHq') / 2;
  sigma = abs(real(eig(C_sym)));
  sigma = sort(sigma, 'descend');
  fro_norm = norm(C_sym, 'fro');
  tol = max(eps, fro_norm * eps);
  sigma_pos = sigma(sigma > tol);
  if isempty(sigma_pos)
    stats.sigma_max = NaN;
    stats.sigma_min = NaN;
    stats.cond_number = NaN;
  else
    stats.sigma_max = sigma_pos(1);
    stats.sigma_min = sigma_pos(end);
    stats.cond_number = stats.sigma_max / stats.sigma_min;
  end
  stats.fro_norm = fro_norm;
end

function local_write_gen_tildeVq_iq(id, iqibz, CCHq, MCHq, l_keep, s_cut)
  st = local_cchq_stats(CCHq);
  sum_trunc_sigma = local_truncated_sigma_sum(CCHq, s_cut);
  fro_M = norm(MCHq, 'fro');
  lines = {
    sprintf('  iqibz: %d', round(double(iqibz)))
    sprintf('  CCHq cond_number: %s', local_fe(st.cond_number))
    sprintf('  CCHq sigma_max: %s', local_fe(st.sigma_max))
    sprintf('  CCHq sigma_min: %s', local_fe(st.sigma_min))
    sprintf('  CCHq fro_norm: %s', local_fe(st.fro_norm))
    sprintf('  CCHq l_keep (SVD truncate): %d', round(double(l_keep)))
    sprintf('  CCHq sum_truncated_sigma: %s', local_fe(sum_trunc_sigma))
    sprintf('  CCHq s_cut: %s', local_fe(s_cut))
    sprintf('  MCHq fro_norm: %s', local_fe(fro_M))
    ''
  };
  local_append(lines);
end

function sum_trunc = local_truncated_sigma_sum(CCHq, s_cut)
  CCHq = (CCHq + CCHq') / 2;
  sigma = abs(eig(CCHq, 'vector'));
  sigma = sort(sigma, 'descend');
  s1 = sigma(1);
  keep = find(sigma > s1 * s_cut);
  drop = setdiff(1:numel(sigma), keep);
  sum_trunc = sum(sigma(drop));
end
