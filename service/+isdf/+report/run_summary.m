% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/07 ZZ

function run_summary(cmd, arg)
%RUN_SUMMARY  Accumulate and write a short ISDF run summary to r-*.
%
%   isdf.report.run_summary('clear')
%   isdf.report.run_summary('seed', struct('desc',..,'seed_id',..,'seed_nmu',..,'method',..))
%   isdf.report.run_summary('adaptive', adaptive_report)   % from adaptiveisdf
%   isdf.report.run_summary('hf', hf_report)               % from isdf_validate_HF
%   isdf.report.run_summary('write')
%
% Detail OF paths point under isdf_report/.

  persistent bag order

  if nargin < 1 || isempty(cmd)
    output.err('run_summary requires a command.');
  end
  cmd = lower(strtrim(char(string(cmd))));

  switch cmd
    case 'clear'
      bag = struct();
      order = {};

    case 'seed'
      key = local_key(arg.desc);
      e = local_get(bag, key);
      e.seed_id = int32(arg.seed_id);
      e.seed_nmu = int32(arg.seed_nmu);
      e.method = char(string(arg.method));
      bag.(key) = e;
      order = local_push_order(order, key);

    case 'adaptive'
      r = arg;
      key = local_key(r.isdf_desc);
      e = local_get(bag, key);
      e.seed_id = int32(r.coarse_isdf_id);
      if isfield(r, 'adaptive_isdf_id')
        e.id_new = int32(r.adaptive_isdf_id);
      end
      e.seed_nmu = int32(r.initial_nisdf);
      e.final_nmu = int32(r.final_nisdf);
      e.added = int32(r.added_nisdf);
      e.rel = double(r.relative_loss);
      e.stop = char(string(r.stop_reason));
      e.iters = int32(r.iterations);
      e.wall = double(r.elapsed_phase1_seconds);
      e.adaptive_file = sprintf(filename_map().adaptive_report, double(r.coarse_isdf_id));
      bag.(key) = e;
      order = local_push_order(order, key);

    case 'hf'
      r = arg;
      desc = 'unknown';
      if isfield(r, 'desc') && strlength(string(r.desc)) > 0
        desc = char(string(r.desc));
      end
      key = local_key(desc);
      e = local_get(bag, key);
      if isfield(r, 'isdf_id')
        e.id_new = int32(r.isdf_id);
      end
      if isfield(r, 'n_mismatch')
        e.hf_mismatch = int32(r.n_mismatch);
      end
      if isfield(r, 'n_samples')
        e.hf_nsamp = int32(r.n_samples);
      end
      if isfield(r, 'max_abs_mismatch')
        e.hf_maxd = double(r.max_abs_mismatch);
      elseif isfield(r, 'max_abs_diff_E')
        e.hf_maxd = double(r.max_abs_diff_E);
      end
      if isfield(r, 'report_file') && ~isempty(r.report_file)
        [~, name, ext] = fileparts(char(string(r.report_file)));
        e.hf_file = [name, ext];
      else
        e.hf_file = sprintf(filename_map().hf_report, double(r.isdf_id));
      end
      bag.(key) = e;
      order = local_push_order(order, key);

    case 'write'
      output.msg('nrs', '----------- ISDF run summary -----------');
      if isempty(order)
        output.msg('rs', ' (no ISDF slots recorded)');
        return;
      end
      for i = 1:numel(order)
        local_emit_entry(bag.(order{i}));
      end
      def = filename_map();
      output.msg('r', '');
      output.msg('rs', ' Cond diagnostics: %s/%s', def.isdf_report_dir, def.cond_report);

    otherwise
      output.err('Unknown run_summary command ''%s''.', cmd);
  end
end

function key = local_key(desc)
  key = lower(strtrim(char(string(desc))));
  if isempty(key)
    key = 'unknown';
  end
  % struct field names cannot start with digit; vc/vn/nn are fine
  if ~isvarname(key)
    key = matlab.lang.makeValidName(key);
  end
end

function e = local_get(bag, key)
  if isstruct(bag) && isfield(bag, key)
    e = bag.(key);
    return;
  end
  e = struct();
  e.desc = key;
  e.method = '';
  e.seed_id = int32(-1);
  e.id_new = int32(-1);
  e.seed_nmu = int32(-1);
  e.final_nmu = int32(-1);
  e.added = int32(-1);
  e.rel = nan;
  e.stop = '';
  e.iters = int32(-1);
  e.wall = nan;
  e.adaptive_file = '';
  e.hf_file = '';
  e.hf_mismatch = int32(-1);
  e.hf_nsamp = int32(-1);
  e.hf_maxd = nan;
end

function order = local_push_order(order, key)
  if isempty(order)
    order = {key};
    return;
  end
  if ~any(strcmp(order, key))
    order{end + 1} = key; 
  end
end

function local_emit_entry(e)
  method = e.method;
  if isempty(method)
    method = '?';
  end
  seed_id = e.seed_id;
  id_new = e.id_new;
  n0 = e.seed_nmu;
  n1 = e.final_nmu;
  if n1 < 0 && n0 >= 0
    n1 = n0;
  end
  dN = int32(0);
  if n0 >= 0 && n1 >= 0
    dN = n1 - n0;
  elseif e.added >= 0
    dN = e.added;
  end

  output.msg('rs', ' %s  method=%-7s  seed id=%d  Nmu=%d  ->  id=%d  Nmu=%d  (%+d)', ...
    e.desc, method, seed_id, n0, id_new, n1, dN);

  if ~isnan(e.rel) || strlength(string(e.stop)) > 0
    if e.iters >= 0 && ~isnan(e.wall)
      output.msg('r', '     adaptive: rel=%.2e  stop=%s  iters=%d  wall=%.1f s', ...
        e.rel, e.stop, e.iters, e.wall);
    elseif ~isnan(e.rel)
      output.msg('r', '     adaptive: rel=%.2e  stop=%s', e.rel, e.stop);
    end
  end
  if e.hf_mismatch >= 0
    output.msg('r', '     HF: max|dEx|=%.2e  mismatches=%d/%d', ...
      e.hf_maxd, e.hf_mismatch, e.hf_nsamp);
  end

  def = filename_map();
  details = {};
  if ~isempty(e.adaptive_file)
    details{end + 1} = fullfile(def.isdf_report_dir, e.adaptive_file); %#ok<AGROW>
  end
  if ~isempty(e.hf_file)
    details{end + 1} = fullfile(def.isdf_report_dir, e.hf_file); %#ok<AGROW>
  end
  if ~isempty(details)
    output.msg('r', '     details: %s', strjoin(details, '  '));
  end
end
