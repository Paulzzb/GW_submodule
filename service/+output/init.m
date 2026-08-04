% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/01 ZZ

function init(varargin)
%INIT  Initialize the output package (report/log files, verbosity).
%
%   output.init()
%   output.init('report', path, 'log', path, 'verbose', n)
%
%   Name-value pairs:
%     report   - path to report file (created/overwritten after rename)
%     log      - path to log file
%     verbose  - 0 quiet / 1 normal / 2 debug (default unchanged if omitted)
%     screen   - logical, write to stdout (default true)

  p = inputParser;
  addParameter(p, 'report', '', @(x) ischar(x) || isstring(x));
  addParameter(p, 'log', '', @(x) ischar(x) || isstring(x));
  addParameter(p, 'verbose', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
  addParameter(p, 'screen', true, @(x) islogical(x) || isnumeric(x));
  parse(p, varargin{:});

  s = state_('get');

  % Close previous destinations but keep verbosity/tags unless reset via free.
  if s.report_fid > 0
    try, fclose(s.report_fid); catch, end %#ok<CTCH>
  end
  if s.log_fid > 0
    try, fclose(s.log_fid); catch, end %#ok<CTCH>
  end
  for i = 1:numel(s.named)
    if s.named(i).fid > 0
      try, fclose(s.named(i).fid); catch, end %#ok<CTCH>
    end
  end
  s.named = struct('name', {}, 'path', {}, 'fid', {});
  s.report_fid = -1;
  s.log_fid = -1;
  s.report_path = '';
  s.log_path = '';
  s.legacy_logfile_only = false;

  if ~isempty(p.Results.verbose)
    s.verbose = max(0, min(2, round(double(p.Results.verbose))));
  end
  s.write_to_screen = logical(p.Results.screen);

  report_path = char(string(p.Results.report));
  log_path = char(string(p.Results.log));

  if ~isempty(report_path)
    local_ensure_parent(report_path);
    local_rename_if_exists(report_path);
    [fid, errmsg] = fopen(report_path, 'w');
    if fid < 0
      error('output:init:OpenFailed', 'Cannot open report file %s: %s', ...
        report_path, errmsg);
    end
    s.report_path = report_path;
    s.report_fid = fid;
    s.write_to_report = true;
  else
    s.write_to_report = false;
  end

  if ~isempty(log_path)
    local_ensure_parent(log_path);
    local_rename_if_exists(log_path);
    [fid, errmsg] = fopen(log_path, 'w');
    if fid < 0
      error('output:init:OpenFailed', 'Cannot open log file %s: %s', ...
        log_path, errmsg);
    end
    s.log_path = log_path;
    s.log_fid = fid;
    s.write_to_log = true;
  else
    s.write_to_log = false;
  end

  s.initialized = true;
  s.depth = -1;
  s.isec = zeros(1, 5);
  s.sec_tics = {};
  state_('set', s);
end

function local_ensure_parent(path)
  d = fileparts(path);
  if ~isempty(d) && ~isfolder(d)
    mkdir(d);
  end
end

function local_rename_if_exists(path)
  if ~isfile(path)
    return
  end
  [folder, name, ext] = fileparts(path);
  for k = 1:99
    cand = fullfile(folder, sprintf('%s_%02d%s', name, k, ext));
    if ~isfile(cand)
      movefile(path, cand);
      return
    end
  end
  warning('output:init:RenameFailed', ...
    'Could not rotate existing file %s; overwriting.', path);
end
