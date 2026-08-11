% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/01 ZZ

function ensure_init(varargin)
%ENSURE_INIT  Initialize output if needed; no-op when already open.
%
%   output.ensure_init('report', path, 'verbose', n)
%
%   If a report file is already open at the same absolute path and the
%   file still exists on disk, only refreshes verbosity. Otherwise
%   forwards to output.init (reopens after free/cd/clean).

  p = inputParser;
  addParameter(p, 'report', '', @(x) ischar(x) || isstring(x));
  addParameter(p, 'log', '', @(x) ischar(x) || isstring(x));
  addParameter(p, 'verbose', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
  addParameter(p, 'screen', true, @(x) islogical(x) || isnumeric(x));
  parse(p, varargin{:});

  s = state_('get');
  report_path = abspath_(p.Results.report);
  log_path = abspath_(p.Results.log);

  already = s.initialized && s.report_fid > 0;
  % File may have been deleted (clean_case_outputs) while the FID was still
  % open — on Linux that leaves a ghost inode; reopen in that case.
  alive = isempty(s.report_path) || isfile(s.report_path);
  same_report = isempty(report_path) || strcmp(s.report_path, report_path);

  if already && alive && same_report
    if ~isempty(p.Results.verbose)
      output.verbose(p.Results.verbose);
    end
    return
  end

  args = {};
  if ~isempty(report_path)
    args = [args, {'report', report_path}]; %#ok<AGROW>
  end
  if ~isempty(log_path)
    args = [args, {'log', log_path}]; %#ok<AGROW>
  end
  if ~isempty(p.Results.verbose)
    args = [args, {'verbose', p.Results.verbose}]; %#ok<AGROW>
  end
  args = [args, {'screen', logical(p.Results.screen)}];
  output.init(args{:});
end
