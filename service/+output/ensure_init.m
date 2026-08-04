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
%   If a report file is already open (same path, or no new path given),
%   only refreshes verbosity. Otherwise forwards to output.init.

  p = inputParser;
  addParameter(p, 'report', '', @(x) ischar(x) || isstring(x));
  addParameter(p, 'log', '', @(x) ischar(x) || isstring(x));
  addParameter(p, 'verbose', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
  addParameter(p, 'screen', true, @(x) islogical(x) || isnumeric(x));
  parse(p, varargin{:});

  s = state_('get');
  report_path = char(string(p.Results.report));

  already = s.initialized && s.report_fid > 0;
  same_report = isempty(report_path) || strcmp(s.report_path, report_path);

  if already && same_report
    if ~isempty(p.Results.verbose)
      output.verbose(p.Results.verbose);
    end
    return
  end

  output.init(varargin{:});
end
