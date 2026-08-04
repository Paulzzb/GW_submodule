% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/01 ZZ

function set_logfile(path)
%SET_LOGFILE  Open a log file in exclusive legacy mode.
%
%   Screen-bound messages are redirected exclusively to this log file.
%   Prefer output.init('log', path) for new code.

  if nargin < 1 || ~(ischar(path) || isstring(path))
    error('output:set_logfile:BadPath', 'path must be a string.');
  end
  path = char(string(path));

  s = state_('get');
  if s.log_fid > 0
    try, fclose(s.log_fid); catch, end %#ok<CTCH>
  end

  d = fileparts(path);
  if ~isempty(d) && ~isfolder(d)
    mkdir(d);
  end

  [fid, errmsg] = fopen(path, 'a');
  if fid < 0
    error('output:set_logfile:OpenFailed', 'Cannot open logfile %s: %s', ...
      path, errmsg);
  end

  s.log_path = path;
  s.log_fid = fid;
  s.write_to_log = true;
  s.legacy_logfile_only = true;
  s.initialized = true;
  state_('set', s);
end
