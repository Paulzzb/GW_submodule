% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/01 ZZ

function out = state_(cmd, varargin)
%STATE_  Persistent state for the output package.
%
%   s = state_('get')
%   state_('set', s)
%   state_('reset')

  persistent S

  if isempty(S)
    S = local_default();
  end

  switch lower(cmd)
    case 'get'
      out = S;

    case 'set'
      if nargin < 2 || ~isstruct(varargin{1})
        error('output:state:BadSet', 'state_(''set'', s) requires a struct.');
      end
      S = varargin{1};
      out = S;

    case 'reset'
      local_close_all(S);
      S = local_default();
      out = S;

    otherwise
      error('output:state:UnknownCommand', 'Unknown command: %s', cmd);
  end
end

function s = local_default()
  s = struct( ...
    'initialized', false, ...
    'verbose', 1, ...
    'showtag', true, ...
    'module_stack', {{}}, ...
    'report_path', '', ...
    'report_fid', -1, ...
    'log_path', '', ...
    'log_fid', -1, ...
    'write_to_screen', true, ...
    'write_to_report', true, ...
    'write_to_log', true, ...
    'legacy_logfile_only', false, ...
    'depth', -1, ...
    'isec', zeros(1, 5), ...
    'sec_tics', {{}}, ...
    'named', struct('name', {}, 'path', {}, 'fid', {}) ...
    );
end

function local_close_all(s)
  if isstruct(s)
    if isfield(s, 'report_fid') && s.report_fid > 0
      try, fclose(s.report_fid); catch, end %#ok<CTCH>
    end
    if isfield(s, 'log_fid') && s.log_fid > 0
      try, fclose(s.log_fid); catch, end %#ok<CTCH>
    end
    if isfield(s, 'named')
      for i = 1:numel(s.named)
        if s.named(i).fid > 0
          try, fclose(s.named(i).fid); catch, end %#ok<CTCH>
        end
      end
    end
  end
end
