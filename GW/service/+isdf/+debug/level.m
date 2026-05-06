function lev = level(varargin)
%LEVEL  Normalized reaction level ('error' or 'warn') for isdf.debug.react.
%
%   lev = isdf.debug.level()  — from session cache (init_from_config).
%   lev = isdf.debug.level(config)  — from config.ISDF only (no cache read).

  if nargin >= 1 && isa(varargin{1}, 'struct') && isfield(varargin{1}, 'ISDF')
    lev = local_level_from_isdf(varargin{1}.ISDF);
    return
  end

  st = isdf.debug.cache('get');
  if st.initialized
    lev = char(string(st.level));
    if ~strcmp(lev, 'warn') && ~strcmp(lev, 'error')
      lev = 'error';
    end
  else
    lev = 'error';
  end
end

function lev = local_level_from_isdf(s)
  lev = 'error';
  if ~isfield(s, 'debug_level') || isempty(s.debug_level)
    return
  end
  v = lower(strtrim(char(string(s.debug_level))));
  if strcmp(v, 'warn') || strcmp(v, 'warning')
    lev = 'warn';
  elseif strcmp(v, 'error')
    lev = 'error';
  else
    warning('isdf:debug:badLevel', ...
      'Unknown config.ISDF.debug_level ''%s''; using ''error''.', char(string(s.debug_level)));
  end
end
