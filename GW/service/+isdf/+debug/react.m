function react(varargin)
%REACT  Apply debug_level when a check found a violation.
%
%   isdf.debug.react(violation, msg)
%   isdf.debug.react(violation, msg, tagForId)
%       Uses session cache from init_from_config for level().
%
%   isdf.debug.react(config, violation, msg [, tagForId])   % optional overload
%       Uses config.ISDF.debug_level only.

  if nargin >= 3 && isa(varargin{1}, 'struct') && isfield(varargin{1}, 'ISDF')
    config = varargin{1};
    violation = varargin{2};
    msg = varargin{3};
    tagForId = '';
    if nargin >= 4
      tagForId = varargin{4};
    end
    lev = isdf.debug.level(config);
  elseif nargin >= 2
    violation = varargin{1};
    msg = varargin{2};
    tagForId = '';
    if nargin >= 3
      tagForId = varargin{3};
    end
    lev = isdf.debug.level();
  else
    error('isdf:debug:reactArgs', 'Use react(violation, msg) or react(config, violation, msg).');
  end

  if isempty(violation) || ~logical(violation(1))
    return
  end
  msg = char(string(msg));
  wid = 'isdf:debug:check';
  if ~isempty(tagForId)
    t = char(string(tagForId));
    t = strtrim(t);
    if strlength(t) > 0
      t = regexprep(t, '[^a-zA-Z0-9_]', '_', 'all');
      wid = ['isdf:debug:' t];
    end
  end

  switch lev
    case 'warn'
      warning(wid, '%s', msg);
    otherwise
      error(wid, '%s', msg);
  end
end
