function check(varargin)
%CHECK  If on(tag), evaluate fh(); react when fh() is true (violation).
%
%   isdftest.debug.check(tag, fh, msg [, tagForId])
%   isdftest.debug.check(config, tag, fh, msg [, tagForId])

  if nargin >= 4 && isa(varargin{1}, 'struct') && isfield(varargin{1}, 'ISDF')
    config = varargin{1};
    tag = varargin{2};
    fh = varargin{3};
    msg = varargin{4};
    if ~isdftest.debug.on(config, tag)
      return
    end
    if ~isa(fh, 'function_handle')
      error('isdf:debug:checkArgs', 'fh must be a function_handle.');
    end
    bad = logical(fh());
    if nargin >= 5
      isdftest.debug.react(config, bad, msg, varargin{5});
    else
      isdftest.debug.react(config, bad, msg, tag);
    end
    return
  end

  if nargin < 3
    error('isdf:debug:checkArgs', ...
      'Use check(tag,fh,msg) or check(config,tag,fh,msg).');
  end
  tag = varargin{1};
  fh = varargin{2};
  msg = varargin{3};
  if ~isdftest.debug.on(tag)
    return
  end
  if ~isa(fh, 'function_handle')
    error('isdf:debug:checkArgs', 'fh must be a function_handle.');
  end
  bad = logical(fh());
  if nargin >= 4
    isdftest.debug.react(bad, msg, varargin{4});
  else
    isdftest.debug.react(bad, msg, tag);
  end
end
