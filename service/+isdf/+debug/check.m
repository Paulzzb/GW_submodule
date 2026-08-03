% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/06 ZZ

function check(varargin)
%CHECK  If on(tag), evaluate fh(); react when fh() is true (violation).
%
%   isdf.debug.check(tag, fh, msg [, tagForId])
%   isdf.debug.check(config, tag, fh, msg [, tagForId])

  if nargin >= 4 && isa(varargin{1}, 'struct') && isfield(varargin{1}, 'ISDF')
    config = varargin{1};
    tag = varargin{2};
    fh = varargin{3};
    msg = varargin{4};
    if ~isdf.debug.on(config, tag)
      return
    end
    if ~isa(fh, 'function_handle')
      error('isdf:debug:checkArgs', 'fh must be a function_handle.');
    end
    bad = logical(fh());
    if nargin >= 5
      isdf.debug.react(config, bad, msg, varargin{5});
    else
      isdf.debug.react(config, bad, msg, tag);
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
  if ~isdf.debug.on(tag)
    return
  end
  if ~isa(fh, 'function_handle')
    error('isdf:debug:checkArgs', 'fh must be a function_handle.');
  end
  bad = logical(fh());
  if nargin >= 4
    isdf.debug.react(bad, msg, varargin{4});
  else
    isdf.debug.react(bad, msg, tag);
  end
end
